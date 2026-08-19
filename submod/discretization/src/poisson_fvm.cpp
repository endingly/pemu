#include <pemu/discretization/poisson_fvm.hpp>

#include <Eigen/SparseCore>

#include <cmath>
#include <concepts>
#include <stdexcept>
#include <type_traits>
#include <variant>
#include <vector>

namespace pemu::discretization {

namespace {

double dot(const mesh::Vec3& a, const mesh::Vec3& b) noexcept {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

}  // namespace

PoissonFvm::PoissonFvm(
    const mesh::IMesh& mesh, const field::CellField<double>& source,
    const double epsilon,
    const boundary::BoundaryConditionSet& boundary_conditions)
    : mesh_(&mesh),
      source_(&source),
      boundary_conditions_(&boundary_conditions),
      epsilon_(epsilon) {
  if (&source.mesh() != &mesh) {
    throw std::invalid_argument(
        "PoissonFvm source field belongs to a different mesh");
  }

  if (epsilon <= 0.0) {
    throw std::invalid_argument("PoissonFvm epsilon must be positive");
  }

  if (mesh.dimension() != 2) {
    throw std::invalid_argument("PoissonFvm currently supports only 2D meshes");
  }
}

// ============================================================
// Internal face coefficient
//
//           epsilon A_f
//     c_f = -----------
//               d_PN
//
// We use the normal projection:
//
//     d_PN = (x_N - x_P) dot n_f
//
// because mesh guarantees:
//
//     n_f : owner -> neighbor
//
// so for a valid FVM face:
//
//     d_PN > 0
//
// This is currently a TPFA discretization.
// ============================================================
double PoissonFvm::internalFaceCoefficient(const mesh::FaceId face) const {
  const auto owner = mesh_->owner(face);

  const auto neighbor = mesh_->neighbor(face);

  if (neighbor == mesh::invalid_cell) {
    throw std::invalid_argument(
        "internalFaceCoefficient called on boundary face");
  }

  const auto xp = mesh_->cellCenter(owner);

  const auto xn = mesh_->cellCenter(neighbor);

  const auto normal = mesh_->faceNormal(face);

  const auto delta = xn - xp;

  const double distance = dot(delta, normal);

  if (distance <= 0.0) {
    throw std::runtime_error("invalid owner-neighbor normal distance");
  }

  return epsilon_ * mesh_->faceArea(face) / distance;
}

// ============================================================
// Dirichlet boundary coefficient
//
// Boundary face:
//
//          epsilon A_f
//     c = -------------
//             d_Pb
//
// where:
//
//     d_Pb =
//       (x_face - x_owner) dot n_f
//
// Boundary flux:
//
//     q = -epsilon grad(phi) dot n
//
// approximated by:
//
//     q A
//
//       = epsilon A / d * (phi_P - phi_b)
//
//       = c phi_P - c phi_b
//
// Therefore:
//
//     A_PP += c
//     b_P  += c phi_b
// ============================================================
double PoissonFvm::boundaryFaceCoefficient(const mesh::FaceId face) const {
  if (!mesh_->isBoundary(face)) {
    throw std::invalid_argument(
        "boundaryFaceCoefficient called on internal face");
  }

  const auto owner = mesh_->owner(face);

  const auto xp = mesh_->cellCenter(owner);

  const auto xf = mesh_->faceCenter(face);

  const auto normal = mesh_->faceNormal(face);

  const auto delta = xf - xp;

  const double distance = dot(delta, normal);

  if (distance <= 0.0) {
    throw std::runtime_error("invalid cell-to-boundary normal distance");
  }

  return epsilon_ * mesh_->faceArea(face) / distance;
}

// ============================================================
// Assembly
//
// Equation:
//
//     -div(epsilon grad(phi)) = rho
//
// Integral over cell P:
//
//     sum_faces q_f A_f = rho_P V_P
//
// where:
//
//     q_f = -epsilon grad(phi) dot n
//
// Internal face:
//
//     q_f A_f
//       ≈ c_f (phi_P - phi_N)
//
// therefore pair contribution:
//
//        [ +c  -c ]
//        [ -c  +c ]
//
// Dirichlet:
//
//     A_PP += c
//     b_P  += c phi_b
//
// Neumann:
//
//     q is already known.
//
//     unknown_flux_sum + q A = rho V
//
// therefore:
//
//     b_P -= q A
// ============================================================

linalg::LinearSystem PoissonFvm::assemble() const {
  using Triplet = Eigen::Triplet<linalg::Scalar, linalg::Index>;

  const auto num_cells = mesh_->numCells();

  linalg::LinearSystem system{
      .A = linalg::SparseMatrix(static_cast<linalg::Index>(num_cells),
                                static_cast<linalg::Index>(num_cells)),
      .b = linalg::Vector::Zero(static_cast<linalg::Index>(num_cells))};

  std::vector<Triplet> entries;

  //
  // Rough estimate:
  //
  // each internal face contributes four triplets,
  // boundary Dirichlet contributes one.
  //
  entries.reserve(mesh_->numFaces() * 4);

  // --------------------------------------------------------
  // Volume source term
  //
  //     b_P = rho_P V_P
  // --------------------------------------------------------

  for (mesh::CellId cell = 0; cell < mesh_->numCells(); ++cell) {

    const auto row = static_cast<linalg::Index>(cell);

    system.b[row] = (*source_)[cell] * mesh_->cellVolume(cell);
  }

  // --------------------------------------------------------
  // Face contributions
  // --------------------------------------------------------

  for (mesh::FaceId face = 0; face < mesh_->numFaces(); ++face) {

    const auto owner = mesh_->owner(face);

    const auto p = static_cast<linalg::Index>(owner);

    // ====================================================
    // Internal face
    // ====================================================

    if (!mesh_->isBoundary(face)) {

      const auto neighbor = mesh_->neighbor(face);

      const auto n = static_cast<linalg::Index>(neighbor);

      const double coefficient = internalFaceCoefficient(face);

      //
      // Owner equation
      //
      entries.emplace_back(p, p, coefficient);

      entries.emplace_back(p, n, -coefficient);

      //
      // Neighbor equation
      //
      entries.emplace_back(n, n, coefficient);

      entries.emplace_back(n, p, -coefficient);

      continue;
    }

    // ====================================================
    // Boundary face
    // ====================================================

    const auto boundary_id = mesh_->boundaryId(face);

    if (boundary_id == mesh::invalid_boundary) {

      throw std::runtime_error("boundary face has no BoundaryId");
    }

    if (!boundary_conditions_->contains(boundary_id)) {

      throw std::runtime_error("boundary face has no boundary condition");
    }

    const auto& condition = boundary_conditions_->at(boundary_id);

    std::visit(
        [&](const auto& bc) {
          using BcType = std::remove_cvref_t<decltype(bc)>;

          // --------------------------------------------
          // Dirichlet
          // --------------------------------------------

          if constexpr (std::same_as<BcType, boundary::Dirichlet>) {

            const double coefficient = boundaryFaceCoefficient(face);

            entries.emplace_back(p, p, coefficient);

            system.b[p] += coefficient * bc.value;
          }

          // --------------------------------------------
          // Neumann
          //
          // Convention:
          //
          //   value =
          //     -epsilon grad(phi) dot n
          //
          // therefore:
          //
          //   b -= q A
          // --------------------------------------------

          else if constexpr (std::same_as<BcType, boundary::Neumann>) {

            system.b[p] -= bc.value * mesh_->faceArea(face);
          }
        },
        condition);
  }

  system.A.setFromTriplets(entries.begin(), entries.end());

  system.A.makeCompressed();

  return system;
}

}  // namespace pemu::discretization