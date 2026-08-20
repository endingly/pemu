#include <Eigen/SparseCore>

#include <pemu/discretization/poisson_fvm.hpp>
#include <pemu/mesh/geometry.hpp>

#include <cmath>
#include <concepts>
#include <stdexcept>
#include <type_traits>
#include <variant>
#include <vector>

namespace pemu::discretization {

void PoissonFvm::assembleRhs(linalg::Vector& b) const {
  diffusion_.assembleBoundaryRhs(b);

  for (mesh::CellId cell = 0; cell < mesh_->numCells(); ++cell) {

    b[static_cast<linalg::Index>(cell)] +=
        (*source_)[cell] * mesh_->cellVolume(cell);
  }
}

void PoissonFvm::assembleMatrix(linalg::SparseMatrix& A) const {
  diffusion_.assembleMatrix(A);
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
  linalg::LinearSystem system;
  assembleMatrix(system.A);
  assembleRhs(system.b);
  return system;
}

bool PoissonFvm::hasDirichletBoundary() const {
  for (mesh::FaceId face = 0; face < mesh_->numFaces(); ++face) {
    if (!mesh_->isBoundary(face)) {
      continue;
    }

    const auto boundary_id = mesh_->boundaryId(face);
    if (boundary_id == mesh::invalid_boundary) {
      throw std::runtime_error("boundary face has no BoundaryId");
    }

    if (std::holds_alternative<boundary::Dirichlet>(bc_->at(boundary_id))) {
      return true;
    }
  }

  return false;
}

}  // namespace pemu::discretization
