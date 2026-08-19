#pragma once

#include <pemu/boundary/boundary_condition_set.hpp>
#include <pemu/linalg/types.hpp>
#include <pemu/mesh/geometry.hpp>
#include <pemu/mesh/i_mesh.hpp>

#include <Eigen/SparseCore>

#include <concepts>
#include <stdexcept>
#include <type_traits>
#include <variant>
#include <vector>

namespace pemu::discretization::operators {

class ScalarDiffusionFvm {
 public:
  ScalarDiffusionFvm(const mesh::IMesh& mesh, double coefficient,
                     const boundary::BoundaryConditionSet& bc)
      : mesh_(&mesh), coefficient_(coefficient), bc_(&bc) {
    if (coefficient < 0.0) {
      throw std::invalid_argument("diffusion coefficient must be non-negative");
    }
  }

  [[nodiscard]]
  const mesh::IMesh& mesh() const noexcept {
    return *mesh_;
  }

  void assembleMatrix(linalg::SparseMatrix& K) const {
    using Triplet = Eigen::Triplet<linalg::Scalar, linalg::Index>;

    const auto n = static_cast<linalg::Index>(mesh_->numCells());

    K.resize(n, n);

    std::vector<Triplet> entries;
    entries.reserve(mesh_->numFaces() * 4);

    for (mesh::FaceId face = 0; face < mesh_->numFaces(); ++face) {

      const auto owner = mesh_->owner(face);

      const auto p = static_cast<linalg::Index>(owner);

      // --------------------------------------------
      // Internal face
      // --------------------------------------------

      if (!mesh_->isBoundary(face)) {

        const auto neighbor = mesh_->neighbor(face);

        const auto q = static_cast<linalg::Index>(neighbor);

        const double c = internalCoefficient(face);

        entries.emplace_back(p, p, +c);

        entries.emplace_back(p, q, -c);

        entries.emplace_back(q, q, +c);

        entries.emplace_back(q, p, -c);

        continue;
      }

      // --------------------------------------------
      // Boundary face
      // --------------------------------------------

      const auto boundary_id = mesh_->boundaryId(face);

      if (boundary_id == mesh::invalid_boundary) {

        throw std::runtime_error("boundary face has no BoundaryId");
      }

      const auto& condition = bc_->at(boundary_id);

      std::visit(
          [&](const auto& value) {
            using BcType = std::remove_cvref_t<decltype(value)>;

            if constexpr (std::same_as<BcType, boundary::Dirichlet>) {

              entries.emplace_back(p, p, boundaryCoefficient(face));
            }

            // Neumann adds no matrix entry.
          },
          condition);
    }

    K.setFromTriplets(entries.begin(), entries.end());

    K.makeCompressed();
  }

  void assembleBoundaryRhs(linalg::Vector& rhs) const {
    rhs = linalg::Vector::Zero(static_cast<linalg::Index>(mesh_->numCells()));

    for (mesh::FaceId face = 0; face < mesh_->numFaces(); ++face) {

      if (!mesh_->isBoundary(face)) {
        continue;
      }

      const auto owner = mesh_->owner(face);

      const auto p = static_cast<linalg::Index>(owner);

      const auto boundary_id = mesh_->boundaryId(face);

      if (boundary_id == mesh::invalid_boundary) {

        throw std::runtime_error("boundary face has no BoundaryId");
      }

      const auto& condition = bc_->at(boundary_id);

      std::visit(
          [&](const auto& value) {
            using BcType = std::remove_cvref_t<decltype(value)>;

            // K u = ... + c u_b
            if constexpr (std::same_as<BcType, boundary::Dirichlet>) {

              rhs[p] += boundaryCoefficient(face) * value.value;
            }

            // q = -k grad(u)·n
            //
            // unknown terms + q A = source
            //
            // => rhs -= q A
            else if constexpr (std::same_as<BcType, boundary::Neumann>) {

              rhs[p] -= value.value * mesh_->faceArea(face);
            }
          },
          condition);
    }
  }

 private:
  [[nodiscard]]
  double internalCoefficient(mesh::FaceId face) const {
    const auto p = mesh_->owner(face);

    const auto n = mesh_->neighbor(face);

    const auto delta = mesh_->cellCenter(n) - mesh_->cellCenter(p);

    const auto normal = mesh_->faceNormal(face);

    const double distance = mesh::dot(delta, normal);

    if (distance <= 0.0) {
      throw std::runtime_error("invalid owner-neighbor normal distance");
    }

    return coefficient_ * mesh_->faceArea(face) / distance;
  }

  [[nodiscard]]
  double boundaryCoefficient(mesh::FaceId face) const {
    const auto p = mesh_->owner(face);

    const auto delta = mesh_->faceCenter(face) - mesh_->cellCenter(p);

    const double distance = mesh::dot(delta, mesh_->faceNormal(face));

    if (distance <= 0.0) {
      throw std::runtime_error("invalid cell-boundary normal distance");
    }

    return coefficient_ * mesh_->faceArea(face) / distance;
  }

 private:
  const mesh::IMesh* mesh_;
  double coefficient_;
  const boundary::BoundaryConditionSet* bc_;
};

}  // namespace pemu::discretization