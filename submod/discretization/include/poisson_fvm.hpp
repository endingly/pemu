#pragma once

#include <pemu/boundary/boundary_condition_set.hpp>
#include <pemu/discretization/operators/scalar_diffusion_fvm.hpp>
#include <pemu/field/cell_field.hpp>
#include <pemu/linalg/types.hpp>
#include <pemu/mesh/i_mesh.hpp>

namespace pemu::discretization {

class PoissonFvm {
 public:
  PoissonFvm(const mesh::IMesh& mesh, const field::CellField<double>& source,
             double epsilon, const boundary::BoundaryConditionSet& bc)
      : mesh_(&mesh),
        source_(&source),
        bc_(&bc),
        diffusion_(mesh, epsilon, bc) {
    if (&source.mesh() != &mesh) {
      throw std::invalid_argument("source field belongs to another mesh");
    }
  }

  void assembleMatrix(linalg::SparseMatrix& A) const;

  void assembleRhs(linalg::Vector& b) const;

  [[nodiscard]]
  linalg::LinearSystem assemble() const;

  /**
   * @brief Reports whether at least one boundary face has a Dirichlet
   * condition.
   *
   * A Poisson problem without any Dirichlet face has a constant null space and
   * therefore needs a gauge condition before it can be solved.
   *
   * @return `true` when a Dirichlet boundary face is present.
   * @throws std::out_of_range if a mesh boundary has no configured condition.
   */
  [[nodiscard]]
  bool hasDirichletBoundary() const;

  [[nodiscard]]
  const mesh::IMesh& mesh() const noexcept {
    return *mesh_;
  }

 private:
  const mesh::IMesh* mesh_;
  const field::CellField<double>* source_;
  const boundary::BoundaryConditionSet* bc_;

  operators::ScalarDiffusionFvm diffusion_;
};

}  // namespace pemu::discretization
