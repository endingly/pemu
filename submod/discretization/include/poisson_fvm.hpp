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
      : mesh_(&mesh), source_(&source), diffusion_(mesh, epsilon, bc) {
    if (&source.mesh() != &mesh) {
      throw std::invalid_argument("source field belongs to another mesh");
    }
  }

  void assembleMatrix(linalg::SparseMatrix& A) const;

  void assembleRhs(linalg::Vector& b) const;

  [[nodiscard]]
  linalg::LinearSystem assemble() const;

  [[nodiscard]]
  const mesh::IMesh& mesh() const noexcept {
    return *mesh_;
  }

 private:
  const mesh::IMesh* mesh_;
  const field::CellField<double>* source_;

  operators::ScalarDiffusionFvm diffusion_;
};

}  // namespace pemu::discretization