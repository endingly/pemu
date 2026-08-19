#pragma once

#include <pemu/discretization/operators/scalar_diffusion_fvm.hpp>
#include <pemu/field/cell_field.hpp>
#include <pemu/linalg/types.hpp>

namespace pemu::discretization {

class BackwardEulerDiffusionFvm {
 public:
  BackwardEulerDiffusionFvm(const mesh::IMesh& mesh, double diffusivity,
                            double dt, const boundary::BoundaryConditionSet& bc)
      : mesh_(&mesh), dt_(dt), diffusion_(mesh, diffusivity, bc) {
    if (dt <= 0.0) {
      throw std::invalid_argument("time step must be positive");
    }
  }

  [[nodiscard]]
  const mesh::IMesh& mesh() const noexcept {
    return *mesh_;
  }

  void assembleMatrix(linalg::SparseMatrix& A) const {
    diffusion_.assembleMatrix(A);

    //
    // A = K + M / dt
    //
    for (mesh::CellId cell = 0; cell < mesh_->numCells(); ++cell) {

      const auto i = static_cast<linalg::Index>(cell);

      A.coeffRef(i, i) += mesh_->cellVolume(cell) / dt_;
    }

    A.makeCompressed();
  }

  void assembleRhs(const field::CellField<double>& state,
                   linalg::Vector& b) const {
    if (&state.mesh() != mesh_) {
      throw std::invalid_argument("state belongs to another mesh");
    }
    diffusion_.assembleBoundaryRhs(b);
    for (mesh::CellId cell = 0; cell < mesh_->numCells(); ++cell) {
      const auto i = static_cast<linalg::Index>(cell);
      b[i] += mesh_->cellVolume(cell) / dt_ * state[cell];
    }
  }

  [[nodiscard]]
  linalg::LinearSystem assemble(const field::CellField<double>& state) const {
    linalg::LinearSystem system;

    assembleMatrix(system.A);
    assembleRhs(state, system.b);

    return system;
  }

 private:
  const mesh::IMesh* mesh_;

  double dt_;

  operators::ScalarDiffusionFvm diffusion_;
};

}  // namespace pemu::discretization