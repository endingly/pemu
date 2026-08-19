#pragma once
#include <pemu/discretization/backward_euler_fvm.hpp>
#include <pemu/linalg/i_solver.hpp>

namespace pemu::equation {

class TransientDiffusionSolver {
 public:
  TransientDiffusionSolver(
      discretization::BackwardEulerDiffusionFvm discretization,
      std::unique_ptr<linalg::ISolver> solver)
      : discretization_(std::move(discretization)), solver_(std::move(solver)) {
    if (!solver_) {
      throw std::invalid_argument("linear solver must not be null");
    }
  }

  linalg::SolverStatus initialize() {
    discretization_.assembleMatrix(A_);

    auto status = solver_->analyzePattern(A_);

    if (status != linalg::SolverStatus::Success) {
      return status;
    }

    status = solver_->factorize(A_);

    if (status != linalg::SolverStatus::Success) {
      return status;
    }

    initialized_ = true;

    return linalg::SolverStatus::Success;
  }

  linalg::SolverResult step(field::CellField<double>& state) {
    if (&state.mesh() != &discretization_.mesh()) {
      throw std::invalid_argument("state belongs to another mesh");
    }
    if (!initialized_) {
      const auto status = initialize();
      if (status != linalg::SolverStatus::Success) {
        return {.status = status};
      }
    }

    //
    // u^n -> RHS
    //
    discretization_.assembleRhs(state, b_);

    //
    // RHS is complete.
    // It is now safe to overwrite state.
    //
    Eigen::Map<linalg::Vector> x(state.data(),
                                 static_cast<linalg::Index>(state.size()));
    return solver_->solve(b_, x);
  }

 private:
  discretization::BackwardEulerDiffusionFvm discretization_;

  std::unique_ptr<linalg::ISolver> solver_;

  linalg::SparseMatrix A_;
  linalg::Vector b_;

  bool initialized_{false};
};

};  // namespace pemu::equation