#include <pemu/equation/poisson_solver.hpp>

#include <Eigen/Core>

#include <stdexcept>
#include <utility>

namespace pemu::equation {

PoissonSolver::PoissonSolver(discretization::PoissonFvm discretization,
                             std::unique_ptr<linalg::ISolver> linear_solver)
    : discretization_(std::move(discretization)),
      linear_solver_(std::move(linear_solver)) {
  if (!linear_solver_) {
    throw std::invalid_argument("PoissonSolver requires a linear solver");
  }
}

linalg::SolverResult PoissonSolver::solve(field::CellField<double>& solution) {
  if (&solution.mesh() != &discretization_.mesh()) {
    throw std::invalid_argument(
        "Poisson solution field "
        "belongs to another mesh");
  }
  if (!initialized_) {
    const auto status = initialize();
    if (status != linalg::SolverStatus::Success) {
      return {.status = status};
    }
  }
  // --------------------------------------------------------
  // Only RHS is rebuilt.
  // --------------------------------------------------------
  discretization_.assembleRhs(b_);
  // --------------------------------------------------------
  // Solve directly into CellField storage.
  // --------------------------------------------------------
  Eigen::Map<linalg::Vector> x(solution.data(),
                               static_cast<linalg::Index>(solution.size()));
  return linear_solver_->solve(b_, x);
}

void PoissonSolver::reset() {
  linear_solver_->reset();
  A_.resize(0, 0);
  b_.resize(0);
  initialized_ = false;
}

linalg::SolverStatus PoissonSolver::initialize() {
  discretization_.assembleMatrix(A_);

  const auto analyze_status = linear_solver_->analyzePattern(A_);

  if (analyze_status != linalg::SolverStatus::Success) {

    initialized_ = false;

    return analyze_status;
  }

  const auto factor_status = linear_solver_->factorize(A_);

  if (factor_status != linalg::SolverStatus::Success) {

    initialized_ = false;

    return factor_status;
  }

  initialized_ = true;

  return linalg::SolverStatus::Success;
}

}  // namespace pemu::equation