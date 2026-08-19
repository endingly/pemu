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
  // --------------------------------------------------------
  // Field / mesh consistency
  // --------------------------------------------------------

  if (&solution.mesh() != &discretization_.mesh()) {

    throw std::invalid_argument(
        "Poisson solution field belongs "
        "to a different mesh");
  }

  // --------------------------------------------------------
  // PDE discretization
  //
  //     -div(epsilon grad(phi)) = rho
  //
  //             ↓
  //
  //            A phi = b
  // --------------------------------------------------------

  auto system = discretization_.assemble();

  // --------------------------------------------------------
  // Symbolic analysis
  // --------------------------------------------------------

  const auto analyze_status = linear_solver_->analyzePattern(system.A);

  if (analyze_status != linalg::SolverStatus::Success) {

    return {.status = analyze_status};
  }

  // --------------------------------------------------------
  // Numeric factorization
  // --------------------------------------------------------

  const auto factor_status = linear_solver_->factorize(system.A);

  if (factor_status != linalg::SolverStatus::Success) {

    return {.status = factor_status};
  }

  // --------------------------------------------------------
  // Map CellField storage directly to Eigen.
  //
  // CellField<double> is contiguous:
  //
  //     solution.data()
  //
  // Therefore no temporary Vector + copy-back is required.
  // --------------------------------------------------------

  Eigen::Map<linalg::Vector> x(solution.data(),
                               static_cast<linalg::Index>(solution.size()));

  // --------------------------------------------------------
  // Linear solve
  // --------------------------------------------------------

  return linear_solver_->solve(system.b, x);
}

void PoissonSolver::reset() {
  linear_solver_->reset();
}

}  // namespace pemu::equation