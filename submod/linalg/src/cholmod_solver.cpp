// linear_solver/suitesparse/cholmod_solver.cpp

#include <pemu/linalg/cholmod_solver.hpp>

namespace pemu::linalg {

SolverStatus CholmodSolver::analyzePattern(const SparseMatrix& A) {
  analyzed_ = false;
  factorized_ = false;

  if (A.rows() == 0 || A.rows() != A.cols()) {
    return SolverStatus::InvalidInput;
  }

  solver_.analyzePattern(A);

  if (solver_.info() != Eigen::Success) {
    return SolverStatus::PatternAnalysisFailed;
  }

  analyzed_ = true;

  return SolverStatus::Success;
}

SolverStatus CholmodSolver::factorize(const SparseMatrix& A) {
  if (!analyzed_) {
    return SolverStatus::NotAnalyzed;
  }

  factorized_ = false;

  solver_.factorize(A);

  if (solver_.info() == Eigen::NumericalIssue) {
    return SolverStatus::NotPositiveDefinite;
  }

  if (solver_.info() != Eigen::Success) {
    return SolverStatus::FactorizationFailed;
  }

  factorized_ = true;

  return SolverStatus::Success;
}

SolverResult CholmodSolver::solve(ConstVectorRef b, VectorRef x) {
  if (!factorized_) {
    return {.status = SolverStatus::NotFactorized};
  }

  if (b.size() != x.size()) {
    return {.status = SolverStatus::InvalidInput};
  }

  x = solver_.solve(b);

  if (solver_.info() != Eigen::Success) {
    return {.status = SolverStatus::SolveFailed};
  }

  return {.status = SolverStatus::Success};
}

void CholmodSolver::reset() {
  // Eigen::CholmodSupernodalLLT has a deleted copy-assignment operator,
  // so re-create the object in place instead of assigning a temporary.
  solver_.~Backend();
  new (&solver_) Backend();

  analyzed_ = false;
  factorized_ = false;
}

}  // namespace pemu::linalg