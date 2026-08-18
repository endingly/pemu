#include "umfpack_solver.hpp"

namespace pemu::linalg {

SolverStatus UmfpackSolver::analyzePattern(const SparseMatrix& A) {
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

SolverStatus UmfpackSolver::factorize(const SparseMatrix& A) {
  if (!analyzed_) {
    return SolverStatus::NotAnalyzed;
  }

  factorized_ = false;

  solver_.factorize(A);

  if (solver_.info() == Eigen::NumericalIssue) {
    return SolverStatus::Singular;
  }

  if (solver_.info() != Eigen::Success) {
    return SolverStatus::FactorizationFailed;
  }

  factorized_ = true;

  return SolverStatus::Success;
}

SolverResult UmfpackSolver::solve(ConstVectorRef b, VectorRef x) {
  if (!factorized_) {
    return {.status = SolverStatus::NotFactorized};
  }

  x = solver_.solve(b);

  if (solver_.info() != Eigen::Success) {
    return {.status = SolverStatus::SolveFailed};
  }

  return {.status = SolverStatus::Success};
}

void UmfpackSolver::reset() {
  solver_.~Backend();
  new (&solver_) Backend();

  analyzed_ = false;
  factorized_ = false;
}

};  // namespace pemu::linalg