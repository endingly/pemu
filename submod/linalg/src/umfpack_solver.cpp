#include <pemu/linalg/umfpack_solver.hpp>

namespace pemu::linalg {
namespace {

[[nodiscard]] SolverResult failure(SolverStatus status, std::string_view name,
                                   std::string_view message) noexcept {
  return {
      .status = status,
      .diagnostic = trace::makeDiagnosticEvent(trace::DiagDomain::linalg,
                                               "umfpack", name, message),
  };
}

}  // namespace

SolverResult UmfpackSolver::analyzePattern(const SparseMatrix& A) {
  analyzed_ = false;
  factorized_ = false;

  if (A.rows() == 0 || A.rows() != A.cols()) {
    return failure(SolverStatus::InvalidInput, "analyze.invalid_matrix",
                   "matrix must be non-empty and square");
  }

  solver_.analyzePattern(A);

  if (solver_.info() != Eigen::Success) {
    return failure(SolverStatus::PatternAnalysisFailed,
                   "analyze.backend_failed",
                   "UMFPACK symbolic analysis failed");
  }

  analyzed_ = true;

  return {.status = SolverStatus::Success};
}

SolverResult UmfpackSolver::factorize(const SparseMatrix& A) {
  if (!analyzed_) {
    return failure(SolverStatus::NotAnalyzed, "factorize.not_analyzed",
                   "analyzePattern must succeed before factorize");
  }

  factorized_ = false;

  solver_.factorize(A);

  if (solver_.info() == Eigen::NumericalIssue) {
    return failure(SolverStatus::Singular, "factorize.singular",
                   "matrix is singular");
  }

  if (solver_.info() != Eigen::Success) {
    return failure(SolverStatus::FactorizationFailed,
                   "factorize.backend_failed",
                   "UMFPACK numerical factorization failed");
  }

  factorized_ = true;

  return {.status = SolverStatus::Success};
}

SolverResult UmfpackSolver::solve(ConstVectorRef b, VectorRef x) {
  if (!factorized_) {
    return failure(SolverStatus::NotFactorized, "solve.not_factorized",
                   "factorize must succeed before solve");
  }

  if (b.size() != x.size()) {
    return failure(SolverStatus::InvalidInput, "solve.size_mismatch",
                   "right-hand side and solution sizes differ");
  }

  x = solver_.solve(b);

  if (solver_.info() != Eigen::Success) {
    return failure(SolverStatus::SolveFailed, "solve.backend_failed",
                   "UMFPACK solve failed");
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
