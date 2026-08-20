// linear_solver/i_solver.hpp
#pragma once

#include "types.hpp"

namespace pemu::linalg {

enum class SolverStatus {
  Success,

  InvalidInput,
  PatternAnalysisFailed,
  FactorizationFailed,
  SolveFailed,

  Singular,
  NotPositiveDefinite,

  /// The right-hand side violates a required equation compatibility condition.
  IncompatibleRhs,

  NotAnalyzed,
  NotFactorized
};

struct SolverResult {
  SolverStatus status{SolverStatus::Success};

  double residual_norm{0.0};
  double relative_residual{0.0};

  [[nodiscard]]
  bool success() const noexcept {
    return status == SolverStatus::Success;
  }
};

class ISolver {
 public:
  virtual ~ISolver() = default;

  virtual SolverStatus analyzePattern(const SparseMatrix& A) = 0;

  virtual SolverStatus factorize(const SparseMatrix& A) = 0;

  virtual SolverResult solve(ConstVectorRef b, VectorRef x) = 0;

  virtual void reset() = 0;

  [[nodiscard]]
  virtual bool isAnalyzed() const noexcept = 0;

  [[nodiscard]]
  virtual bool isFactorized() const noexcept = 0;
};

}  // namespace pemu::linalg
