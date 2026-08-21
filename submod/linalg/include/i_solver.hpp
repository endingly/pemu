// linear_solver/i_solver.hpp
#pragma once

#include <pemu/trace/trace.hpp>

#include "types.hpp"

#include <optional>
#include <string_view>

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
  std::optional<trace::TraceEvent> diagnostic{};

  [[nodiscard]]
  bool success() const noexcept {
    return status == SolverStatus::Success;
  }
};

class ISolver {
 public:
  virtual ~ISolver() = default;

  virtual SolverResult analyzePattern(const SparseMatrix& A) = 0;

  virtual SolverResult factorize(const SparseMatrix& A) = 0;

  virtual SolverResult solve(ConstVectorRef b, VectorRef x) = 0;

  virtual void reset() = 0;

  [[nodiscard]]
  virtual bool isAnalyzed() const noexcept = 0;

  [[nodiscard]]
  virtual bool isFactorized() const noexcept = 0;
};

}  // namespace pemu::linalg

namespace pemu {

[[nodiscard]] constexpr std::string_view to_string(
    linalg::SolverStatus status) noexcept {
  switch (status) {
    case linalg::SolverStatus::Success:
      return "Success";
    case linalg::SolverStatus::InvalidInput:
      return "InvalidInput";
    case linalg::SolverStatus::PatternAnalysisFailed:
      return "PatternAnalysisFailed";
    case linalg::SolverStatus::FactorizationFailed:
      return "FactorizationFailed";
    case linalg::SolverStatus::SolveFailed:
      return "SolveFailed";
    case linalg::SolverStatus::Singular:
      return "Singular";
    case linalg::SolverStatus::NotPositiveDefinite:
      return "NotPositiveDefinite";
    case linalg::SolverStatus::IncompatibleRhs:
      return "IncompatibleRhs";
    case linalg::SolverStatus::NotAnalyzed:
      return "NotAnalyzed";
    case linalg::SolverStatus::NotFactorized:
      return "NotFactorized";
  }
  return "Unknown";
}

}  // namespace pemu
