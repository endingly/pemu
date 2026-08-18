#pragma once

#include <Eigen/UmfPackSupport>
#include "i_solver.hpp"
#include "types.hpp"

namespace pemu::linalg {

class UmfpackSolver final : public ISolver {
 public:
  using Backend = Eigen::UmfPackLU<SparseMatrix>;

  SolverStatus analyzePattern(const SparseMatrix& A) override;

  SolverStatus factorize(const SparseMatrix& A) override;

  SolverResult solve(ConstVectorRef b, VectorRef x) override;

  void reset() override;

  [[nodiscard]]
  bool isAnalyzed() const noexcept override {
    return analyzed_;
  }

  [[nodiscard]]
  bool isFactorized() const noexcept override {
    return factorized_;
  }

 private:
  Backend solver_;

  bool analyzed_{false};
  bool factorized_{false};
};

}  // namespace pemu::linalg