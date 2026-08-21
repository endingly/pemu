// linear_solver/suitesparse/cholmod_solver.hpp
#pragma once

#include <Eigen/CholmodSupport>
#include "i_solver.hpp"

namespace pemu::linalg {

class CholmodSolver final : public ISolver {
 public:
  using Backend = Eigen::CholmodSupernodalLLT<SparseMatrix, Eigen::Lower>;

  CholmodSolver() = default;

  SolverResult analyzePattern(const SparseMatrix& A) override;

  SolverResult factorize(const SparseMatrix& A) override;

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
