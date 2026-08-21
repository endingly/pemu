#include <gtest/gtest.h>

#include <pemu/linalg/cholmod_solver.hpp>
#include <pemu/linalg/types.hpp>
#include <pemu/linalg/umfpack_solver.hpp>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <cmath>
#include <string_view>
#include <vector>

namespace pemu::linalg::test {

namespace {

using Triplet = Eigen::Triplet<Scalar, Index>;

constexpr double kTolerance = 1e-10;

// ------------------------------------------------------------
// Helper
// ------------------------------------------------------------

double relativeResidual(const SparseMatrix& A, const Vector& x,
                        const Vector& b) {
  const Vector r = A * x - b;

  const double b_norm = b.norm();

  if (b_norm == 0.0) {
    return r.norm();
  }

  return r.norm() / b_norm;
}

void expectDiagnostic(const SolverResult& result, std::string_view category,
                      std::string_view name) {
  ASSERT_TRUE(result.diagnostic.has_value());
  EXPECT_EQ(result.diagnostic->kind, trace::EventKind::Diagnostic);
  EXPECT_EQ(result.diagnostic->domain, trace::DiagDomain::linalg);
  EXPECT_EQ(result.diagnostic->category, category);
  EXPECT_EQ(result.diagnostic->name, name);
  EXPECT_FALSE(result.diagnostic->message.empty());
}

// ------------------------------------------------------------
// SPD matrix:
//
//     [ 4 -1  0 ]
// A = [-1  4 -1 ]
//     [ 0 -1  3 ]
//
// This matrix is symmetric positive definite.
// ------------------------------------------------------------

SparseMatrix makeSpdMatrix() {
  SparseMatrix A(3, 3);

  std::vector<Triplet> triplets{
      {0, 0, 4.0},  {0, 1, -1.0},

      {1, 0, -1.0}, {1, 1, 4.0},  {1, 2, -1.0},

      {2, 1, -1.0}, {2, 2, 3.0},
  };

  A.setFromTriplets(triplets.begin(), triplets.end());

  A.makeCompressed();

  return A;
}

// ------------------------------------------------------------
// General non-symmetric matrix:
//
//     [ 3  1  0 ]
// A = [ 0  4  2 ]
//     [ 1  0  5 ]
//
// Used for UMFPACK.
// ------------------------------------------------------------

SparseMatrix makeGeneralMatrix() {
  SparseMatrix A(3, 3);

  std::vector<Triplet> triplets{
      {0, 0, 3.0}, {0, 1, 1.0},

      {1, 1, 4.0}, {1, 2, 2.0},

      {2, 0, 1.0}, {2, 2, 5.0},
  };

  A.setFromTriplets(triplets.begin(), triplets.end());

  A.makeCompressed();

  return A;
}

// ------------------------------------------------------------
// Instead of inventing b directly, choose an exact solution
//
//     x_exact
//
// and compute
//
//     b = A * x_exact.
//
// Then solver output can be compared against known answer.
// ------------------------------------------------------------

Vector makeExactSolution() {
  Vector x(3);

  x << 1.0, 2.0, -1.0;

  return x;
}

}  // namespace

// ============================================================
// CHOLMOD
// ============================================================

TEST(CholmodSolverTest, SolvesSpdSystem) {
  const SparseMatrix A = makeSpdMatrix();

  const Vector x_expected = makeExactSolution();

  const Vector b = A * x_expected;

  Vector x = Vector::Zero(A.rows());

  CholmodSolver solver;

  EXPECT_FALSE(solver.isAnalyzed());
  EXPECT_FALSE(solver.isFactorized());

  ASSERT_EQ(solver.analyzePattern(A).status, SolverStatus::Success);

  EXPECT_TRUE(solver.isAnalyzed());
  EXPECT_FALSE(solver.isFactorized());

  ASSERT_EQ(solver.factorize(A).status, SolverStatus::Success);

  EXPECT_TRUE(solver.isAnalyzed());
  EXPECT_TRUE(solver.isFactorized());

  const SolverResult result = solver.solve(b, x);

  ASSERT_TRUE(result.success());
  EXPECT_FALSE(result.diagnostic.has_value());

  EXPECT_TRUE(x.isApprox(x_expected, kTolerance));

  EXPECT_LT(relativeResidual(A, x, b), kTolerance);
}

// ------------------------------------------------------------
// factorize() before analyzePattern() should be rejected.
// ------------------------------------------------------------

TEST(CholmodSolverTest, FactorizeWithoutAnalyzeFails) {
  const SparseMatrix A = makeSpdMatrix();

  CholmodSolver solver;

  const auto result = solver.factorize(A);

  EXPECT_EQ(result.status, SolverStatus::NotAnalyzed);
  expectDiagnostic(result, "cholmod", "factorize.not_analyzed");
}

// ------------------------------------------------------------
// solve() before factorization should be rejected.
// ------------------------------------------------------------

TEST(CholmodSolverTest, SolveWithoutFactorizationFails) {
  const SparseMatrix A = makeSpdMatrix();

  const Vector x_expected = makeExactSolution();

  const Vector b = A * x_expected;

  Vector x = Vector::Zero(A.rows());

  CholmodSolver solver;

  ASSERT_EQ(solver.analyzePattern(A).status, SolverStatus::Success);

  const auto result = solver.solve(b, x);

  EXPECT_EQ(result.status, SolverStatus::NotFactorized);
  expectDiagnostic(result, "cholmod", "solve.not_factorized");
}

TEST(CholmodSolverTest, RejectsMismatchedSolveDimensionsWithDiagnostic) {
  const SparseMatrix A = makeSpdMatrix();
  const Vector b = A * makeExactSolution();
  Vector x = Vector::Zero(A.rows() - 1);
  CholmodSolver solver;
  ASSERT_TRUE(solver.analyzePattern(A).success());
  ASSERT_TRUE(solver.factorize(A).success());

  const auto result = solver.solve(b, x);

  EXPECT_EQ(result.status, SolverStatus::InvalidInput);
  expectDiagnostic(result, "cholmod", "solve.size_mismatch");
}

// ------------------------------------------------------------
// Important PDE use case:
//
// sparsity pattern remains unchanged,
// but matrix values change.
//
// analyzePattern() should therefore only be needed once.
// ------------------------------------------------------------

TEST(CholmodSolverTest, ReusesPatternForNewMatrixValues) {
  SparseMatrix A = makeSpdMatrix();

  CholmodSolver solver;

  ASSERT_EQ(solver.analyzePattern(A).status, SolverStatus::Success);

  // ---------------- first system ----------------

  {
    const Vector x_expected = makeExactSolution();

    const Vector b = A * x_expected;

    Vector x = Vector::Zero(A.rows());

    ASSERT_EQ(solver.factorize(A).status, SolverStatus::Success);

    ASSERT_TRUE(solver.solve(b, x).success());

    EXPECT_TRUE(x.isApprox(x_expected, kTolerance));
  }

  // ---------------- change values only ----------------
  //
  // Same sparsity pattern:
  //
  // [ 6 -1  0 ]
  // [-1  5 -1 ]
  // [ 0 -1  4 ]
  //

  A.coeffRef(0, 0) = 6.0;
  A.coeffRef(1, 1) = 5.0;
  A.coeffRef(2, 2) = 4.0;

  // No analyzePattern() here.

  ASSERT_EQ(solver.factorize(A).status, SolverStatus::Success);

  const Vector x_expected = makeExactSolution();

  const Vector b = A * x_expected;

  Vector x = Vector::Zero(A.rows());

  ASSERT_TRUE(solver.solve(b, x).success());

  EXPECT_TRUE(x.isApprox(x_expected, kTolerance));

  EXPECT_LT(relativeResidual(A, x, b), kTolerance);
}

// ------------------------------------------------------------
// reset() should invalidate previous solver state.
// ------------------------------------------------------------

TEST(CholmodSolverTest, ResetClearsSolverState) {
  const SparseMatrix A = makeSpdMatrix();

  CholmodSolver solver;

  ASSERT_EQ(solver.analyzePattern(A).status, SolverStatus::Success);

  ASSERT_EQ(solver.factorize(A).status, SolverStatus::Success);

  ASSERT_TRUE(solver.isAnalyzed());
  ASSERT_TRUE(solver.isFactorized());

  solver.reset();

  EXPECT_FALSE(solver.isAnalyzed());
  EXPECT_FALSE(solver.isFactorized());
}

// ============================================================
// UMFPACK
// ============================================================

TEST(UmfpackSolverTest, SolvesGeneralSparseSystem) {
  const SparseMatrix A = makeGeneralMatrix();

  const Vector x_expected = makeExactSolution();

  const Vector b = A * x_expected;

  Vector x = Vector::Zero(A.rows());

  UmfpackSolver solver;

  ASSERT_EQ(solver.analyzePattern(A).status, SolverStatus::Success);

  ASSERT_EQ(solver.factorize(A).status, SolverStatus::Success);

  const SolverResult result = solver.solve(b, x);

  ASSERT_TRUE(result.success());
  EXPECT_FALSE(result.diagnostic.has_value());

  EXPECT_TRUE(x.isApprox(x_expected, kTolerance));

  EXPECT_LT(relativeResidual(A, x, b), kTolerance);
}

TEST(UmfpackSolverTest, FactorizeWithoutAnalyzeFails) {
  const SparseMatrix A = makeGeneralMatrix();

  UmfpackSolver solver;

  const auto result = solver.factorize(A);

  EXPECT_EQ(result.status, SolverStatus::NotAnalyzed);
  expectDiagnostic(result, "umfpack", "factorize.not_analyzed");
}

TEST(UmfpackSolverTest, SolveWithoutFactorizationFails) {
  const SparseMatrix A = makeGeneralMatrix();

  const Vector x_expected = makeExactSolution();

  const Vector b = A * x_expected;

  Vector x = Vector::Zero(A.rows());

  UmfpackSolver solver;

  ASSERT_EQ(solver.analyzePattern(A).status, SolverStatus::Success);

  const SolverResult result = solver.solve(b, x);

  EXPECT_EQ(result.status, SolverStatus::NotFactorized);
  expectDiagnostic(result, "umfpack", "solve.not_factorized");
}

TEST(UmfpackSolverTest, RejectsMismatchedSolveDimensionsWithDiagnostic) {
  const SparseMatrix A = makeGeneralMatrix();
  const Vector b = A * makeExactSolution();
  Vector x = Vector::Zero(A.rows() - 1);
  UmfpackSolver solver;
  ASSERT_TRUE(solver.analyzePattern(A).success());
  ASSERT_TRUE(solver.factorize(A).success());

  const auto result = solver.solve(b, x);

  EXPECT_EQ(result.status, SolverStatus::InvalidInput);
  expectDiagnostic(result, "umfpack", "solve.size_mismatch");
}

// ------------------------------------------------------------
// Same pattern, different values.
// ------------------------------------------------------------

TEST(UmfpackSolverTest, ReusesPatternForNewMatrixValues) {
  SparseMatrix A = makeGeneralMatrix();

  UmfpackSolver solver;

  ASSERT_EQ(solver.analyzePattern(A).status, SolverStatus::Success);

  // First numerical factorization.

  ASSERT_EQ(solver.factorize(A).status, SolverStatus::Success);

  {
    const Vector expected = makeExactSolution();

    const Vector b = A * expected;

    Vector x = Vector::Zero(A.rows());

    ASSERT_TRUE(solver.solve(b, x).success());

    EXPECT_TRUE(x.isApprox(expected, kTolerance));
  }

  // Change numerical values while keeping structure.

  A.coeffRef(0, 0) = 5.0;
  A.coeffRef(0, 1) = 2.0;

  A.coeffRef(1, 1) = 6.0;
  A.coeffRef(1, 2) = 1.0;

  A.coeffRef(2, 0) = 3.0;
  A.coeffRef(2, 2) = 7.0;

  ASSERT_EQ(solver.factorize(A).status, SolverStatus::Success);

  const Vector expected = makeExactSolution();

  const Vector b = A * expected;

  Vector x = Vector::Zero(A.rows());

  ASSERT_TRUE(solver.solve(b, x).success());

  EXPECT_TRUE(x.isApprox(expected, kTolerance));

  EXPECT_LT(relativeResidual(A, x, b), kTolerance);
}

TEST(UmfpackSolverTest, ResetClearsSolverState) {
  const SparseMatrix A = makeGeneralMatrix();

  UmfpackSolver solver;

  ASSERT_EQ(solver.analyzePattern(A).status, SolverStatus::Success);

  ASSERT_EQ(solver.factorize(A).status, SolverStatus::Success);

  solver.reset();

  EXPECT_FALSE(solver.isAnalyzed());
  EXPECT_FALSE(solver.isFactorized());
}

}  // namespace pemu::linalg::test
