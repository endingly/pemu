#include <gtest/gtest.h>

#include <pemu/linalg/cholmod_solver.hpp>
#include <pemu/linalg/types.hpp>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <cmath>
#include <numbers>
#include <vector>

namespace pemu::test {

using pemu::linalg::CholmodSolver;
using pemu::linalg::Index;
using pemu::linalg::Scalar;
using pemu::linalg::SolverStatus;
using pemu::linalg::SparseMatrix;
using pemu::linalg::Vector;

using Triplet = Eigen::Triplet<Scalar, Index>;

namespace {

constexpr double pi = std::numbers::pi_v<double>;

// ------------------------------------------------------------
// Exact solution:
//
// u(x,y) = sin(pi x) sin(pi y)
//
// -laplacian(u)
//     = 2 pi^2 sin(pi x) sin(pi y)
// ------------------------------------------------------------

double exactSolution(double x, double y) {
  return std::sin(pi * x) * std::sin(pi * y);
}

double sourceTerm(double x, double y) {
  return 2.0 * pi * pi * exactSolution(x, y);
}

// ------------------------------------------------------------
// Unknowns only contain interior nodes.
//
// Grid:
//
// 0 ---- h ---- ... ---- 1
//
// total grid points per dimension:
//
//     n + 2
//
// interior unknowns:
//
//     n x n
//
// h = 1 / (n + 1)
//
// Interior (i,j):
//
//     i,j = 1 ... n
//
// Algebraic index:
//
//     k = (j-1) * n + (i-1)
// ------------------------------------------------------------

Index dof(Index i, Index j, Index n) {
  return (j - 1) * n + (i - 1);
}

struct PoissonSystem {
  SparseMatrix A;
  Vector b;

  Index n;
  double h;
};

// ------------------------------------------------------------
// Assemble:
//
// -laplacian(u) = f
//
// using the standard five-point stencil:
//
//     4 u_ij
//   - u_i-1,j
//   - u_i+1,j
//   - u_i,j-1
//   - u_i,j+1
// ------------------ = h^2 f_ij
//
// Here we multiply the whole equation by h^2.
//
// Therefore:
//
// diag = 4
// neighbor = -1
// rhs = h^2 f
//
// Boundary:
//     u = 0
//
// so boundary values contribute nothing to rhs.
// ------------------------------------------------------------

PoissonSystem assemblePoisson(Index n) {
  const Index N = n * n;

  const double h = 1.0 / static_cast<double>(n + 1);

  SparseMatrix A(N, N);
  Vector b = Vector::Zero(N);

  std::vector<Triplet> triplets;

  // About 5 nonzeros per row.
  triplets.reserve(static_cast<std::size_t>(5 * N));

  for (Index j = 1; j <= n; ++j) {
    for (Index i = 1; i <= n; ++i) {

      const Index p = dof(i, j, n);

      const double x = static_cast<double>(i) * h;

      const double y = static_cast<double>(j) * h;

      // Center coefficient.
      triplets.emplace_back(p, p, 4.0);

      // West
      if (i > 1) {
        triplets.emplace_back(p, dof(i - 1, j, n), -1.0);
      }

      // East
      if (i < n) {
        triplets.emplace_back(p, dof(i + 1, j, n), -1.0);
      }

      // South
      if (j > 1) {
        triplets.emplace_back(p, dof(i, j - 1, n), -1.0);
      }

      // North
      if (j < n) {
        triplets.emplace_back(p, dof(i, j + 1, n), -1.0);
      }

      b[p] = h * h * sourceTerm(x, y);
    }
  }

  A.setFromTriplets(triplets.begin(), triplets.end());

  A.makeCompressed();

  return {.A = std::move(A), .b = std::move(b), .n = n, .h = h};
}

Vector makeExactVector(Index n, double h) {
  Vector u(n * n);

  for (Index j = 1; j <= n; ++j) {
    for (Index i = 1; i <= n; ++i) {

      const Index p = dof(i, j, n);

      const double x = static_cast<double>(i) * h;

      const double y = static_cast<double>(j) * h;

      u[p] = exactSolution(x, y);
    }
  }

  return u;
}

// ------------------------------------------------------------
// Discrete L2 error:
//
// sqrt(
//     h^2 * sum_i |u_h - u_exact|^2
// )
//
// This approximates the continuous L2 norm.
// ------------------------------------------------------------

double l2Error(const Vector& u, const Vector& exact, double h) {
  const Vector error = u - exact;

  return h * error.norm();
}

double relativeResidual(const SparseMatrix& A, const Vector& x,
                        const Vector& b) {
  const Vector r = A * x - b;

  return r.norm() / b.norm();
}

double solvePoissonAndGetError(Index n) {
  auto system = assemblePoisson(n);

  Vector x = Vector::Zero(system.A.rows());

  CholmodSolver solver;

  EXPECT_EQ(solver.analyzePattern(system.A), SolverStatus::Success);

  EXPECT_EQ(solver.factorize(system.A), SolverStatus::Success);

  const auto result = solver.solve(system.b, x);

  EXPECT_TRUE(result.success());

  EXPECT_LT(relativeResidual(system.A, x, system.b), 1e-10);

  const Vector exact = makeExactVector(system.n, system.h);

  return l2Error(x, exact, system.h);
}

}  // namespace

TEST(PoissonTest, SolvesManufacturedSolution) {
  constexpr Index n = 32;

  auto system = assemblePoisson(n);

  Vector numerical = Vector::Zero(system.A.rows());

  CholmodSolver solver;

  ASSERT_EQ(solver.analyzePattern(system.A), SolverStatus::Success);

  ASSERT_EQ(solver.factorize(system.A), SolverStatus::Success);

  const auto result = solver.solve(system.b, numerical);

  ASSERT_TRUE(result.success());

  const Vector exact = makeExactVector(system.n, system.h);

  const double error = l2Error(numerical, exact, system.h);

  const double residual = relativeResidual(system.A, numerical, system.b);

  EXPECT_LT(residual, 1e-10);

  EXPECT_LT(error, 2e-3);
}

TEST(PoissonTest, HasSecondOrderSpatialConvergence) {
  const double e8 = solvePoissonAndGetError(8);

  const double e16 = solvePoissonAndGetError(16);

  const double e32 = solvePoissonAndGetError(32);

  const double e64 = solvePoissonAndGetError(64);

  const double p1 = std::log(e8 / e16) / std::log(2.0);

  const double p2 = std::log(e16 / e32) / std::log(2.0);

  const double p3 = std::log(e32 / e64) / std::log(2.0);

  //
  // All observed orders should already
  // clearly indicate second-order behavior.
  //
  EXPECT_GT(p1, 1.8);
  EXPECT_GT(p2, 1.9);
  EXPECT_GT(p3, 1.95);

  //
  // The observed order should converge toward 2.
  //
  EXPECT_LT(std::abs(2.0 - p2), std::abs(2.0 - p1));

  EXPECT_LT(std::abs(2.0 - p3), std::abs(2.0 - p2));
}

TEST(PoissonTest, MatrixIsSymmetric) {
  auto system = assemblePoisson(16);

  SparseMatrix diff = system.A - SparseMatrix(system.A.transpose());

  EXPECT_LT(diff.norm(), 1e-14);
}

TEST(PoissonTest, MatrixCanBeCholeskyFactorized) {
  auto system = assemblePoisson(16);

  CholmodSolver solver;

  ASSERT_EQ(solver.analyzePattern(system.A), SolverStatus::Success);

  EXPECT_EQ(solver.factorize(system.A), SolverStatus::Success);
}

TEST(PoissonTest, HasExpectedSparsityPattern) {
  constexpr Index n = 16;

  auto system = assemblePoisson(n);

  EXPECT_EQ(system.A.rows(), n * n);

  EXPECT_EQ(system.A.cols(), n * n);

  EXPECT_EQ(system.A.nonZeros(), 5 * n * n - 4 * n);
}

}  // namespace pemu::test