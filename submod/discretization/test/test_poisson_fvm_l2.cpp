#include <gtest/gtest.h>
#include <cmath>
#include <filesystem>
#include <numbers>

#include <pemu/boundary/boundary_condition_set.hpp>
#include <pemu/discretization/poisson_fvm.hpp>
#include <pemu/field/cell_field.hpp>
#include <pemu/linalg/cholmod_solver.hpp>
#include <pemu/mesh/moab_mesh.hpp>

namespace pemu::discretization::test {

namespace {

double exactSolution(const double x, const double y) {
  constexpr double pi = std::numbers::pi_v<double>;
  return std::sin(pi * x) * std::sin(pi * y);
}

double sourceTerm(const double x, const double y) {
  constexpr double pi = std::numbers::pi_v<double>;
  return 2.0 * pi * pi * exactSolution(x, y);
}

double convergenceRate(const double coarse_error, const double fine_error) {
  return std::log(coarse_error / fine_error) / std::log(2.0);
}

struct PoissonSolveResult {
  double l2_error{};
  double relative_residual{};
};

std::filesystem::path poissonMeshPath(const std::size_t n) {
  const auto path =
      std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} /
      ("poisson_" + std::to_string(n) + "x" + std::to_string(n) + ".msh");

  if (!std::filesystem::exists(path)) {
    throw std::runtime_error("Poisson test mesh does not exist: " +
                             path.string());
  }

  return path;
}

PoissonSolveResult solveManufacturedPoisson(const std::size_t n) {
  // ========================================================
  // Mesh
  //
  // Read a fixed Gmsh test fixture.
  // ========================================================

  mesh::MoabMesh mesh(poissonMeshPath(n).string());

  // ========================================================
  // Basic fixture sanity checks.
  // ========================================================

  const std::size_t expected_cells = n * n;

  if (mesh.numCells() != expected_cells) {
    throw std::runtime_error("unexpected number of cells in Poisson test mesh");
  }

  // ========================================================
  // Source
  //
  // -Delta(phi) = f
  //
  // phi = sin(pi x) sin(pi y)
  //
  // f = 2 pi^2 sin(pi x) sin(pi y)
  // ========================================================

  field::CellField<double> source(mesh);

  for (mesh::CellId cell = 0; cell < mesh.numCells(); ++cell) {

    const auto center = mesh.cellCenter(cell);

    source[cell] = sourceTerm(center.x, center.y);
  }

  // ========================================================
  // Boundary conditions
  //
  // Physical groups:
  //
  // left   = 1
  // right  = 2
  // bottom = 3
  // top    = 4
  //
  // Exact solution is zero on the entire boundary.
  // ========================================================

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(mesh::BoundaryId{1}, 0.0);

  bc.setDirichlet(mesh::BoundaryId{2}, 0.0);

  bc.setDirichlet(mesh::BoundaryId{3}, 0.0);

  bc.setDirichlet(mesh::BoundaryId{4}, 0.0);

  // ========================================================
  // Discretization
  // ========================================================

  discretization::PoissonFvm poisson(mesh, source, 1.0, bc);

  auto system = poisson.assemble();

  // ========================================================
  // Linear solver
  // ========================================================

  linalg::CholmodSolver solver;

  const auto analyze_result = solver.analyzePattern(system.A);

  if (!analyze_result.success()) {

    throw std::runtime_error("CHOLMOD pattern analysis failed");
  }

  const auto factorize_result = solver.factorize(system.A);

  if (!factorize_result.success()) {

    throw std::runtime_error("CHOLMOD factorization failed");
  }

  linalg::Vector solution(system.b.size());

  const auto solve_result = solver.solve(system.b, solution);

  if (!solve_result.success()) {
    throw std::runtime_error("CHOLMOD solve failed");
  }

  // ========================================================
  // L2 discretization error
  //
  // ||e||_L2
  //
  //   ~= sqrt(
  //        sum_P e_P^2 V_P
  //      )
  // ========================================================

  double error_squared = 0.0;

  for (mesh::CellId cell = 0; cell < mesh.numCells(); ++cell) {

    const auto center = mesh.cellCenter(cell);

    const double exact = exactSolution(center.x, center.y);

    const auto row = static_cast<linalg::Index>(cell);

    const double error = solution[row] - exact;

    error_squared += error * error * mesh.cellVolume(cell);
  }

  const double l2_error = std::sqrt(error_squared);

  // ========================================================
  // Linear algebra residual
  // ========================================================

  const linalg::Vector residual = system.A * solution - system.b;

  const double rhs_norm = system.b.norm();

  const double relative_residual =
      rhs_norm > 0.0 ? residual.norm() / rhs_norm : residual.norm();

  return {.l2_error = l2_error,

          .relative_residual = relative_residual};
}

};  // namespace

TEST(PoissonFvmConvergenceTest, MeshFixturesHaveExpectedCellCounts) {
  const std::array<std::size_t, 4> resolutions{8, 16, 32, 64};

  for (const auto n : resolutions) {

    mesh::MoabMesh mesh(poissonMeshPath(n).string());

    EXPECT_EQ(mesh.numCells(), n * n);
  }
}

TEST(PoissonFvmConvergenceTest, MeshFixturesHaveExpectedBoundaryFaceCounts) {
  const std::array<std::size_t, 4> resolutions{8, 16, 32, 64};

  for (const auto n : resolutions) {

    mesh::MoabMesh mesh(poissonMeshPath(n).string());

    std::size_t boundary_faces = 0;

    for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {

      if (mesh.isBoundary(face)) {
        ++boundary_faces;
      }
    }

    EXPECT_EQ(boundary_faces, 4 * n);
  }
}

TEST(PoissonFvmConvergenceTest, HasSecondOrderSpatialConvergence) {
  const auto r8 = solveManufacturedPoisson(8);

  const auto r16 = solveManufacturedPoisson(16);

  const auto r32 = solveManufacturedPoisson(32);

  const auto r64 = solveManufacturedPoisson(64);

  const double p8_16 = convergenceRate(r8.l2_error, r16.l2_error);

  const double p16_32 = convergenceRate(r16.l2_error, r32.l2_error);

  const double p32_64 = convergenceRate(r32.l2_error, r64.l2_error);

  EXPECT_GT(p8_16, 1.9);

  EXPECT_GT(p16_32, 1.95);

  EXPECT_GT(p32_64, 1.98);

  EXPECT_NEAR(p32_64, 2.0, 0.05);
}

TEST(PoissonFvmConvergenceTest, LinearSolverResidualIsSmall) {
  const auto result = solveManufacturedPoisson(32);

  EXPECT_LT(result.relative_residual, 1e-10);
}

TEST(PoissonFvmConvergenceTest, ErrorDecreasesUnderMeshRefinement) {
  const auto r8 = solveManufacturedPoisson(8);

  const auto r16 = solveManufacturedPoisson(16);

  const auto r32 = solveManufacturedPoisson(32);

  const auto r64 = solveManufacturedPoisson(64);

  EXPECT_GT(r8.l2_error, r16.l2_error);

  EXPECT_GT(r16.l2_error, r32.l2_error);

  EXPECT_GT(r32.l2_error, r64.l2_error);
}

TEST(PoissonFvmConvergenceTest, ErrorReducesApproximatelyByFactorFour) {
  const auto r16 = solveManufacturedPoisson(16);

  const auto r32 = solveManufacturedPoisson(32);

  const auto r64 = solveManufacturedPoisson(64);

  EXPECT_NEAR(r16.l2_error / r32.l2_error, 4.0, 0.15);

  EXPECT_NEAR(r32.l2_error / r64.l2_error, 4.0, 0.10);
}

}  // namespace pemu::discretization::test
