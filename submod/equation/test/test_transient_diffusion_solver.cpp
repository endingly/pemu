#include <gtest/gtest.h>
#include <pemu/equation/transient_diffusion_solver.hpp>
#include <pemu/linalg/cholmod_solver.hpp>
#include <pemu/mesh/moab_mesh.hpp>

#include <filesystem>

namespace pemu::equation::test {

namespace {

std::filesystem::path testMeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
}

class TransientDiffusionTest : public ::testing::Test {
 protected:
  TransientDiffusionTest() : mesh_(testMeshPath().string()) {}

  mesh::MoabMesh mesh_;
};

double totalMass(const mesh::IMesh& mesh, const field::CellField<double>& u) {
  double mass = 0.0;
  for (mesh::CellId cell = 0; cell < mesh.numCells(); ++cell) {
    mass += u[cell] * mesh.cellVolume(cell);
  }
  return mass;
}

double exactTransientSolution(double x, double y, double t,
                              double diffusivity) {
  constexpr double pi = std::numbers::pi_v<double>;

  return std::exp(-2.0 * pi * pi * diffusivity * t) * std::sin(pi * x) *
         std::sin(pi * y);
}

struct TransientSolveResult {
  double l2_error{};
  double final_time{};
};

TransientSolveResult solveTransientManufactured(std::size_t n, double dt,
                                                std::size_t num_steps) {
  constexpr double diffusivity = 1.0;

  // ========================================================
  // Select mesh according to requested resolution.
  //
  // n = 8  -> poisson_8x8.msh
  // n = 16 -> poisson_16x16.msh
  // ...
  // ========================================================

  const auto mesh_path =
      std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} /
      ("poisson_" + std::to_string(n) + "x" + std::to_string(n) + ".msh");
  if (!std::filesystem::exists(mesh_path)) {
    throw std::runtime_error("convergence mesh does not exist: " +
                             mesh_path.string());
  }

  mesh::MoabMesh mesh(mesh_path);

  // ========================================================
  // Initial state:
  //
  //     u(x,y,0) = sin(pi x) sin(pi y)
  // ========================================================

  field::CellField<double> u(mesh, 0.0);

  for (mesh::CellId cell = 0; cell < mesh.numCells(); ++cell) {

    const auto center = mesh.cellCenter(cell);

    u[cell] = exactTransientSolution(center.x, center.y, 0.0, diffusivity);
  }

  // ========================================================
  // Homogeneous Dirichlet BC
  // ========================================================

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(mesh::BoundaryId{1}, 0.0);

  bc.setDirichlet(mesh::BoundaryId{2}, 0.0);

  bc.setDirichlet(mesh::BoundaryId{3}, 0.0);

  bc.setDirichlet(mesh::BoundaryId{4}, 0.0);

  // ========================================================
  // Backward Euler FVM
  // ========================================================

  discretization::BackwardEulerDiffusionFvm fvm(mesh, diffusivity, dt, bc);

  equation::TransientDiffusionSolver solver(
      std::move(fvm), std::make_unique<linalg::CholmodSolver>());

  // ========================================================
  // Time integration
  // ========================================================

  for (std::size_t step = 0; step < num_steps; ++step) {

    const auto result = solver.step(u);

    if (!result.success()) {
      throw std::runtime_error("transient diffusion solve failed");
    }
  }

  const double final_time = static_cast<double>(num_steps) * dt;

  // ========================================================
  // Volume-weighted L2 error
  //
  //       sqrt(sum e_P^2 V_P)
  // ========================================================

  double error_squared = 0.0;

  for (mesh::CellId cell = 0; cell < mesh.numCells(); ++cell) {

    const auto center = mesh.cellCenter(cell);

    const double exact =
        exactTransientSolution(center.x, center.y, final_time, diffusivity);

    const double error = u[cell] - exact;

    error_squared += error * error * mesh.cellVolume(cell);
  }

  return {.l2_error = std::sqrt(error_squared),

          .final_time = final_time};
}

double convergenceRate(double coarse_error, double fine_error) {
  return std::log(coarse_error / fine_error) / std::log(2.0);
}

}  // namespace

TEST_F(TransientDiffusionTest, ConstantSolutionRemainsConstant) {
  field::CellField<double> u(mesh_, 3.0);

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(1, 3.0);
  bc.setDirichlet(2, 3.0);
  bc.setDirichlet(3, 3.0);
  bc.setDirichlet(4, 3.0);

  discretization::BackwardEulerDiffusionFvm fvm(mesh_, 1.0, 0.1, bc);

  equation::TransientDiffusionSolver solver(
      std::move(fvm), std::make_unique<linalg::CholmodSolver>());

  ASSERT_TRUE(solver.step(u).success());

  for (const double value : u) {
    EXPECT_NEAR(value, 3.0, 1e-12);
  }
}

TEST_F(TransientDiffusionTest, ZeroNeumannBoundaryConservesMass) {
  field::CellField<double> u(mesh_, 0.0);

  u[0] = 1.0;
  u[1] = 4.0;

  const double initial_mass = totalMass(mesh_, u);

  boundary::BoundaryConditionSet bc;

  bc.setNeumann(1, 0.0);
  bc.setNeumann(2, 0.0);
  bc.setNeumann(3, 0.0);
  bc.setNeumann(4, 0.0);

  discretization::BackwardEulerDiffusionFvm fvm(mesh_, 1.0, 0.1, bc);

  equation::TransientDiffusionSolver solver(
      std::move(fvm), std::make_unique<linalg::CholmodSolver>());

  for (int step = 0; step < 20; ++step) {

    ASSERT_TRUE(solver.step(u).success());
  }

  const double final_mass = totalMass(mesh_, u);

  EXPECT_NEAR(final_mass, initial_mass, 1e-11);
}

TEST_F(TransientDiffusionTest, DiffusionReducesCellDifference) {
  // --------------------------------------------------------
  // Initial state
  //
  // two_quads.msh has two unit cells.
  //
  // Let one cell start at 1 and the other at 4.
  // Diffusion should reduce their difference.
  // --------------------------------------------------------

  field::CellField<double> u(mesh_, 0.0);

  u[0] = 1.0;
  u[1] = 4.0;

  const double initial_difference = std::abs(u[1] - u[0]);

  // --------------------------------------------------------
  // Zero-Neumann boundary:
  //
  //     -D grad(u) · n = 0
  //
  // No material leaves the domain.
  // --------------------------------------------------------

  boundary::BoundaryConditionSet bc;

  bc.setNeumann(mesh::BoundaryId{1}, 0.0);

  bc.setNeumann(mesh::BoundaryId{2}, 0.0);

  bc.setNeumann(mesh::BoundaryId{3}, 0.0);

  bc.setNeumann(mesh::BoundaryId{4}, 0.0);

  // --------------------------------------------------------
  // Spatial + temporal discretization:
  //
  //     du/dt - div(D grad u) = 0
  //
  // with:
  //
  //     D  = 1
  //     dt = 0.1
  //
  // Backward Euler:
  //
  //     (M/dt + K) u^{n+1}
  //       =
  //     (M/dt) u^n
  // --------------------------------------------------------

  discretization::BackwardEulerDiffusionFvm fvm(mesh_,
                                                1.0,  // diffusivity D
                                                0.1,  // time step dt
                                                bc);

  // --------------------------------------------------------
  // Linear equation solver
  //
  // TransientDiffusionSolver owns the linear solver backend.
  //
  // On the first step it performs:
  //
  //     assembleMatrix
  //     analyzePattern
  //     factorize
  //     assembleRhs
  //     solve
  //
  // Later steps reuse the same factorization.
  // --------------------------------------------------------

  equation::TransientDiffusionSolver solver(
      std::move(fvm), std::make_unique<linalg::CholmodSolver>());

  // --------------------------------------------------------
  // Advance one time step.
  //
  // u is read as u^n while assembling RHS,
  // then overwritten in-place with u^{n+1}.
  // --------------------------------------------------------

  const auto result = solver.step(u);

  ASSERT_TRUE(result.success());

  // --------------------------------------------------------
  // Diffusion should smooth the field.
  // --------------------------------------------------------
  const double final_difference = std::abs(u[1] - u[0]);
  EXPECT_LT(final_difference, initial_difference);
  EXPECT_GT(u[0], 1.0);
  EXPECT_LT(u[1], 4.0);
  EXPECT_NEAR(u[0] + u[1], 5.0, 1e-12);
}

TEST(TransientDiffusionConvergenceTest, HasSecondOrderSpatialConvergence) {
  constexpr double dt = 1e-7;

  constexpr std::size_t steps = 100;

  const auto r8 = solveTransientManufactured(8, dt, steps);

  const auto r16 = solveTransientManufactured(16, dt, steps);

  const auto r32 = solveTransientManufactured(32, dt, steps);

  const auto r64 = solveTransientManufactured(64, dt, steps);

  const double p8_16 = convergenceRate(r8.l2_error, r16.l2_error);

  const double p16_32 = convergenceRate(r16.l2_error, r32.l2_error);

  const double p32_64 = convergenceRate(r32.l2_error, r64.l2_error);

  EXPECT_GT(p8_16, 1.9);

  EXPECT_GT(p16_32, 1.9);

  EXPECT_GT(p32_64, 1.9);

  EXPECT_NEAR(p32_64, 2.0, 0.1);

  EXPECT_GT(r8.l2_error, r16.l2_error);

  EXPECT_GT(r16.l2_error, r32.l2_error);

  EXPECT_GT(r32.l2_error, r64.l2_error);
}

TEST(TransientDiffusionConvergenceTest, HasFirstOrderTemporalConvergence) {
  constexpr std::size_t n = 64;

  const auto r1 = solveTransientManufactured(n, 1e-2, 10);

  const auto r2 = solveTransientManufactured(n, 5e-3, 20);

  const auto r3 = solveTransientManufactured(n, 2.5e-3, 40);

  EXPECT_NEAR(r1.final_time, 0.1, 1e-14);

  EXPECT_NEAR(r2.final_time, 0.1, 1e-14);

  EXPECT_NEAR(r3.final_time, 0.1, 1e-14);

  const double p12 = convergenceRate(r1.l2_error, r2.l2_error);

  const double p23 = convergenceRate(r2.l2_error, r3.l2_error);

  EXPECT_GT(p12, 0.85);

  EXPECT_GT(p23, 0.85);

  EXPECT_LT(p23, 1.15);

  EXPECT_GT(r1.l2_error, r2.l2_error);

  EXPECT_GT(r2.l2_error, r3.l2_error);
}

};  // namespace pemu::equation::test