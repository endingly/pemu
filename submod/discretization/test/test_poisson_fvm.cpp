#include <gtest/gtest.h>
#include <filesystem>
#include <pemu/discretization/poisson_fvm.hpp>
#include <pemu/linalg/cholmod_solver.hpp>
#include <pemu/mesh/moab_mesh.hpp>

namespace pemu::discretization::test {

namespace {

std::filesystem::path testMeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
}

class PoissonFvmTest : public ::testing::Test {
 protected:
  PoissonFvmTest() : mesh_(testMeshPath().string()) {}

  mesh::MoabMesh mesh_;
};

}  // namespace

TEST_F(PoissonFvmTest, AssemblesExpectedMatrixForTwoQuads) {
  field::CellField<double> source(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(1, 0.0);
  bc.setDirichlet(2, 0.0);
  bc.setDirichlet(3, 0.0);
  bc.setDirichlet(4, 0.0);

  PoissonFvm discretization(mesh_, source, 1.0, bc);

  const auto system = discretization.assemble();

  ASSERT_EQ(system.A.rows(), 2);

  ASSERT_EQ(system.A.cols(), 2);

  EXPECT_NEAR(system.A.coeff(0, 0), 7.0, 1e-12);

  EXPECT_NEAR(system.A.coeff(1, 1), 7.0, 1e-12);

  EXPECT_NEAR(system.A.coeff(0, 1), -1.0, 1e-12);

  EXPECT_NEAR(system.A.coeff(1, 0), -1.0, 1e-12);
}

TEST_F(PoissonFvmTest, AssemblesVolumeSource) {
  field::CellField<double> source(mesh_, 0.0);

  source[0] = 2.0;
  source[1] = 3.0;

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(1, 0.0);
  bc.setDirichlet(2, 0.0);
  bc.setDirichlet(3, 0.0);
  bc.setDirichlet(4, 0.0);

  PoissonFvm discretization(mesh_, source, 1.0, bc);

  const auto system = discretization.assemble();

  EXPECT_NEAR(system.b[0], 2.0, 1e-12);

  EXPECT_NEAR(system.b[1], 3.0, 1e-12);
}

TEST_F(PoissonFvmTest, AppliesNonZeroDirichletBoundary) {
  field::CellField<double> source(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(1, 1.0);
  bc.setDirichlet(2, 0.0);
  bc.setDirichlet(3, 0.0);
  bc.setDirichlet(4, 0.0);

  PoissonFvm discretization(mesh_, source, 1.0, bc);

  const auto system = discretization.assemble();

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    const auto center = mesh_.cellCenter(cell);

    if (center.x < 1.0) {

      EXPECT_NEAR(system.b[cell], 2.0, 1e-12);

    } else {

      EXPECT_NEAR(system.b[cell], 0.0, 1e-12);
    }
  }
}

TEST_F(PoissonFvmTest, AppliesNeumannBoundaryFlux) {
  field::CellField<double> source(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  bc.setNeumann(1, 1.0);

  bc.setDirichlet(2, 0.0);

  bc.setDirichlet(3, 0.0);

  bc.setDirichlet(4, 0.0);

  PoissonFvm discretization(mesh_, source, 1.0, bc);

  const auto system = discretization.assemble();

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    const auto center = mesh_.cellCenter(cell);

    if (center.x < 1.0) {

      EXPECT_NEAR(system.b[cell], -1.0, 1e-12);

    } else {

      EXPECT_NEAR(system.b[cell], 0.0, 1e-12);
    }
  }
}

TEST_F(PoissonFvmTest, DetectsPureNeumannBoundary) {
  field::CellField<double> source(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;
  bc.setNeumann(1, 0.0);
  bc.setNeumann(2, 0.0);
  bc.setNeumann(3, 0.0);
  bc.setNeumann(4, 0.0);

  PoissonFvm discretization(mesh_, source, 1.0, bc);

  EXPECT_FALSE(discretization.hasDirichletBoundary());
}

TEST_F(PoissonFvmTest, DetectsMixedBoundaryAsHavingDirichlet) {
  field::CellField<double> source(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;
  bc.setNeumann(1, 0.0);
  bc.setDirichlet(2, 0.0);
  bc.setNeumann(3, 0.0);
  bc.setNeumann(4, 0.0);

  PoissonFvm discretization(mesh_, source, 1.0, bc);

  EXPECT_TRUE(discretization.hasDirichletBoundary());
}

TEST_F(PoissonFvmTest, MatrixIsSymmetric) {
  field::CellField<double> source(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  bc.setNeumann(1, 1.0);

  bc.setDirichlet(2, 0.0);

  bc.setDirichlet(3, 0.0);

  bc.setDirichlet(4, 0.0);

  PoissonFvm discretization(mesh_, source, 1.0, bc);

  auto system = discretization.assemble();

  linalg::SparseMatrix transpose = system.A.transpose();

  linalg::SparseMatrix difference = system.A - transpose;

  EXPECT_LT(difference.norm(), 1e-12);
}

TEST_F(PoissonFvmTest, IntergrationTest) {
  field::CellField<double> source(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  bc.setNeumann(1, 1.0);

  bc.setDirichlet(2, 0.0);

  bc.setDirichlet(3, 0.0);

  bc.setDirichlet(4, 0.0);

  PoissonFvm discretization(mesh_, source, 1.0, bc);

  auto system = discretization.assemble();

  linalg::CholmodSolver solver;

  ASSERT_EQ(solver.analyzePattern(system.A).status,
            linalg::SolverStatus::Success);

  ASSERT_EQ(solver.factorize(system.A).status, linalg::SolverStatus::Success);

  linalg::Vector x(system.b.size());

  const auto result = solver.solve(system.b, x);

  ASSERT_TRUE(result.success());
}

};  // namespace pemu::discretization::test
