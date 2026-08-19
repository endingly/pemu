#include <gtest/gtest.h>
#include <filesystem>
#include <pemu/equation/poisson_solver.hpp>
#include <pemu/linalg/cholmod_solver.hpp>
#include <pemu/mesh/moab_mesh.hpp>

namespace pemu::equation::test {

namespace {

std::filesystem::path testMeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "poisson_64x64.msh";
}

class PoissonSolverTest : public ::testing::Test {
 protected:
  PoissonSolverTest() : mesh_(testMeshPath().string()) {}

  mesh::MoabMesh mesh_;
};

class CountingSolver final : public linalg::ISolver {
 public:
  int analyze_count{};
  int factorize_count{};
  int solve_count{};

  linalg::SolverStatus analyzePattern(const linalg::SparseMatrix&) override {
    ++analyze_count;

    analyzed_ = true;

    return linalg::SolverStatus::Success;
  }

  linalg::SolverStatus factorize(const linalg::SparseMatrix&) override {
    ++factorize_count;

    factorized_ = true;

    return linalg::SolverStatus::Success;
  }

  linalg::SolverResult solve(linalg::ConstVectorRef b,
                             linalg::VectorRef x) override {
    ++solve_count;

    //
    // This mock does not actually solve.
    //
    x = b;

    return {.status = linalg::SolverStatus::Success};
  }

  void reset() override {
    analyzed_ = false;
    factorized_ = false;
  }

  bool isAnalyzed() const noexcept override { return analyzed_; }

  bool isFactorized() const noexcept override { return factorized_; }

 private:
  bool analyzed_{false};
  bool factorized_{false};
};

}  // namespace

TEST_F(PoissonSolverTest, SolvesConstantDirichletSolution) {
  field::CellField<double> source(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(1, 1.0);
  bc.setDirichlet(2, 1.0);
  bc.setDirichlet(3, 1.0);
  bc.setDirichlet(4, 1.0);

  discretization::PoissonFvm poisson(mesh_, source, 1.0, bc);

  field::CellField<double> phi(mesh_, 0.0);

  equation::PoissonSolver solver(std::move(poisson),
                                 std::make_unique<linalg::CholmodSolver>());

  const auto result = solver.solve(phi);

  ASSERT_TRUE(result.success());

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(phi[cell], 1.0, 1e-12);
  }
}

TEST_F(PoissonSolverTest, WritesSolutionDirectlyIntoField) {
  field::CellField<double> source(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(1, 2.0);
  bc.setDirichlet(2, 2.0);
  bc.setDirichlet(3, 2.0);
  bc.setDirichlet(4, 2.0);

  discretization::PoissonFvm poisson(mesh_, source, 1.0, bc);

  field::CellField<double> phi(mesh_, -999.0);

  equation::PoissonSolver solver(std::move(poisson),
                                 std::make_unique<linalg::CholmodSolver>());

  ASSERT_TRUE(solver.solve(phi).success());

  for (const double value : phi) {
    EXPECT_NEAR(value, 2.0, 1e-12);
  }
}

TEST_F(PoissonSolverTest, RejectsSolutionFieldFromDifferentMesh) {
  mesh::MoabMesh mesh_a(testMeshPath().string());

  mesh::MoabMesh mesh_b(std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} /
                        "poisson_32x32.msh");

  field::CellField<double> source(mesh_a, 0.0);

  field::CellField<double> phi(mesh_b, 0.0);

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(1, 0.0);
  bc.setDirichlet(2, 0.0);
  bc.setDirichlet(3, 0.0);
  bc.setDirichlet(4, 0.0);

  discretization::PoissonFvm poisson(mesh_a, source, 1.0, bc);

  equation::PoissonSolver solver(std::move(poisson),
                                 std::make_unique<linalg::CholmodSolver>());

  EXPECT_THROW(solver.solve(phi), std::invalid_argument);
}

TEST_F(PoissonSolverTest, RejectsNullLinearSolver) {
  field::CellField<double> source(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(1, 0.0);
  bc.setDirichlet(2, 0.0);
  bc.setDirichlet(3, 0.0);
  bc.setDirichlet(4, 0.0);

  discretization::PoissonFvm poisson(mesh_, source, 1.0, bc);

  EXPECT_THROW(equation::PoissonSolver(std::move(poisson), nullptr),
               std::invalid_argument);
}

TEST_F(PoissonSolverTest, ReusesFactorizationAcrossSolves) {
  field::CellField<double> source(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(1, 0.0);
  bc.setDirichlet(2, 0.0);
  bc.setDirichlet(3, 0.0);
  bc.setDirichlet(4, 0.0);

  discretization::PoissonFvm poisson(mesh_, source, 1.0, bc);

  auto backend = std::make_unique<CountingSolver>();

  auto* raw_backend = backend.get();

  equation::PoissonSolver solver(std::move(poisson), std::move(backend));

  field::CellField<double> phi(mesh_, 0.0);

  ASSERT_TRUE(solver.solve(phi).success());

  ASSERT_TRUE(solver.solve(phi).success());

  EXPECT_EQ(raw_backend->analyze_count, 1);

  EXPECT_EQ(raw_backend->factorize_count, 1);

  EXPECT_EQ(raw_backend->solve_count, 2);
}

};  // namespace pemu::equation::test