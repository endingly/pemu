#include <gtest/gtest.h>
#include <filesystem>
#include <pemu/equation/explicit_species_continuity_stepper.hpp>
#include <pemu/mesh/moab_mesh.hpp>

namespace pemu::equation::test {

namespace {

std::filesystem::path testMeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
}

class SpeciesContinuityTest : public ::testing::Test {
 protected:
  SpeciesContinuityTest() : mesh_(testMeshPath().string()) {}

  mesh::MoabMesh mesh_;
};

double totalParticles(const mesh::IMesh& mesh,
                      const field::CellField<double>& density) {
  double total = 0.0;

  for (mesh::CellId cell = 0; cell < mesh.numCells(); ++cell) {

    total += density[cell] * mesh.cellVolume(cell);
  }

  return total;
}

};  // namespace

TEST_F(SpeciesContinuityTest, ConstantDensityRemainsConstant) {
  constexpr double density_value = 3.0;

  field::CellField<double> density(mesh_, density_value);

  field::CellField<double> source(mesh_, 0.0);

  field::FaceField<double> velocity(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(mesh::BoundaryId{1}, density_value);

  bc.setDirichlet(mesh::BoundaryId{2}, density_value);

  bc.setDirichlet(mesh::BoundaryId{3}, density_value);

  bc.setDirichlet(mesh::BoundaryId{4}, density_value);

  equation::ExplicitSpeciesContinuityStepper stepper(mesh_, velocity, 1.0, 0.01,
                                                     bc);

  if (stepper.maxTransportCfl() > 1.0) {
    throw std::runtime_error(
        "explicit SG transport CFL "
        "condition violated");
  }

  stepper.step(density, source);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(density[cell], density_value, 1e-12);
  }
}

TEST_F(SpeciesContinuityTest, SourceIncreasesDensity) {
  field::CellField<double> density(mesh_, 1.0);

  field::CellField<double> source(mesh_, 2.0);

  field::FaceField<double> velocity(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(1, 1.0);
  bc.setDirichlet(2, 1.0);
  bc.setDirichlet(3, 1.0);
  bc.setDirichlet(4, 1.0);

  equation::ExplicitSpeciesContinuityStepper stepper(mesh_, velocity, 1.0, 0.1,
                                                     bc);

  if (stepper.maxTransportCfl() > 1.0) {
    throw std::runtime_error(
        "explicit SG transport CFL "
        "condition violated");
  }
  stepper.step(density, source);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(density[cell], 1.2, 1e-12);
  }
}

TEST_F(SpeciesContinuityTest, SourceChangesTotalParticlesByIntegratedSource) {
  field::CellField<double> density(mesh_, 1.0);

  field::CellField<double> source(mesh_, 2.0);

  field::FaceField<double> velocity(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(1, 1.0);
  bc.setDirichlet(2, 1.0);
  bc.setDirichlet(3, 1.0);
  bc.setDirichlet(4, 1.0);

  constexpr double dt = 0.1;

  equation::ExplicitSpeciesContinuityStepper stepper(mesh_, velocity, 1.0, dt,
                                                     bc);

  const double before = totalParticles(mesh_, density);

  double integrated_source = 0.0;

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    integrated_source += source[cell] * mesh_.cellVolume(cell);
  }

  if (stepper.maxTransportCfl() > 1.0) {
    throw std::runtime_error(
        "explicit SG transport CFL "
        "condition violated");
  }

  stepper.step(density, source);

  const double after = totalParticles(mesh_, density);

  EXPECT_NEAR(after - before, dt * integrated_source, 1e-12);
}

};  // namespace pemu::equation::test