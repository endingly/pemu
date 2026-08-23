#include <gtest/gtest.h>
#include <filesystem>
#include <pemu/equation/adaptive_step_explicit_species_continuity_stepper.hpp>
#include <pemu/equation/fixed_step_explicit_species_continuity_stepper.hpp>
#include <pemu/mesh/moab_mesh.hpp>

namespace pemu::equation::test {

namespace {

std::filesystem::path testMeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
}

class FixedStepExplicitSpeciesContinuityStepperTest : public ::testing::Test {
 protected:
  FixedStepExplicitSpeciesContinuityStepperTest()
      : mesh_(testMeshPath().string()) {}

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

TEST_F(FixedStepExplicitSpeciesContinuityStepperTest,
       ConstantDensityRemainsConstant) {
  constexpr double density_value = 3.0;

  field::CellField<double> density(mesh_, density_value);

  field::CellField<double> source(mesh_, 0.0);

  field::FaceField<double> velocity(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(mesh::BoundaryId{1}, density_value);

  bc.setDirichlet(mesh::BoundaryId{2}, density_value);

  bc.setDirichlet(mesh::BoundaryId{3}, density_value);

  bc.setDirichlet(mesh::BoundaryId{4}, density_value);

  equation::FixedStepExplicitSpeciesContinuityStepper stepper(mesh_, velocity,
                                                              1.0, 0.01, bc);

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

TEST_F(FixedStepExplicitSpeciesContinuityStepperTest, SourceIncreasesDensity) {
  field::CellField<double> density(mesh_, 1.0);

  field::CellField<double> source(mesh_, 2.0);

  field::FaceField<double> velocity(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(1, 1.0);
  bc.setDirichlet(2, 1.0);
  bc.setDirichlet(3, 1.0);
  bc.setDirichlet(4, 1.0);

  equation::FixedStepExplicitSpeciesContinuityStepper stepper(mesh_, velocity,
                                                              1.0, 0.1, bc);

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

TEST_F(FixedStepExplicitSpeciesContinuityStepperTest,
       SourceChangesTotalParticlesByIntegratedSource) {
  field::CellField<double> density(mesh_, 1.0);

  field::CellField<double> source(mesh_, 2.0);

  field::FaceField<double> velocity(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(1, 1.0);
  bc.setDirichlet(2, 1.0);
  bc.setDirichlet(3, 1.0);
  bc.setDirichlet(4, 1.0);

  constexpr double dt = 0.1;

  equation::FixedStepExplicitSpeciesContinuityStepper stepper(mesh_, velocity,
                                                              1.0, dt, bc);

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

TEST_F(FixedStepExplicitSpeciesContinuityStepperTest,
       NeumannOutflowContributesToFixedAndAdaptiveLossRates) {
  constexpr double density_value = 1.0;
  constexpr double outward_velocity = 2.0;
  constexpr double prescribed_outward_flux = 0.5;
  constexpr double diffusivity = 0.1;
  constexpr double dt = 0.1;
  field::CellField<double> density(mesh_, density_value);
  field::CellField<double> source(mesh_, 0.0);
  field::FaceField<double> velocity(mesh_, 0.0);
  boundary::BoundaryConditionSet conditions;
  for (mesh::BoundaryId boundary = 1; boundary <= 4; ++boundary) {
    conditions.setNeumann(boundary, prescribed_outward_flux);
  }
  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {
    if (mesh_.isBoundary(face)) {
      velocity[face] = outward_velocity;
    }
  }
  FixedStepExplicitSpeciesContinuityStepper fixed(mesh_, velocity, diffusivity,
                                                  dt, conditions);
  AdaptiveStepExplicitSpeciesContinuityStepper adaptive(
      mesh_, velocity, diffusivity, conditions);
  field::CellField<double> adaptive_loss_rate(mesh_, 0.0);
  field::CellField<double> prescribed_sink_rate(mesh_, 0.0);

  adaptive.computeTransportLossRate(adaptive_loss_rate);
  adaptive.computePrescribedBoundarySinkRate(prescribed_sink_rate);
  fixed.step(density, source);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {
    double boundary_area = 0.0;
    for (const auto face : mesh_.cellFaces(cell)) {
      if (mesh_.isBoundary(face)) {
        boundary_area += mesh_.faceArea(face);
      }
    }
    const double expected_loss_rate =
        outward_velocity * boundary_area / mesh_.cellVolume(cell);
    const double expected_prescribed_sink =
        prescribed_outward_flux * boundary_area / mesh_.cellVolume(cell);
    EXPECT_DOUBLE_EQ(adaptive_loss_rate[cell], expected_loss_rate + 0.1);
    EXPECT_DOUBLE_EQ(prescribed_sink_rate[cell], expected_prescribed_sink);
    EXPECT_DOUBLE_EQ(density[cell],
                     density_value * (1.0 - dt * expected_loss_rate) -
                         dt * expected_prescribed_sink);
  }
  EXPECT_DOUBLE_EQ(fixed.maxTransportCfl(), 0.61);
}

};  // namespace pemu::equation::test
