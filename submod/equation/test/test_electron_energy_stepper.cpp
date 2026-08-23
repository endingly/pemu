#include <gtest/gtest.h>

#include <pemu/equation/explicit_electron_energy_stepper.hpp>
#include <pemu/mesh/moab_mesh.hpp>
#include <pemu/physics/electron_energy.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>

namespace pemu::equation::test {
namespace {

/** @brief Returns the two-cell mesh used by electron-energy unit tests. */
[[nodiscard]] std::filesystem::path testMeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
}

/** @brief Creates the same Dirichlet energy density on every boundary. */
[[nodiscard]] boundary::BoundaryConditionSet uniformBoundary(double value) {
  boundary::BoundaryConditionSet conditions;
  conditions.setDirichlet(mesh::BoundaryId{1}, value);
  conditions.setDirichlet(mesh::BoundaryId{2}, value);
  conditions.setDirichlet(mesh::BoundaryId{3}, value);
  conditions.setDirichlet(mesh::BoundaryId{4}, value);
  return conditions;
}

/** @brief Creates one outward diffusive energy flux on every boundary. */
[[nodiscard]] boundary::BoundaryConditionSet uniformNeumannBoundary(
    double value = 0.0) {
  boundary::BoundaryConditionSet conditions;
  for (mesh::BoundaryId boundary = 1; boundary <= 4; ++boundary) {
    conditions.setNeumann(boundary, value);
  }
  return conditions;
}

/** @brief Creates symmetric boundary data for cell energies one and three. */
[[nodiscard]] boundary::BoundaryConditionSet symmetricGradientBoundary() {
  boundary::BoundaryConditionSet conditions;
  conditions.setDirichlet(mesh::BoundaryId{1}, 1.0);
  conditions.setDirichlet(mesh::BoundaryId{2}, 3.0);
  conditions.setDirichlet(mesh::BoundaryId{3}, 2.0);
  conditions.setDirichlet(mesh::BoundaryId{4}, 2.0);
  return conditions;
}

/** @brief Integrates energy density over all control volumes. */
[[nodiscard]] double integratedEnergy(
    const mesh::IMesh& mesh, const field::CellField<double>& energy_density) {
  double result{};
  for (mesh::CellId cell = 0; cell < mesh.numCells(); ++cell) {
    result += energy_density[cell] * mesh.cellVolume(cell);
  }
  return result;
}

class ExplicitElectronEnergyStepperTest : public ::testing::Test {
 protected:
  /** @brief Loads the shared orthogonal two-cell fixture mesh. */
  ExplicitElectronEnergyStepperTest() : mesh_(testMeshPath().string()) {}

  mesh::MoabMesh mesh_;
};

}  // namespace

TEST_F(ExplicitElectronEnergyStepperTest,
       ConvertsBetweenEnergyDensityAndMeanEnergy) {
  field::CellField<double> number_density(mesh_, 0.0);
  number_density[0] = 2.0;
  number_density[1] = 0.0;
  field::CellField<double> mean_energy(mesh_, 0.0);
  mean_energy[0] = 3.0;
  mean_energy[1] = 7.0;
  field::CellField<double> energy_density(mesh_, 0.0);

  physics::computeElectronEnergyDensity(number_density, mean_energy,
                                        energy_density);
  EXPECT_DOUBLE_EQ(energy_density[0], 6.0);
  EXPECT_DOUBLE_EQ(energy_density[1], 0.0);

  energy_density[1] = 5.0;
  physics::computeElectronMeanEnergy(energy_density, number_density, 1e-12,
                                     mean_energy);
  EXPECT_DOUBLE_EQ(mean_energy[0], 3.0);
  EXPECT_DOUBLE_EQ(mean_energy[1], 0.0);
}

TEST_F(ExplicitElectronEnergyStepperTest,
       MaxwellianClosureScalesCompleteTransportOperator) {
  field::FaceField<double> electron_velocity(mesh_, 0.0);
  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {
    electron_velocity[face] = 0.4 * mesh_.faceNormal(face).x;
  }
  field::CellField<double> energy_density(mesh_, 0.0);
  energy_density[0] = 1.0;
  energy_density[1] = 3.0;
  field::CellField<double> source(mesh_, 0.0);
  field::CellField<double> unit_increment(mesh_, 0.0);
  field::CellField<double> maxwellian_increment(mesh_, 0.0);
  constexpr double electron_diffusivity = 0.2;
  constexpr double dt = 0.01;
  ExplicitElectronEnergyStepper unit_stepper(mesh_, electron_velocity,
                                             electron_diffusivity,
                                             uniformBoundary(2.0), 1.0);
  ExplicitElectronEnergyStepper maxwellian_stepper(
      mesh_, electron_velocity, electron_diffusivity, uniformBoundary(2.0));

  unit_stepper.computeIncrement(energy_density, source, dt, unit_increment);
  maxwellian_stepper.computeIncrement(energy_density, source, dt,
                                      maxwellian_increment);

  EXPECT_DOUBLE_EQ(maxwellian_stepper.energyTransportFactor(), 5.0 / 3.0);
  EXPECT_DOUBLE_EQ(maxwellian_stepper.energyDiffusivity(),
                   (5.0 / 3.0) * electron_diffusivity);
  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {
    EXPECT_NEAR(maxwellian_increment[cell], (5.0 / 3.0) * unit_increment[cell],
                1e-14);
  }
}

TEST_F(ExplicitElectronEnergyStepperTest,
       RejectsInvalidEnergyBoundaryAtConstruction) {
  field::FaceField<double> electron_velocity(mesh_, 0.0);
  auto boundary_conditions = uniformBoundary(1.0);
  boundary_conditions.setDirichlet(mesh::BoundaryId{1}, -1.0);

  EXPECT_THROW(ExplicitElectronEnergyStepper(mesh_, electron_velocity, 0.1,
                                             boundary_conditions),
               std::invalid_argument);
}

TEST_F(ExplicitElectronEnergyStepperTest,
       ConstantEnergyDensityRemainsConstant) {
  constexpr double value = 4.0;
  field::FaceField<double> electron_velocity(mesh_, 0.0);
  field::CellField<double> energy_density(mesh_, value);
  field::CellField<double> source(mesh_, 0.0);
  ExplicitElectronEnergyStepper stepper(mesh_, electron_velocity, 0.1,
                                        uniformBoundary(value));
  const double dt = 0.5 * stepper.maxStableTransportTimeStep();

  stepper.step(energy_density, source, dt);

  for (const double cell_energy : energy_density) {
    EXPECT_NEAR(cell_energy, value, 1e-14);
  }
}

TEST_F(ExplicitElectronEnergyStepperTest,
       HomogeneousNeumannPreservesConstantZeroDriftEnergy) {
  constexpr double value = 4.0;
  field::FaceField<double> electron_velocity(mesh_, 0.0);
  field::CellField<double> energy_density(mesh_, value);
  field::CellField<double> source(mesh_, 0.0);
  ExplicitElectronEnergyStepper stepper(mesh_, electron_velocity, 0.1,
                                        uniformNeumannBoundary());
  const double dt = 0.5 * stepper.maxStableTransportTimeStep();

  stepper.step(energy_density, source, dt);

  for (const double cell_energy : energy_density) {
    EXPECT_NEAR(cell_energy, value, 1e-14);
  }
}

TEST_F(ExplicitElectronEnergyStepperTest,
       PrescribedNeumannOutflowLimitsPositiveTimeStep) {
  field::FaceField<double> electron_velocity(mesh_, 0.0);
  field::CellField<double> energy_density(mesh_, 1.0);
  field::CellField<double> source(mesh_, 0.0);
  ExplicitElectronEnergyStepper stepper(mesh_, electron_velocity, 0.1,
                                        uniformNeumannBoundary(0.5));

  // Each unit cell has three boundary edges, so prescribed depletion is 1.5.
  // The Maxwellian energy diffusivity contributes the conservative internal
  // diagonal 1/6, giving dt_positive = 1 / (1.5 + 1/6) = 0.6.
  EXPECT_NEAR(stepper.maxPositiveTimeStep(energy_density, source), 0.6, 1e-14);
}

TEST_F(ExplicitElectronEnergyStepperTest,
       SourceChangesIntegratedEnergyByIntegratedSource) {
  constexpr double dt = 0.1;
  field::FaceField<double> electron_velocity(mesh_, 0.0);
  field::CellField<double> energy_density(mesh_, 2.0);
  field::CellField<double> source(mesh_, 3.0);
  ExplicitElectronEnergyStepper stepper(mesh_, electron_velocity, 0.1,
                                        uniformBoundary(2.0));
  const double before = integratedEnergy(mesh_, energy_density);
  const double integrated_source = integratedEnergy(mesh_, source);

  stepper.step(energy_density, source, dt);

  EXPECT_NEAR(integratedEnergy(mesh_, energy_density) - before,
              dt * integrated_source, 1e-13);
}

TEST_F(ExplicitElectronEnergyStepperTest,
       DiffusionRedistributesSymmetricEnergyConservatively) {
  field::FaceField<double> electron_velocity(mesh_, 0.0);
  field::CellField<double> energy_density(mesh_, 0.0);
  energy_density[0] = 1.0;
  energy_density[1] = 3.0;
  field::CellField<double> source(mesh_, 0.0);
  ExplicitElectronEnergyStepper stepper(mesh_, electron_velocity, 0.1,
                                        symmetricGradientBoundary());
  const double before = integratedEnergy(mesh_, energy_density);
  const double dt = 0.25 * stepper.maxStableTransportTimeStep();

  stepper.step(energy_density, source, dt);

  EXPECT_GT(energy_density[0], 1.0);
  EXPECT_LT(energy_density[1], 3.0);
  EXPECT_NEAR(integratedEnergy(mesh_, energy_density), before, 1e-13);
}

TEST_F(ExplicitElectronEnergyStepperTest,
       RejectsTransportCflViolationWithoutChangingEnergy) {
  field::FaceField<double> electron_velocity(mesh_, 0.0);
  field::CellField<double> energy_density(mesh_, 1.0);
  field::CellField<double> source(mesh_, 0.0);
  ExplicitElectronEnergyStepper stepper(mesh_, electron_velocity, 1.0,
                                        uniformBoundary(1.0));
  const double unstable_dt = 1.01 * stepper.maxStableTransportTimeStep();

  EXPECT_THROW(stepper.step(energy_density, source, unstable_dt),
               std::runtime_error);
  for (const double cell_energy : energy_density) {
    EXPECT_DOUBLE_EQ(cell_energy, 1.0);
  }
}

TEST_F(ExplicitElectronEnergyStepperTest,
       RejectsNegativeSourceUpdateWithoutChangingEnergy) {
  field::FaceField<double> electron_velocity(mesh_, 0.0);
  field::CellField<double> energy_density(mesh_, 1.0);
  field::CellField<double> source(mesh_, 0.0);
  ExplicitElectronEnergyStepper stepper(mesh_, electron_velocity, 0.01,
                                        uniformBoundary(1.0));
  const double dt = std::min(0.1, 0.5 * stepper.maxStableTransportTimeStep());
  source.fill(-2.0 / dt);
  const double positivity_limit =
      stepper.maxPositiveTimeStep(energy_density, source);

  EXPECT_LT(positivity_limit, dt);
  EXPECT_THROW(stepper.step(energy_density, source, dt), std::runtime_error);
  for (const double cell_energy : energy_density) {
    EXPECT_DOUBLE_EQ(cell_energy, 1.0);
  }
}

TEST_F(ExplicitElectronEnergyStepperTest,
       StableIncrementSupportsValidationBeforeCoupledCommit) {
  field::FaceField<double> electron_velocity(mesh_, 0.0);
  field::CellField<double> energy_density(mesh_, 1.0);
  field::CellField<double> source(mesh_, 2.0);
  field::CellField<double> increment(mesh_, 0.0);
  ExplicitElectronEnergyStepper stepper(mesh_, electron_velocity, 0.1,
                                        uniformBoundary(1.0));
  constexpr double dt = 0.1;

  stepper.computeStableIncrement(energy_density, source, dt, increment);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {
    EXPECT_DOUBLE_EQ(energy_density[cell], 1.0);
    EXPECT_NEAR(increment[cell], 0.2, 1e-14);
  }
}

}  // namespace pemu::equation::test
