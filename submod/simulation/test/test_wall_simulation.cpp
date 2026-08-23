#include <gtest/gtest.h>

#include <pemu/equation/adaptive_step_multi_species_drift_diffusion_stepper.hpp>
#include <pemu/equation/fixed_step_multi_species_drift_diffusion_stepper.hpp>
#include <pemu/linalg/cholmod_solver.hpp>
#include <pemu/mesh/moab_mesh.hpp>
#include <pemu/physics/reaction/network.hpp>
#include <pemu/physics/wall/types.hpp>
#include <pemu/simulation/adaptive_step_plasma_simulation.hpp>
#include <pemu/simulation/fixed_step_plasma_simulation.hpp>

#include <filesystem>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pemu::simulation::test {
namespace {

/** @brief Returns the two-cell end-to-end wall test mesh. */
[[nodiscard]] std::filesystem::path testMeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
}

/** @brief Finds one face by physical boundary ID. */
[[nodiscard]] mesh::FaceId findBoundaryFace(const mesh::IMesh& mesh,
                                            mesh::BoundaryId boundary) {
  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {
    if (mesh.isBoundary(face) && mesh.boundaryId(face) == boundary) {
      return face;
    }
  }
  throw std::runtime_error("test boundary face not found");
}

/** @brief Creates grounded electrostatic boundaries. */
[[nodiscard]] boundary::BoundaryConditionSet potentialBoundaries() {
  boundary::BoundaryConditionSet conditions;
  for (mesh::BoundaryId boundary = 1; boundary <= 4; ++boundary) {
    conditions.setDirichlet(boundary, 0.0);
  }
  return conditions;
}

/** @brief Creates matching SG boundaries outside the left plasma wall. */
[[nodiscard]] std::vector<boundary::BoundaryConditionSet> speciesBoundaries(
    std::size_t count, double value) {
  std::vector<boundary::BoundaryConditionSet> conditions(count);
  for (auto& species_conditions : conditions) {
    species_conditions.setDirichlet(mesh::BoundaryId{2}, value);
    species_conditions.setDirichlet(mesh::BoundaryId{3}, value);
    species_conditions.setDirichlet(mesh::BoundaryId{4}, value);
  }
  return conditions;
}

/** @brief Creates matching energy boundaries outside the left plasma wall. */
[[nodiscard]] boundary::BoundaryConditionSet energyBoundaries(double value) {
  boundary::BoundaryConditionSet conditions;
  conditions.setDirichlet(mesh::BoundaryId{2}, value);
  conditions.setDirichlet(mesh::BoundaryId{3}, value);
  conditions.setDirichlet(mesh::BoundaryId{4}, value);
  return conditions;
}

/** @brief Integrates one cell field over physical cell volumes. */
[[nodiscard]] double integratedField(const mesh::IMesh& mesh,
                                     const field::CellField<double>& values) {
  double result = 0.0;
  for (mesh::CellId cell = 0; cell < mesh.numCells(); ++cell) {
    result += values[cell] * mesh.cellVolume(cell);
  }
  return result;
}

}  // namespace

TEST(FixedStepWallSimulationTest,
     AdvancesParticleAndElectronEnergyWallFluxInOneAtomicWorkflow) {
  mesh::MoabMesh mesh(testMeshPath().string());
  physics::SpeciesSet species;
  const auto electron = species.add(
      {.name = "e",
       .charge = -1.0,
       .mobility = 1.0,
       .diffusivity = 0.1,
       .transport_model = physics::SpeciesTransportModel::DriftDiffusion});
  const auto ion = species.add(
      {.name = "ion",
       .charge = +1.0,
       .mobility = 0.5,
       .diffusivity = 0.1,
       .transport_model = physics::SpeciesTransportModel::DriftDiffusion});
  constexpr double density_value = 7.0;
  constexpr double energy_density_value = 14.0;
  constexpr double dt = 1e-4;
  physics::SpeciesCellFields density(mesh, species.size(), density_value);
  field::CellField<double> electron_energy(mesh, energy_density_value);
  physics::reaction::ReactionNetwork reactions(species);
  physics::wall::WallBoundarySet wall_boundaries;
  wall_boundaries.set(mesh::BoundaryId{1},
                      {.losses = {{electron, 2.0, 5.0}, {ion, 3.0, 0.0}},
                       .secondary_emissions = {{ion, electron, 0.1, 4.0}}});
  equation::FixedStepMultiSpeciesDriftDiffusionStepper transport(
      mesh, species, 1.0, dt, potentialBoundaries(),
      speciesBoundaries(species.size(), density_value),
      std::make_unique<linalg::CholmodSolver>(), {}, std::nullopt,
      std::move(wall_boundaries));
  ElectronEnergyConfiguration energy_configuration{
      .energy_density = electron_energy,
      .electron = electron,
      .boundary_conditions = energyBoundaries(energy_density_value),
      .additional_source_evaluator =
          [](const ElectronEnergySourceContext&,
             field::CellField<double>& source) { source.fill(0.0); },
      .density_floor = 0.0,
      .allow_unitless_raw_values = true};
  FixedStepPlasmaSimulation simulation(
      density, reactions, transport,
      [](const PlasmaReactionRateContext&,
         physics::reaction::ReactionRateFields&) {},
      std::move(energy_configuration), FixedStepClock(dt, 1));
  const auto wall_face = findBoundaryFace(mesh, mesh::BoundaryId{1});
  const double electron_before = integratedField(mesh, density[electron]);
  const double ion_before = integratedField(mesh, density[ion]);
  const double energy_before = integratedField(mesh, electron_energy);

  const auto result = simulation.advance();

  ASSERT_TRUE(result.success());
  const double area = mesh.faceArea(wall_face);
  EXPECT_NEAR(integratedField(mesh, density[electron]) - electron_before,
              -dt * 11.9 * area, 1e-12);
  EXPECT_NEAR(integratedField(mesh, density[ion]) - ion_before,
              -dt * 21.0 * area, 1e-12);
  EXPECT_NEAR(integratedField(mesh, electron_energy) - energy_before,
              -dt * (5.0 * energy_density_value - 8.4) * area, 1e-12);
  EXPECT_TRUE(simulation.finished());
}

TEST(AdaptiveStepWallSimulationTest,
     AdvancesParticleAndElectronEnergyWallFluxInOneAtomicWorkflow) {
  mesh::MoabMesh mesh(testMeshPath().string());
  physics::SpeciesSet species;
  const auto electron = species.add(
      {.name = "e",
       .charge = -1.0,
       .mobility = 1.0,
       .diffusivity = 0.1,
       .transport_model = physics::SpeciesTransportModel::DriftDiffusion});
  const auto ion = species.add(
      {.name = "ion",
       .charge = +1.0,
       .mobility = 0.5,
       .diffusivity = 0.1,
       .transport_model = physics::SpeciesTransportModel::DriftDiffusion});
  constexpr double density_value = 7.0;
  constexpr double energy_density_value = 14.0;
  constexpr double end_time = 1e-4;
  physics::SpeciesCellFields density(mesh, species.size(), density_value);
  field::CellField<double> electron_energy(mesh, energy_density_value);
  physics::reaction::ReactionNetwork reactions(species);
  physics::wall::WallBoundarySet wall_boundaries;
  wall_boundaries.set(mesh::BoundaryId{1},
                      {.losses = {{electron, 2.0, 5.0}, {ion, 3.0, 0.0}},
                       .secondary_emissions = {{ion, electron, 0.1, 4.0}}});
  equation::AdaptiveStepMultiSpeciesDriftDiffusionStepper transport(
      mesh, species, 1.0, potentialBoundaries(),
      speciesBoundaries(species.size(), density_value),
      std::make_unique<linalg::CholmodSolver>(),
      {.safety = 1.0, .min_dt = 1e-8, .max_dt = end_time, .max_growth = 2.0},
      {}, std::nullopt, std::move(wall_boundaries));
  ElectronEnergyConfiguration energy_configuration{
      .energy_density = electron_energy,
      .electron = electron,
      .boundary_conditions = energyBoundaries(energy_density_value),
      .additional_source_evaluator =
          [](const ElectronEnergySourceContext&,
             field::CellField<double>& source) { source.fill(0.0); },
      .density_floor = 0.0,
      .allow_unitless_raw_values = true};
  AdaptiveStepPlasmaSimulation simulation(
      density, reactions, transport,
      [](const PlasmaReactionRateContext&,
         physics::reaction::ReactionRateFields&) {},
      std::move(energy_configuration), AdaptiveTimeClock(end_time));
  const auto wall_face = findBoundaryFace(mesh, mesh::BoundaryId{1});
  const double electron_before = integratedField(mesh, density[electron]);
  const double ion_before = integratedField(mesh, density[ion]);
  const double energy_before = integratedField(mesh, electron_energy);

  const auto result = simulation.advance();

  ASSERT_TRUE(result.success());
  ASSERT_DOUBLE_EQ(simulation.lastTimeStep(), end_time);
  const double area = mesh.faceArea(wall_face);
  EXPECT_NEAR(integratedField(mesh, density[electron]) - electron_before,
              -end_time * 11.9 * area, 1e-12);
  EXPECT_NEAR(integratedField(mesh, density[ion]) - ion_before,
              -end_time * 21.0 * area, 1e-12);
  EXPECT_NEAR(integratedField(mesh, electron_energy) - energy_before,
              -end_time * (5.0 * energy_density_value - 8.4) * area, 1e-12);
  EXPECT_TRUE(simulation.finished());
}

}  // namespace pemu::simulation::test
