#include <gtest/gtest.h>

#include <pemu/equation/fixed_step_multi_species_drift_diffusion_stepper.hpp>
#include <pemu/linalg/cholmod_solver.hpp>
#include <pemu/mesh/moab_mesh.hpp>
#include <pemu/physics/electron_energy.hpp>
#include <pemu/physics/reaction/network.hpp>
#include <pemu/physics/reaction/rate_coefficient.hpp>
#include <pemu/physics/wall/types.hpp>
#include <pemu/simulation/fixed_step_plasma_simulation.hpp>
#include <pemu/simulation/tabulated_electron_impact_evaluator.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

namespace pemu::simulation::test {
namespace {

/** @brief Returns the symmetric two-cell mesh used by temporal benchmarks. */
[[nodiscard]] std::filesystem::path benchmarkMeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
}

/** @brief Creates grounded potential conditions on every mesh boundary. */
[[nodiscard]] boundary::BoundaryConditionSet groundedPotentialBoundaries() {
  boundary::BoundaryConditionSet conditions;
  for (mesh::BoundaryId boundary = 1; boundary <= 4; ++boundary) {
    conditions.setDirichlet(boundary, 0.0);
  }
  return conditions;
}

/** @brief Creates empty species BC sets because wall laws own every face. */
[[nodiscard]] std::vector<boundary::BoundaryConditionSet>
wallOwnedSpeciesBoundaries(std::size_t species_count) {
  return std::vector<boundary::BoundaryConditionSet>(species_count);
}

/** @brief Returns the exterior area-to-volume ratio of one test cell. */
[[nodiscard]] double boundaryAreaToVolume(const mesh::IMesh& mesh,
                                          mesh::CellId cell) {
  double boundary_area = 0.0;
  for (const auto face : mesh.cellFaces(cell)) {
    if (mesh.isBoundary(face)) {
      boundary_area += mesh.faceArea(face);
    }
  }
  return boundary_area / mesh.cellVolume(cell);
}

/** @brief Requires every cell of a homogeneous benchmark to share one value. */
void expectHomogeneous(const field::CellField<double>& field, double expected,
                       double tolerance) {
  for (const double value : field) {
    EXPECT_NEAR(value, expected, tolerance);
  }
}

/** @brief Supplies newborn electrons with the current mean electron energy. */
struct IsoenergeticIonizationSource {
  physics::reaction::ReactionId ionization;

  /** @brief Writes w-source = mean-energy times ionization rate. */
  void operator()(const ElectronEnergySourceContext& context,
                  field::CellField<double>& source) const {
    for (mesh::CellId cell = 0; cell < source.mesh().numCells(); ++cell) {
      source[cell] =
          context.mean_energy[cell] * context.reaction_rates[ionization][cell];
    }
  }
};

/** @brief Cancels mandatory field work in the isolated wall ODE benchmark. */
struct CancelFieldPowerSource {
  /** @brief Reconstructs the same field power and writes its opposite. */
  void operator()(const ElectronEnergySourceContext& context,
                  field::CellField<double>& source) const {
    physics::computeElectronFieldPowerDensity(
        context.electron_particle_flux_normal, context.electric_field_normal,
        source);
    for (double& value : source) {
      value = -value;
    }
  }
};

struct HomogeneousIonizationResult {
  double electron_density{};
  double ion_density{};
  double energy_density{};
  double mean_energy{};
};

/** @brief Runs one source-only ionization/energy temporal refinement level. */
[[nodiscard]] HomogeneousIonizationResult runHomogeneousIonization(
    double dt, std::size_t steps) {
  constexpr double initial_density = 2.0;
  constexpr double initial_mean_energy = 3.0;
  constexpr double target_density = 4.0;
  constexpr double rate_coefficient = 0.125;

  mesh::MoabMesh mesh(benchmarkMeshPath().string());
  physics::SpeciesSet species;
  const auto electron = species.add(
      {.name = "e",
       .charge = -1.0,
       .mobility = 0.0,
       .diffusivity = 0.1,
       .transport_model = physics::SpeciesTransportModel::DriftDiffusion});
  const auto ion = species.add(
      {.name = "ion",
       .charge = +1.0,
       .mobility = 0.0,
       .diffusivity = 0.1,
       .transport_model = physics::SpeciesTransportModel::DriftDiffusion});
  physics::SpeciesCellFields density(mesh, species.size(), initial_density);
  field::CellField<double> electron_energy(
      mesh, initial_density * initial_mean_energy);
  field::CellField<double> target(mesh, target_density);

  physics::reaction::ReactionNetwork reactions(species);
  const auto ionization =
      reactions.addReaction({.name = "homogeneous isoenergetic ionization",
                             .stoichiometry = {{electron, +1.0}, {ion, +1.0}}});

  physics::wall::WallBoundarySet zero_flux_walls;
  for (mesh::BoundaryId boundary = 1; boundary <= 4; ++boundary) {
    zero_flux_walls.set(boundary,
                        {.losses = {{electron, 0.0, 0.0}, {ion, 0.0, 0.0}}});
  }
  equation::FixedStepMultiSpeciesDriftDiffusionStepper transport(
      mesh, species, 1.0, dt, groundedPotentialBoundaries(),
      wallOwnedSpeciesBoundaries(species.size()),
      std::make_unique<linalg::CholmodSolver>(), {}, std::nullopt,
      std::move(zero_flux_walls));
  TabulatedElectronImpactEvaluator rate_evaluator(
      electron, ionization, target,
      ElectronImpactRateCoordinate::electron_temperature_ev,
      physics::reaction::TabulatedRateCoefficient(
          {0.0, 4.0}, {rate_coefficient, rate_coefficient},
          physics::reaction::RateInterpolation::linear),
      true);
  ElectronEnergyConfiguration energy_configuration{
      .energy_density = electron_energy,
      .electron = electron,
      .boundary_conditions = {},
      .additional_source_evaluator = IsoenergeticIonizationSource{ionization},
      .density_floor = 0.0,
      .allow_unitless_raw_values = true};
  FixedStepPlasmaSimulation simulation(
      density, reactions, transport, std::move(rate_evaluator),
      std::move(energy_configuration), FixedStepClock(dt, steps));

  simulation.run();

  EXPECT_TRUE(simulation.finished());
  EXPECT_EQ(simulation.step(), steps);
  return {.electron_density = density[electron][mesh::CellId{0}],
          .ion_density = density[ion][mesh::CellId{0}],
          .energy_density = electron_energy[mesh::CellId{0}],
          .mean_energy = simulation.electronMeanEnergy()[mesh::CellId{0}]};
}

struct WallEmissionResult {
  double electron_density{};
  double ion_density{};
  double energy_density{};
};

/** @brief Runs one wall-loss/secondary-emission temporal refinement level. */
[[nodiscard]] WallEmissionResult runWallEmission(double dt, std::size_t steps) {
  constexpr double initial_electron_density = 5.0;
  constexpr double initial_ion_density = 7.0;
  constexpr double initial_energy_density = 11.0;
  constexpr double electron_loss_velocity = 0.2;
  constexpr double ion_loss_velocity = 0.1;
  constexpr double energy_loss_velocity = 0.3;
  constexpr double secondary_yield = 0.25;
  constexpr double emitted_mean_energy = 2.0;

  mesh::MoabMesh mesh(benchmarkMeshPath().string());
  physics::SpeciesSet species;
  const auto electron = species.add(
      {.name = "e",
       .charge = -1.0,
       .mobility = 0.0,
       .diffusivity = 0.1,
       .transport_model = physics::SpeciesTransportModel::DriftDiffusion});
  const auto ion = species.add(
      {.name = "ion",
       .charge = +1.0,
       .mobility = 0.0,
       .diffusivity = 0.1,
       .transport_model = physics::SpeciesTransportModel::DriftDiffusion});
  physics::SpeciesCellFields density(mesh, species.size(), 0.0);
  density[electron].fill(initial_electron_density);
  density[ion].fill(initial_ion_density);
  field::CellField<double> electron_energy(mesh, initial_energy_density);
  physics::reaction::ReactionNetwork reactions(species);

  physics::wall::WallBoundarySet walls;
  for (mesh::BoundaryId boundary = 1; boundary <= 4; ++boundary) {
    walls.set(
        boundary,
        {.losses = {{electron, electron_loss_velocity, energy_loss_velocity},
                    {ion, ion_loss_velocity, 0.0}},
         .secondary_emissions = {
             {ion, electron, secondary_yield, emitted_mean_energy}}});
  }
  equation::FixedStepMultiSpeciesDriftDiffusionStepper transport(
      mesh, species, 1.0, dt, groundedPotentialBoundaries(),
      wallOwnedSpeciesBoundaries(species.size()),
      std::make_unique<linalg::CholmodSolver>(), {}, std::nullopt,
      std::move(walls));
  ElectronEnergyConfiguration energy_configuration{
      .energy_density = electron_energy,
      .electron = electron,
      .boundary_conditions = {},
      .additional_source_evaluator = CancelFieldPowerSource{},
      .density_floor = 0.0,
      .allow_unitless_raw_values = true};
  FixedStepPlasmaSimulation simulation(
      density, reactions, transport,
      [](const PlasmaReactionRateContext&,
         physics::reaction::ReactionRateFields&) {},
      std::move(energy_configuration), FixedStepClock(dt, steps));

  simulation.run();

  EXPECT_TRUE(simulation.finished());
  expectHomogeneous(density[electron], density[electron][mesh::CellId{0}],
                    1e-12);
  expectHomogeneous(density[ion], density[ion][mesh::CellId{0}], 1e-12);
  expectHomogeneous(electron_energy, electron_energy[mesh::CellId{0}], 1e-12);
  return {.electron_density = density[electron][mesh::CellId{0}],
          .ion_density = density[ion][mesh::CellId{0}],
          .energy_density = electron_energy[mesh::CellId{0}]};
}

/** @brief Returns the observed order between two halved timestep errors. */
[[nodiscard]] double observedOrder(double coarse_error, double fine_error) {
  return std::log(coarse_error / fine_error) / std::log(2.0);
}

}  // namespace

TEST(PlasmaBenchmarkTest,
     HomogeneousIonizationEnergyHasFirstOrderTemporalConvergence) {
  constexpr double initial_density = 2.0;
  constexpr double initial_mean_energy = 3.0;
  constexpr double growth_rate = 0.5;
  constexpr double end_time = 0.4;
  constexpr std::array time_steps{0.1, 0.05, 0.025};
  std::array<double, time_steps.size()> density_errors{};

  for (std::size_t level = 0; level < time_steps.size(); ++level) {
    const double dt = time_steps[level];
    const auto steps = static_cast<std::size_t>(std::llround(end_time / dt));
    const auto result = runHomogeneousIonization(dt, steps);
    const double discrete_density =
        initial_density * std::pow(1.0 + growth_rate * dt, steps);
    const double continuous_density =
        initial_density * std::exp(growth_rate * end_time);

    EXPECT_NEAR(result.electron_density, discrete_density, 1e-12);
    EXPECT_NEAR(result.ion_density, discrete_density, 1e-12);
    EXPECT_NEAR(result.energy_density, initial_mean_energy * discrete_density,
                1e-11);
    EXPECT_NEAR(result.mean_energy, initial_mean_energy, 1e-12);
    density_errors[level] =
        std::abs(result.electron_density - continuous_density);
  }

  EXPECT_GT(observedOrder(density_errors[0], density_errors[1]), 0.9);
  EXPECT_GT(observedOrder(density_errors[1], density_errors[2]), 0.9);
}

TEST(PlasmaBenchmarkTest,
     WallLossAndSecondaryEmissionHaveFirstOrderTemporalConvergence) {
  constexpr double initial_electron_density = 5.0;
  constexpr double initial_ion_density = 7.0;
  constexpr double initial_energy_density = 11.0;
  constexpr double secondary_yield = 0.25;
  constexpr double emitted_mean_energy = 2.0;
  constexpr double end_time = 0.4;
  constexpr std::array time_steps{0.05, 0.025, 0.0125};

  mesh::MoabMesh mesh(benchmarkMeshPath().string());
  const double area_to_volume = boundaryAreaToVolume(mesh, mesh::CellId{0});
  EXPECT_DOUBLE_EQ(area_to_volume, boundaryAreaToVolume(mesh, mesh::CellId{1}));
  const double electron_loss_rate = 0.2 * area_to_volume;
  const double ion_loss_rate = 0.1 * area_to_volume;
  const double energy_loss_rate = 0.3 * area_to_volume;
  const double electron_inflow_rate = secondary_yield * ion_loss_rate;
  const double energy_inflow_rate = emitted_mean_energy * electron_inflow_rate;

  const double exact_ion =
      initial_ion_density * std::exp(-ion_loss_rate * end_time);
  const double exact_electron =
      initial_electron_density * std::exp(-electron_loss_rate * end_time) +
      electron_inflow_rate * initial_ion_density *
          (std::exp(-ion_loss_rate * end_time) -
           std::exp(-electron_loss_rate * end_time)) /
          (electron_loss_rate - ion_loss_rate);
  const double exact_energy =
      initial_energy_density * std::exp(-energy_loss_rate * end_time) +
      energy_inflow_rate * initial_ion_density *
          (std::exp(-ion_loss_rate * end_time) -
           std::exp(-energy_loss_rate * end_time)) /
          (energy_loss_rate - ion_loss_rate);
  std::array<double, time_steps.size()> combined_errors{};

  for (std::size_t level = 0; level < time_steps.size(); ++level) {
    const double dt = time_steps[level];
    const auto steps = static_cast<std::size_t>(std::llround(end_time / dt));
    const auto result = runWallEmission(dt, steps);
    double discrete_electron = initial_electron_density;
    double discrete_ion = initial_ion_density;
    double discrete_energy = initial_energy_density;
    for (std::size_t step = 0; step < steps; ++step) {
      const double previous_ion = discrete_ion;
      discrete_ion -= dt * ion_loss_rate * previous_ion;
      discrete_electron += dt * (-electron_loss_rate * discrete_electron +
                                 electron_inflow_rate * previous_ion);
      discrete_energy += dt * (-energy_loss_rate * discrete_energy +
                               energy_inflow_rate * previous_ion);
    }

    EXPECT_NEAR(result.electron_density, discrete_electron, 1e-12);
    EXPECT_NEAR(result.ion_density, discrete_ion, 1e-12);
    EXPECT_NEAR(result.energy_density, discrete_energy, 1e-12);
    combined_errors[level] =
        std::max({std::abs(result.electron_density - exact_electron),
                  std::abs(result.ion_density - exact_ion),
                  std::abs(result.energy_density - exact_energy)});
  }

  EXPECT_GT(observedOrder(combined_errors[0], combined_errors[1]), 0.85);
  EXPECT_GT(observedOrder(combined_errors[1], combined_errors[2]), 0.85);
}

}  // namespace pemu::simulation::test
