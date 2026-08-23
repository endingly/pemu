#pragma once

#include <pemu/boundary/boundary_condition_set.hpp>
#include <pemu/equation/explicit_electron_energy_stepper.hpp>
#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>
#include <pemu/physics/reaction.hpp>
#include <pemu/physics/species.hpp>

#include <functional>

namespace pemu::simulation {

/** @brief Read-only plasma state available while assembling the energy source. */
struct ElectronEnergySourceContext {
  const physics::SpeciesCellFields& density;
  const field::CellField<double>& potential;
  const field::FaceField<double>& electric_field_normal;
  const field::FaceField<double>& electron_particle_flux_normal;
  const physics::ReactionRateFields& reaction_rates;
  const physics::SpeciesCellFields& species_source;
  const field::CellField<double>& mean_energy;
};

/**
 * @brief Builds additional collisional or externally imposed energy sources.
 *
 * Simulation assembles electric-field power from the same SG particle flux as
 * the electron continuity equation. This callback must overwrite every target
 * cell with only the remaining source terms, such as collisional exchange. It
 * is evaluated after electrostatics, reaction rates, and species sources.
 */
using ElectronEnergyAdditionalSourceEvaluator = std::function<void(
    const ElectronEnergySourceContext&, field::CellField<double>&)>;

/**
 * @brief Mandatory state and closure configuration for electron energy.
 *
 * Simulation advances the referenced energy-density field in place. The
 * caller retains ownership, matching the ownership model of species density.
 * Physical metadata is required by default. Unitless synthetic workflows must
 * explicitly opt into the raw-value convention, where potential is measured in
 * volts and electron energy in electronvolts.
 */
struct ElectronEnergyConfiguration {
  field::CellField<double>& energy_density;
  physics::SpeciesId electron;
  boundary::BoundaryConditionSet boundary_conditions;
  ElectronEnergyAdditionalSourceEvaluator additional_source_evaluator;
  double density_floor{};
  double transport_factor{
      equation::ExplicitElectronEnergyStepper::maxwellian_transport_factor};
  bool allow_unitless_raw_values{false};
};

}  // namespace pemu::simulation
