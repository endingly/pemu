#pragma once

#include <pemu/output/dump/types.hpp>
#include <pemu/physics/species.hpp>
#include <pemu/simulation/field_output.hpp>

#include <cstdint>
#include <filesystem>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace pemu::simulation::detail {

/**
 * @brief Determines whether a configured output schedule selects a state.
 * @param options Field-output configuration.
 * @param step Current completed-state step index.
 * @param terminal Whether the state is the simulation's terminal state.
 * @return `true` when the state should be written.
 */
[[nodiscard]]
inline bool shouldWritePlasmaFieldSnapshot(const FieldOutputOptions& options,
                                           std::size_t step,
                                           bool terminal) noexcept {
  if (!options.enabled()) {
    return false;
  }
  if (terminal) {
    return options.write_final ||
           (options.every_steps != 0 && step % options.every_steps == 0);
  }
  if (step == 0) {
    return options.write_initial;
  }
  return options.every_steps != 0 && step % options.every_steps == 0;
}

/**
 * @brief Builds the deterministic path for one plasma time series.
 * @param options Field-output configuration.
 * @return Output path below the configured directory.
 */
[[nodiscard]]
inline std::filesystem::path plasmaFieldOutputPath(
    const FieldOutputOptions& options) {
  return options.directory / (options.file_stem + ".vtkhdf");
}

/**
 * @brief Writes the synchronized density and electrostatic fields of a plasma
 * state.
 *
 * Species fields receive stable, species-ID-prefixed names so homogeneous
 * number-density metadata cannot cause duplicate VTK array names.
 *
 * @tparam Stepper A multi-species electrostatic transport stepper.
 * @param density Species densities synchronized with the stepper fields.
 * @param stepper Stepper supplying charge, potential, electric, and drift
 * fields.
 * @param options Enabled output configuration.
 * @param series Open series receiving the state.
 * @param stamp Step and time associated with the state.
 * @return Record returned by the configured writer.
 * @throws std::logic_error if output is disabled.
 */
template <typename Stepper>
[[nodiscard]]
output::OutputRecord writePlasmaFieldSnapshot(
    const physics::SpeciesCellFields& density, const Stepper& stepper,
    const FieldOutputOptions& options, output::dump::ISeries& series,
    output::OutputStamp stamp) {
  if (!options.enabled()) {
    throw std::logic_error(
        "cannot write a snapshot with field output disabled");
  }

  std::vector<output::dump::CellFieldSelection> cell_fields;
  cell_fields.reserve(density.size() + 2);
  std::vector<output::dump::FaceFieldSelection> face_fields;
  face_fields.reserve(density.size() + 1);

  for (std::size_t index = 0; index < density.size(); ++index) {
    const physics::SpeciesId id{static_cast<std::uint32_t>(index)};
    const auto& species = stepper.species().at(id);
    const std::string prefix =
        "species_" + std::to_string(id.value) + "_" + species.name;
    cell_fields.push_back(
        {.field = &density[id], .name = prefix + "_number_density"});
    face_fields.push_back(
        {.field = &stepper.driftVelocityNormal(id),
         .name = prefix + "_normal_drift_velocity",
         .include_cell_centered_visualization =
             options.include_face_cell_visualization,
         .cell_centered_name =
             prefix + "_normal_drift_velocity_cell_centered"});
  }

  cell_fields.push_back(
      {.field = &stepper.chargeDensity(), .name = "charge_density"});
  cell_fields.push_back(
      {.field = &stepper.potential(), .name = "electric_potential"});
  face_fields.push_back(
      {.field = &stepper.electricFieldNormal(),
       .name = "normal_electric_field",
       .include_cell_centered_visualization =
           options.include_face_cell_visualization,
       .cell_centered_name = "normal_electric_field_cell_centered"});

  return series.append(
      {.mesh = &stepper.mesh(),
       .path = plasmaFieldOutputPath(options),
       .stamp = stamp,
       .cell_field_selections =
           std::span<const output::dump::CellFieldSelection>{cell_fields},
       .face_fields =
           std::span<const output::dump::FaceFieldSelection>{face_fields},
       .overwrite = options.overwrite});
}

}  // namespace pemu::simulation::detail
