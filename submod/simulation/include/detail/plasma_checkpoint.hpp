#pragma once

#include <pemu/output/checkpoint/types.hpp>
#include <pemu/physics/species.hpp>

#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace pemu::simulation::detail {

inline constexpr std::string_view plasma_checkpoint_schema_key =
    "plasma_workflow_schema_version";
inline constexpr std::string_view plasma_checkpoint_kind_key =
    "plasma_workflow_kind";
inline constexpr double plasma_checkpoint_schema_version = 2.0;
inline constexpr double fixed_step_checkpoint_kind = 1.0;
inline constexpr double adaptive_step_checkpoint_kind = 2.0;

[[nodiscard]] inline std::string checkpointHexToken(std::string_view value) {
  constexpr char digits[] = "0123456789abcdef";
  std::string encoded;
  encoded.reserve(value.size() * 2);
  for (const unsigned char character : value) {
    encoded.push_back(digits[character >> 4]);
    encoded.push_back(digits[character & 0x0f]);
  }
  return encoded;
}

[[nodiscard]] inline std::string speciesDensityCheckpointKey(
    physics::SpeciesId id, std::string_view species_name) {
  return "species_" + std::to_string(id.value) + "_name_" +
         checkpointHexToken(species_name) + "_number_density";
}

/** @brief Builds the restart key that binds energy to one electron species. */
[[nodiscard]] inline std::string electronEnergyDensityCheckpointKey(
    physics::SpeciesId id, std::string_view species_name) {
  return "electron_species_" + std::to_string(id.value) + "_name_" +
         checkpointHexToken(species_name) + "_energy_density";
}

[[nodiscard]] inline std::vector<output::checkpoint::CellFieldSource>
plasmaCheckpointSources(const physics::SpeciesCellFields& density,
                        const physics::SpeciesSet& species) {
  if (density.size() != species.size()) {
    throw std::invalid_argument(
        "checkpoint density and species counts must match");
  }
  std::vector<output::checkpoint::CellFieldSource> sources;
  sources.reserve(density.size());
  for (std::size_t index = 0; index < density.size(); ++index) {
    const physics::SpeciesId id{static_cast<std::uint32_t>(index)};
    sources.push_back(
        {.field = &density[id],
         .key = speciesDensityCheckpointKey(id, species.at(id).name)});
  }
  return sources;
}

[[nodiscard]] inline std::vector<output::checkpoint::CellFieldTarget>
plasmaCheckpointTargets(physics::SpeciesCellFields& density,
                        const physics::SpeciesSet& species) {
  if (density.size() != species.size()) {
    throw std::invalid_argument(
        "checkpoint density and species counts must match");
  }
  std::vector<output::checkpoint::CellFieldTarget> targets;
  targets.reserve(density.size());
  for (std::size_t index = 0; index < density.size(); ++index) {
    const physics::SpeciesId id{static_cast<std::uint32_t>(index)};
    targets.push_back(
        {.field = &density[id],
         .key = speciesDensityCheckpointKey(id, species.at(id).name)});
  }
  return targets;
}

}  // namespace pemu::simulation::detail
