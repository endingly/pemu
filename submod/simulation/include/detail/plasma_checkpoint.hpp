#pragma once

#include <pemu/output/checkpoint/types.hpp>
#include <pemu/physics/species.hpp>

#include <cstdint>
#include <string>
#include <vector>

namespace pemu::simulation::detail {

[[nodiscard]] inline std::string speciesDensityCheckpointKey(
    physics::SpeciesId id) {
  return "species_" + std::to_string(id.value) + "_number_density";
}

[[nodiscard]] inline std::vector<output::checkpoint::CellFieldSource>
plasmaCheckpointSources(const physics::SpeciesCellFields& density) {
  std::vector<output::checkpoint::CellFieldSource> sources;
  sources.reserve(density.size());
  for (std::size_t index = 0; index < density.size(); ++index) {
    const physics::SpeciesId id{static_cast<std::uint32_t>(index)};
    sources.push_back(
        {.field = &density[id], .key = speciesDensityCheckpointKey(id)});
  }
  return sources;
}

[[nodiscard]] inline std::vector<output::checkpoint::CellFieldTarget>
plasmaCheckpointTargets(physics::SpeciesCellFields& density) {
  std::vector<output::checkpoint::CellFieldTarget> targets;
  targets.reserve(density.size());
  for (std::size_t index = 0; index < density.size(); ++index) {
    const physics::SpeciesId id{static_cast<std::uint32_t>(index)};
    targets.push_back(
        {.field = &density[id], .key = speciesDensityCheckpointKey(id)});
  }
  return targets;
}

}  // namespace pemu::simulation::detail
