#pragma once

#include <pemu/field/field_set.hpp>
#include <pemu/physics/species.hpp>

#include <cstdint>
#include <string>
#include <vector>

namespace pemu::physics::reaction {

/** @brief Dense stable identifier for one reaction in a network. */
struct ReactionId {
  std::uint32_t value{};
  friend constexpr bool operator==(ReactionId, ReactionId) noexcept = default;
};

/** @brief Homogeneous cell fields indexed by ReactionId. */
using ReactionRateFields = field::CellFieldSet<double, ReactionId>;

/** @brief Net species production per occurrence of one reaction. */
struct StoichiometricTerm {
  physics::SpeciesId species;
  double coefficient{};
};

/**
 * @brief Species exponent in a mass-action rate law.
 *
 * Kinetic order is deliberately independent of net stoichiometry. A catalyst
 * or third body can have non-zero kinetic order while having zero net source.
 */
struct KineticOrderTerm {
  physics::SpeciesId species;
  double order{};
};

/** @brief Structural and kinetic definition of one chemical reaction. */
struct Reaction {
  std::string name;
  std::vector<StoichiometricTerm> stoichiometry;
  std::vector<KineticOrderTerm> kinetic_orders;
};

}  // namespace pemu::physics::reaction
