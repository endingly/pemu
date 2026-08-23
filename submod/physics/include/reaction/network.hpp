#pragma once

#include <pemu/physics/reaction/types.hpp>

#include <algorithm>
#include <cmath>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pemu::physics::reaction {

/** @brief Owns reaction definitions and maps reaction rates to species sources. */
class ReactionNetwork {
 public:
  /** @brief Binds the species registry used by every reaction definition. */
  explicit ReactionNetwork(const physics::SpeciesSet& species)
      : species_(&species) {}

  /** @brief Validates and appends one reaction, returning its dense ID. */
  [[nodiscard]] ReactionId addReaction(Reaction reaction) {
    validateReaction(reaction);
    const ReactionId id{static_cast<std::uint32_t>(reactions_.size())};
    reactions_.push_back(std::move(reaction));
    return id;
  }

  /** @brief Returns the number of reactions in insertion order. */
  [[nodiscard]] std::size_t size() const noexcept { return reactions_.size(); }

  /** @brief Returns one validated reaction definition. */
  [[nodiscard]] const Reaction& at(ReactionId id) const {
    return reactions_.at(id.value);
  }

  /** @brief Returns the species registry bound to this network. */
  [[nodiscard]] const physics::SpeciesSet& species() const noexcept {
    return *species_;
  }

  /**
   * @brief Accumulates S_s += sum_r nu_sr R_r into existing species sources.
   */
  void accumulateSources(
      std::span<const field::CellField<double>> reaction_rates,
      physics::SpeciesCellFields& source) const {
    if (reaction_rates.size() != reactions_.size()) {
      throw std::invalid_argument(
          "reaction rate count does not match network size");
    }
    if (source.size() != species_->size()) {
      throw std::invalid_argument(
          "species source count does not match species set");
    }

    for (std::size_t reaction_index = 0; reaction_index < reactions_.size();
         ++reaction_index) {
      const auto& rate = reaction_rates[reaction_index];
      if (&rate.mesh() != &source.mesh()) {
        throw std::invalid_argument("reaction rate belongs to another mesh");
      }
      for (const double value : rate) {
        if (!std::isfinite(value) || value < 0.0) {
          throw std::invalid_argument(
              "reaction rates must be finite and non-negative");
        }
      }
    }

    for (std::size_t reaction_index = 0; reaction_index < reactions_.size();
         ++reaction_index) {
      const auto& rate = reaction_rates[reaction_index];
      for (const auto& term : reactions_[reaction_index].stoichiometry) {
        auto& species_source = source[term.species];
        for (mesh::CellId cell = 0; cell < source.mesh().numCells(); ++cell) {
          species_source[cell] += term.coefficient * rate[cell];
        }
      }
    }
  }

  /** @brief Returns the net electric charge created by one reaction event. */
  [[nodiscard]] double netChargePerReaction(ReactionId id) const {
    double charge = 0.0;
    for (const auto& term : at(id).stoichiometry) {
      charge += species_->at(term.species).charge * term.coefficient;
    }
    return charge;
  }

  /** @brief Reports charge conservation within an absolute tolerance. */
  [[nodiscard]] bool conservesCharge(ReactionId id,
                                     double tolerance = 1e-12) const {
    if (!std::isfinite(tolerance) || tolerance < 0.0) {
      throw std::invalid_argument(
          "charge-conservation tolerance must be finite and non-negative");
    }
    return std::abs(netChargePerReaction(id)) <= tolerance;
  }

 private:
  /** @brief Validates names, species IDs, coefficients, orders, and duplicates. */
  void validateReaction(const Reaction& reaction) const {
    if (reaction.name.empty()) {
      throw std::invalid_argument("reaction name must not be empty");
    }
    if (reaction.stoichiometry.empty()) {
      throw std::invalid_argument("reaction stoichiometry must not be empty");
    }

    std::vector<physics::SpeciesId> stoichiometric_species;
    stoichiometric_species.reserve(reaction.stoichiometry.size());
    for (const auto& term : reaction.stoichiometry) {
      static_cast<void>(species_->at(term.species));
      if (!std::isfinite(term.coefficient) || term.coefficient == 0.0) {
        throw std::invalid_argument(
            "stoichiometric coefficient must be finite and non-zero");
      }
      if (std::ranges::find(stoichiometric_species, term.species) !=
          stoichiometric_species.end()) {
        throw std::invalid_argument(
            "reaction contains duplicate stoichiometric species");
      }
      stoichiometric_species.push_back(term.species);
    }

    std::vector<physics::SpeciesId> kinetic_species;
    kinetic_species.reserve(reaction.kinetic_orders.size());
    for (const auto& term : reaction.kinetic_orders) {
      static_cast<void>(species_->at(term.species));
      if (!std::isfinite(term.order) || term.order <= 0.0) {
        throw std::invalid_argument(
            "kinetic order must be finite and positive");
      }
      if (std::ranges::find(kinetic_species, term.species) !=
          kinetic_species.end()) {
        throw std::invalid_argument(
            "reaction contains duplicate kinetic-order species");
      }
      kinetic_species.push_back(term.species);
    }
  }

  const physics::SpeciesSet* species_;
  std::vector<Reaction> reactions_;
};

}  // namespace pemu::physics::reaction
