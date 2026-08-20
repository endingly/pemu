#pragma once

#include <cstdint>
#include <pemu/field/field_set.hpp>
#include <pemu/physics/ionization_reaction.hpp>
#include <pemu/physics/species.hpp>
#include <string>

namespace pemu::physics {

struct ReactionId {
  std::uint32_t value{};
  friend constexpr bool operator==(ReactionId, ReactionId) noexcept = default;
};

using ReactionRateFields = field::CellFieldSet<double, ReactionId>;

struct StoichiometricTerm {
  physics::SpeciesId species;
  double coefficient{};
};

struct Reaction {
  std::string name;
  std::vector<StoichiometricTerm> stoichiometry;
};

class ReactionNetwork {
 public:
  explicit ReactionNetwork(const physics::SpeciesSet& species)
      : species_(&species) {}

  [[nodiscard]]
  ReactionId addReaction(Reaction reaction) {
    if (reaction.name.empty()) {
      throw std::invalid_argument("reaction name must not be empty");
    }

    if (reaction.stoichiometry.empty()) {
      throw std::invalid_argument(
          "reaction stoichiometry "
          "must not be empty");
    }

    for (const auto& term : reaction.stoichiometry) {

      //
      // validates SpeciesId
      //
      auto _ = species_->at(term.species);

      if (term.coefficient == 0.0) {
        throw std::invalid_argument(
            "zero stoichiometric "
            "coefficient");
      }
    }

    const ReactionId id{static_cast<std::uint32_t>(reactions_.size())};

    reactions_.push_back(std::move(reaction));

    return id;
  }

  [[nodiscard]]
  std::size_t size() const noexcept {
    return reactions_.size();
  }

  [[nodiscard]]
  const Reaction& at(ReactionId id) const {
    return reactions_.at(id.value);
  }

  // ========================================================
  // S_s += sum_r nu_sr R_r
  //
  // reaction_rates[r] is a CellField storing R_r.
  //
  // This function ACCUMULATES into source.
  // Caller is responsible for source.fill(0) if desired.
  // ========================================================

  void accumulateSources(
      std::span<const field::CellField<double>> reaction_rates,

      physics::SpeciesCellFields& source) const {
    if (reaction_rates.size() != reactions_.size()) {

      throw std::invalid_argument(
          "reaction rate count does "
          "not match network size");
    }

    if (source.size() != species_->size()) {

      throw std::invalid_argument(
          "species source count does "
          "not match species set");
    }

    for (std::size_t r = 0; r < reactions_.size(); ++r) {

      const auto& rate = reaction_rates[r];

      if (&rate.mesh() != &source.mesh()) {

        throw std::invalid_argument(
            "reaction rate belongs "
            "to another mesh");
      }

      for (const auto& term : reactions_[r].stoichiometry) {

        auto& species_source = source[term.species];

        for (mesh::CellId cell = 0; cell < source.mesh().numCells(); ++cell) {

          species_source[cell] += term.coefficient * rate[cell];
        }
      }
    }
  }

  // ========================================================
  // Charge produced per reaction event:
  //
  //     Delta Q_r
  //
  //       = sum_s q_s nu_sr
  //
  // ========================================================

  [[nodiscard]]
  double netChargePerReaction(ReactionId id) const {
    const auto& reaction = at(id);

    double charge = 0.0;

    for (const auto& term : reaction.stoichiometry) {

      charge += species_->at(term.species).charge * term.coefficient;
    }

    return charge;
  }

  [[nodiscard]]
  bool conservesCharge(ReactionId id, double tolerance = 1e-12) const {
    return std::abs(netChargePerReaction(id)) <= tolerance;
  }

 private:
  const physics::SpeciesSet* species_;

  std::vector<Reaction> reactions_;
};

struct ElectronImpactIonizationEvaluator {
  physics::SpeciesId electron;
  physics::ReactionId ionization;
  double neutral_density{};
  double rate_coefficient{};
  void operator()(const physics::SpeciesCellFields& density,
                  const field::CellField<double>& /* potential */,
                  const field::FaceField<double>& /* electric_field */,
                  physics::ReactionRateFields& rates) const {
    physics::reaction::electronImpactIonizationRate(
        density[electron], neutral_density, rate_coefficient,
        rates[ionization]);
  }
};

}  // namespace pemu::physics
