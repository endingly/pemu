#pragma once

#include <pemu/physics/reaction/mass_action.hpp>
#include <pemu/physics/reaction/network.hpp>
#include <pemu/unit/quantity_metadata.hpp>

#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace pemu::physics::reaction {

/**
 * @brief Pre-resolves one ReactionNetwork into reusable mass-action factors.
 *
 * The assembler allocates only during construction. Each evaluate() call uses
 * the current referenced density and coefficient values, supports any number
 * of kinetic factors, and writes one rate field per reaction.
 */
class MassActionReactionAssembler {
 public:
  /** @brief Binds a network to stable species-density fields. */
  MassActionReactionAssembler(const ReactionNetwork& network,
                              const physics::SpeciesCellFields& density,
                              bool allow_unitless_raw_values = false)
      : network_(&network),
        density_(&density),
        allow_unitless_raw_values_(allow_unitless_raw_values) {
    if (density.size() != network.species().size()) {
      throw std::invalid_argument(
          "mass-action density and species counts differ");
    }
    resolved_terms_.reserve(network.size());
    for (std::size_t reaction_index = 0; reaction_index < network.size();
         ++reaction_index) {
      const auto& definition =
          network.at(ReactionId{static_cast<std::uint32_t>(reaction_index)});
      if (definition.kinetic_orders.empty()) {
        throw std::invalid_argument(
            "mass-action reaction requires explicit kinetic orders");
      }
      auto& terms = resolved_terms_.emplace_back();
      terms.reserve(definition.kinetic_orders.size());
      for (const auto& kinetic : definition.kinetic_orders) {
        terms.emplace_back(density[kinetic.species], kinetic.order);
      }
    }
  }

  /** @brief Evaluates all configured reactions from current field values. */
  void evaluate(std::span<const field::CellField<double>> rate_coefficients,
                ReactionRateFields& reaction_rates) const {
    if (rate_coefficients.size() != network_->size() ||
        reaction_rates.size() != network_->size()) {
      throw std::invalid_argument(
          "mass-action coefficient or rate count differs from network");
    }
    if (&reaction_rates.mesh() != &density_->mesh()) {
      throw std::invalid_argument("mass-action rates belong to another mesh");
    }
    validateMetadata(rate_coefficients, reaction_rates);

    // Validate the complete batch first so one invalid later reaction cannot
    // leave earlier rate fields partially updated.
    for (std::size_t reaction_index = 0; reaction_index < network_->size();
         ++reaction_index) {
      const auto id = ReactionId{static_cast<std::uint32_t>(reaction_index)};
      validateMassActionReactionRate(
          resolved_terms_[reaction_index],
          RateCoefficientView{rate_coefficients[reaction_index]},
          reaction_rates[id]);
    }
    for (std::size_t reaction_index = 0; reaction_index < network_->size();
         ++reaction_index) {
      detail::writeValidatedMassActionReactionRate(
          resolved_terms_[reaction_index],
          RateCoefficientView{rate_coefficients[reaction_index]},
          reaction_rates[ReactionId{
              static_cast<std::uint32_t>(reaction_index)}]);
    }
  }

 private:
  /** @brief Returns whether one field declares a runtime physical quantity. */
  [[nodiscard]] static bool hasQuantity(
      const field::FieldMetadata& metadata) noexcept {
    return metadata.physical_quantity.has_value();
  }

  /** @brief Requires one field to declare the requested semantic kind. */
  static void requireKind(const field::FieldMetadata& metadata,
                          unit::QuantityKind kind, std::string_view name) {
    if (!metadata.physical_quantity.has_value() ||
        metadata.physical_quantity->kind() != kind) {
      throw std::invalid_argument(std::string{name} +
                                  " metadata quantity kind differs");
    }
  }

  /** @brief Converts a validated integral kinetic order to an integer. */
  [[nodiscard]] static int integralOrder(double order) {
    const double rounded = std::round(order);
    if (std::abs(order - rounded) > 1e-12 || rounded < 1.0 ||
        rounded > static_cast<double>(std::numeric_limits<int>::max())) {
      throw std::invalid_argument(
          "physical metadata requires positive integral kinetic orders");
    }
    return static_cast<int>(rounded);
  }

  /** @brief Validates coefficient dimensions implied by every kinetic order. */
  void validateMetadata(
      std::span<const field::CellField<double>> rate_coefficients,
      const ReactionRateFields& reaction_rates) const {
    bool any_quantity = false;
    for (const auto& density : *density_) {
      any_quantity = any_quantity || hasQuantity(density.metadata());
    }
    for (std::size_t reaction_index = 0; reaction_index < network_->size();
         ++reaction_index) {
      any_quantity =
          any_quantity ||
          hasQuantity(rate_coefficients[reaction_index].metadata()) ||
          hasQuantity(reaction_rates[ReactionId{static_cast<std::uint32_t>(
                                         reaction_index)}]
                          .metadata());
    }

    if (!any_quantity) {
      if (!allow_unitless_raw_values_) {
        throw std::invalid_argument(
            "mass-action assembly requires physical metadata or an explicit "
            "unitless raw-value opt-in");
      }
      return;
    }

    for (std::size_t reaction_index = 0; reaction_index < network_->size();
         ++reaction_index) {
      const auto id = ReactionId{static_cast<std::uint32_t>(reaction_index)};
      const auto& rate = reaction_rates[id];
      const auto& coefficient = rate_coefficients[reaction_index];
      requireKind(rate.metadata(), unit::QuantityKind::reaction_rate_density,
                  "reaction rate");
      requireKind(coefficient.metadata(),
                  unit::QuantityKind::reaction_rate_coefficient,
                  "reaction rate coefficient");

      auto expected_coefficient_unit =
          rate.metadata().physical_quantity->unit();
      for (const auto& kinetic : network_->at(id).kinetic_orders) {
        const auto& species_density = (*density_)[kinetic.species];
        requireKind(species_density.metadata(),
                    unit::QuantityKind::particle_number_density,
                    "reactant density");
        expected_coefficient_unit =
            expected_coefficient_unit /
            species_density.metadata().physical_quantity->unit().pow(
                integralOrder(kinetic.order));
      }
      if (coefficient.metadata().physical_quantity->unit() !=
          expected_coefficient_unit) {
        throw std::invalid_argument(
            "reaction rate coefficient unit differs from kinetic order");
      }
    }
  }

  const ReactionNetwork* network_;
  const physics::SpeciesCellFields* density_;
  bool allow_unitless_raw_values_;
  std::vector<std::vector<MassActionTerm>> resolved_terms_;
};

}  // namespace pemu::physics::reaction
