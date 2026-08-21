#pragma once

#include <llnl-units/units.hpp>

#include <cstdint>
#include <string_view>

namespace pemu::unit {

enum class QuantityKind : std::uint8_t {
  dimensionless,
  length,
  inverse_length,
  frequency,
  speed,
  particle_number_density,
  particle_number_density_rate,
  reaction_rate_density,
  electric_charge_density,
  electric_potential,
  electric_potential_difference,
  normal_electric_field_strength,
  normal_drift_velocity,
};

// Runtime metadata produced once from an mp-units reference at a Field
// construction boundary. No mp-units type and no presentation string is kept.
class PhysicalQuantityMetadata {
 public:
  constexpr PhysicalQuantityMetadata(QuantityKind kind,
                                     units::precise_unit unit) noexcept
      : kind_(kind), unit_(unit) {}

  [[nodiscard]] constexpr QuantityKind kind() const noexcept { return kind_; }

  [[nodiscard]] constexpr units::precise_unit unit() const noexcept {
    return unit_;
  }

  friend bool operator==(const PhysicalQuantityMetadata& lhs,
                         const PhysicalQuantityMetadata& rhs) {
    return lhs.kind_ == rhs.kind_ && lhs.unit_ == rhs.unit_;
  }

 private:
  QuantityKind kind_;
  units::precise_unit unit_;
};

}  // namespace pemu::unit

namespace pemu {

[[nodiscard]] constexpr std::string_view to_string(
    unit::QuantityKind kind) noexcept {
  using enum unit::QuantityKind;
  switch (kind) {
    case dimensionless:
      return "dimensionless";
    case length:
      return "length";
    case inverse_length:
      return "inverse_length";
    case frequency:
      return "frequency";
    case speed:
      return "speed";
    case particle_number_density:
      return "particle_number_density";
    case particle_number_density_rate:
      return "particle_number_density_rate";
    case reaction_rate_density:
      return "reaction_rate_density";
    case electric_charge_density:
      return "electric_charge_density";
    case electric_potential:
      return "electric_potential";
    case electric_potential_difference:
      return "electric_potential_difference";
    case normal_electric_field_strength:
      return "normal_electric_field_strength";
    case normal_drift_velocity:
      return "normal_drift_velocity";
  }
  return "unknown";
}

}  // namespace pemu
