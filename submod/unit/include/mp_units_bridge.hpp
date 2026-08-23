#pragma once

#include <pemu/unit/plasma_quantities.hpp>
#include <pemu/unit/quantity_metadata.hpp>

#include <mp-units/framework/reference.h>
#include <mp-units/framework/unit.h>
#include <mp-units/systems/si.h>

#include <type_traits>

namespace pemu::unit {
namespace detail {

template <typename>
inline constexpr bool dependent_false = false;

template <typename QuantitySpec>
[[nodiscard]] consteval QuantityKind quantityKindOf(QuantitySpec) {
  using QS = std::remove_cvref_t<QuantitySpec>;
  using namespace mp_units;

  if constexpr (std::same_as<QS,
                             std::remove_cvref_t<decltype(dimensionless)>>) {
    return QuantityKind::dimensionless;
  } else if constexpr (std::same_as<
                           QS, std::remove_cvref_t<decltype(isq::length)>>) {
    return QuantityKind::length;
  } else if constexpr (std::same_as<
                           QS, std::remove_cvref_t<decltype(isq::repetency)>>) {
    return QuantityKind::inverse_length;
  } else if constexpr (std::same_as<
                           QS, std::remove_cvref_t<decltype(isq::frequency)>>) {
    return QuantityKind::frequency;
  } else if constexpr (std::same_as<
                           QS, std::remove_cvref_t<decltype(isq::speed)>>) {
    return QuantityKind::speed;
  } else if constexpr (
      std::same_as<QS,
                   std::remove_cvref_t<
                       decltype(plasma_quantity::particle_number_density)>>) {
    return QuantityKind::particle_number_density;
  } else if constexpr (
      std::same_as<
          QS, std::remove_cvref_t<
                  decltype(plasma_quantity::particle_number_density_rate)>>) {
    return QuantityKind::particle_number_density_rate;
  } else if constexpr (
      std::same_as<QS, std::remove_cvref_t<
                           decltype(plasma_quantity::reaction_rate_density)>>) {
    return QuantityKind::reaction_rate_density;
  } else if constexpr (std::same_as<
                           QS, std::remove_cvref_t<
                                   decltype(isq::electric_charge_density)>>) {
    return QuantityKind::electric_charge_density;
  } else if constexpr (std::same_as<QS,
                                    std::remove_cvref_t<
                                        decltype(isq::electric_potential)>>) {
    return QuantityKind::electric_potential;
  } else if constexpr (std::same_as<
                           QS,
                           std::remove_cvref_t<
                               decltype(isq::electric_potential_difference)>>) {
    return QuantityKind::electric_potential_difference;
  } else if constexpr (
      std::same_as<
          QS, std::remove_cvref_t<
                  decltype(plasma_quantity::normal_electric_field_strength)>>) {
    return QuantityKind::normal_electric_field_strength;
  } else if constexpr (
      std::same_as<QS, std::remove_cvref_t<
                           decltype(plasma_quantity::normal_drift_velocity)>>) {
    return QuantityKind::normal_drift_velocity;
  } else if constexpr (
      std::same_as<QS, std::remove_cvref_t<
                           decltype(plasma_quantity::electron_mean_energy)>>) {
    return QuantityKind::electron_mean_energy;
  } else if constexpr (
      std::same_as<QS,
                   std::remove_cvref_t<
                       decltype(plasma_quantity::electron_energy_density)>>) {
    return QuantityKind::electron_energy_density;
  } else if constexpr (
      std::same_as<
          QS, std::remove_cvref_t<
                  decltype(plasma_quantity::electron_energy_density_rate)>>) {
    return QuantityKind::electron_energy_density_rate;
  } else {
    static_assert(dependent_false<QS>,
                  "mp-units quantity specification is not bridged to "
                  "pemu::unit::QuantityKind");
  }
}

template <typename MpUnit>
[[nodiscard]] consteval units::precise_unit preciseUnitOf(MpUnit) {
  using U = std::remove_cvref_t<MpUnit>;
  using namespace mp_units;
  using namespace mp_units::si::unit_symbols;

  if constexpr (std::same_as<U, std::remove_cvref_t<decltype(one)>>) {
    return units::precise::one;
  } else if constexpr (std::same_as<U, std::remove_cvref_t<decltype(m)>>) {
    return units::precise::m;
  } else if constexpr (std::same_as<U, std::remove_cvref_t<decltype(cm)>>) {
    return units::precise::cm;
  } else if constexpr (std::same_as<U, std::remove_cvref_t<decltype(mm)>>) {
    return units::precise::mm;
  } else if constexpr (std::same_as<U, std::remove_cvref_t<decltype(s)>>) {
    return units::precise::s;
  } else if constexpr (std::same_as<U, std::remove_cvref_t<decltype(C)>>) {
    return units::precise::C;
  } else if constexpr (std::same_as<U, std::remove_cvref_t<decltype(V)>>) {
    return units::precise::V;
  } else if constexpr (std::same_as<U,
                                    std::remove_cvref_t<decltype(one / m)>>) {
    return units::precise::one / units::precise::m;
  } else if constexpr (std::same_as<U,
                                    std::remove_cvref_t<decltype(one / s)>>) {
    return units::precise::one / units::precise::s;
  } else if constexpr (std::same_as<U, std::remove_cvref_t<decltype(m / s)>>) {
    return units::precise::m / units::precise::s;
  } else if constexpr (std::same_as<U, std::remove_cvref_t<decltype(cm / s)>>) {
    return units::precise::cm / units::precise::s;
  } else if constexpr (std::same_as<
                           U, std::remove_cvref_t<decltype(one / cubic(cm))>>) {
    return units::precise::one / units::precise::cm.pow(3);
  } else if constexpr (std::same_as<U, std::remove_cvref_t<
                                           decltype(one / (cubic(cm) * s))>>) {
    return units::precise::one /
           (units::precise::cm.pow(3) * units::precise::s);
  } else if constexpr (std::same_as<
                           U, std::remove_cvref_t<decltype(C / cubic(cm))>>) {
    return units::precise::C / units::precise::cm.pow(3);
  } else if constexpr (std::same_as<U, std::remove_cvref_t<decltype(V / cm)>>) {
    return units::precise::V / units::precise::cm;
  } else if constexpr (std::same_as<U, std::remove_cvref_t<decltype(eV)>>) {
    return units::precise::energy::eV;
  } else if constexpr (std::same_as<
                           U, std::remove_cvref_t<decltype(eV / cubic(cm))>>) {
    return units::precise::energy::eV / units::precise::cm.pow(3);
  } else if constexpr (std::same_as<U, std::remove_cvref_t<
                                           decltype(eV / (cubic(cm) * s))>>) {
    return units::precise::energy::eV /
           (units::precise::cm.pow(3) * units::precise::s);
  } else {
    static_assert(dependent_false<U>,
                  "mp-units unit is not bridged to LLNL units::precise_unit");
  }
}

}  // namespace detail

template <mp_units::Reference R>
[[nodiscard]] consteval PhysicalQuantityMetadata bridgeReference(R) {
  return {detail::quantityKindOf(mp_units::get_quantity_spec(R{})),
          detail::preciseUnitOf(mp_units::get_unit(R{}))};
}

template <mp_units::Reference R>
[[nodiscard]] consteval units::precise_unit bridgeUnit(R) {
  return detail::preciseUnitOf(mp_units::get_unit(R{}));
}

}  // namespace pemu::unit
