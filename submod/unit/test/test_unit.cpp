#include <gtest/gtest.h>

#include <pemu/unit/mp_units_bridge.hpp>
#include <pemu/unit/plasma_quantities.hpp>
#include <pemu/unit/quantity_metadata.hpp>

#include <mp-units/framework/quantity.h>
#include <mp-units/systems/si.h>

#include <type_traits>

namespace pemu::unit::test {

template <typename QuantitySpec, typename Unit>
concept FormsQuantityReference =
    requires(std::remove_cvref_t<QuantitySpec> quantity_spec,
             std::remove_cvref_t<Unit> unit) { quantity_spec[unit]; };

TEST(PlasmaQuantitiesTest, DefineDimensionallyConsistentReferences) {
  using namespace mp_units;
  using namespace mp_units::si::unit_symbols;

  static_assert(
      FormsQuantityReference<decltype(plasma_quantity::particle_number_density),
                             decltype(one / cubic(cm))>);
  static_assert(FormsQuantityReference<
                decltype(plasma_quantity::particle_number_density_rate),
                decltype(one / (cubic(cm) * s))>);
  static_assert(
      FormsQuantityReference<decltype(plasma_quantity::reaction_rate_density),
                             decltype(one / (cubic(cm) * s))>);
  static_assert(FormsQuantityReference<
                decltype(plasma_quantity::normal_electric_field_strength),
                decltype(V / cm)>);
  static_assert(
      FormsQuantityReference<decltype(plasma_quantity::normal_drift_velocity),
                             decltype(cm / s)>);
  static_assert(
      FormsQuantityReference<decltype(plasma_quantity::electron_mean_energy),
                             decltype(eV)>);
  static_assert(
      FormsQuantityReference<decltype(plasma_quantity::electron_energy_density),
                             decltype(eV / cubic(cm))>);
  static_assert(FormsQuantityReference<
                decltype(plasma_quantity::electron_energy_density_rate),
                decltype(eV / (cubic(cm) * s))>);

  static_assert(
      !FormsQuantityReference<
          decltype(plasma_quantity::particle_number_density), decltype(V)>);
  static_assert(
      !FormsQuantityReference<decltype(plasma_quantity::normal_drift_velocity),
                              decltype(V / cm)>);
}

TEST(MpUnitsBridgeTest, ProducesRuntimeMetadataWithoutStringParsing) {
  using namespace mp_units;
  using namespace mp_units::si::unit_symbols;

  constexpr auto density = bridgeReference(
      plasma_quantity::particle_number_density[one / cubic(cm)]);
  constexpr auto charge =
      bridgeReference(isq::electric_charge_density[C / cubic(cm)]);
  constexpr auto potential = bridgeReference(isq::electric_potential[V]);
  constexpr auto electric_field =
      bridgeReference(plasma_quantity::normal_electric_field_strength[V / cm]);
  constexpr auto energy_density =
      bridgeReference(plasma_quantity::electron_energy_density[eV / cubic(cm)]);

  static_assert(density.kind() == QuantityKind::particle_number_density);
  static_assert(charge.kind() == QuantityKind::electric_charge_density);
  static_assert(potential.kind() == QuantityKind::electric_potential);
  static_assert(electric_field.kind() ==
                QuantityKind::normal_electric_field_strength);
  static_assert(energy_density.kind() == QuantityKind::electron_energy_density);
  static_assert(std::is_trivially_copyable_v<PhysicalQuantityMetadata>);

  EXPECT_EQ(density.unit(), units::precise::one / units::precise::cm.pow(3));
  EXPECT_EQ(charge.unit(), units::precise::C / units::precise::cm.pow(3));
  EXPECT_EQ(potential.unit(), units::precise::V);
  EXPECT_EQ(electric_field.unit(), units::precise::V / units::precise::cm);
  EXPECT_EQ(energy_density.unit(),
            units::precise::energy::eV / units::precise::cm.pow(3));
  EXPECT_EQ(pemu::to_string(electric_field.kind()),
            "normal_electric_field_strength");
}

TEST(RuntimeUnitStringTest, UsesStablePlasmaEnergySpellings) {
  const auto energy_density =
      units::precise::energy::eV / units::precise::cm.pow(3);
  const auto energy_source =
      energy_density / units::precise::s;

  EXPECT_EQ(pemu::to_string(energy_density), "eV/cm^3");
  EXPECT_EQ(pemu::to_string(energy_source), "eV/(cm^3*s)");
  EXPECT_EQ(pemu::to_string(units::precise::energy::eV / units::precise::s),
            "eV/s");
  EXPECT_EQ(pemu::to_string(units::precise::V / units::precise::cm),
            units::to_string(units::precise::V / units::precise::cm));
}

}  // namespace pemu::unit::test
