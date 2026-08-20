#pragma once

#include <pemu/field/mp_units_metadata.hpp>
#include <pemu/unit/plasma_quantities.hpp>

#include <mp-units/systems/si.h>

namespace pemu::field {

struct PlasmaFieldMetadata {
  FieldMetadata number_density;
  FieldMetadata number_density_source;
  FieldMetadata reaction_rate;
  FieldMetadata charge_density;
  FieldMetadata electric_potential;
  FieldMetadata electric_field;
  FieldMetadata drift_velocity;
  FieldMetadata inverse_time;
};

// Canonical metadata for meshes whose coordinate values are centimetres. The
// associated Field<double> objects still store plain numerical values expressed
// in these units.
[[nodiscard]] inline PlasmaFieldMetadata centimetrePlasmaFieldMetadata() {
  using namespace mp_units;
  using namespace mp_units::si::unit_symbols;

  constexpr auto number_density_unit = one / cubic(cm);
  constexpr auto rate_density_unit = one / (cubic(cm) * s);

  return {
      .number_density = makeFieldMetadata(
          "number density", pemu::unit::plasma_quantity::particle_number_density
                                [number_density_unit]),
      .number_density_source = makeFieldMetadata(
          "number density source",
          pemu::unit::plasma_quantity::particle_number_density_rate
              [rate_density_unit]),
      .reaction_rate = makeFieldMetadata(
          "reaction rate", pemu::unit::plasma_quantity::reaction_rate_density
                               [rate_density_unit]),
      .charge_density = makeFieldMetadata(
          "charge density", isq::electric_charge_density[C / cubic(cm)]),
      .electric_potential =
          makeFieldMetadata("electric potential", isq::electric_potential[V]),
      .electric_field = makeFieldMetadata(
          "normal electric field",
          pemu::unit::plasma_quantity::normal_electric_field_strength[V / cm]),
      .drift_velocity = makeFieldMetadata(
          "normal drift velocity",
          pemu::unit::plasma_quantity::normal_drift_velocity[cm / s]),
      .inverse_time =
          makeFieldMetadata("inverse time", isq::frequency[one / s]),
  };
}

}  // namespace pemu::field
