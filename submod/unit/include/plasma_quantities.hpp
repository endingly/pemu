#pragma once

#include <mp-units/compat_macros.h>
#include <mp-units/framework/quantity_spec.h>
#include <mp-units/systems/isq/electromagnetism.h>
#include <mp-units/systems/isq/space_and_time.h>

namespace pemu::unit::plasma_quantity {

QUANTITY_SPEC(particle_number_density,
              mp_units::inverse(mp_units::isq::volume));
QUANTITY_SPEC(particle_number_density_rate,
              particle_number_density / mp_units::isq::time);
QUANTITY_SPEC(reaction_rate_density,
              particle_number_density / mp_units::isq::time);
QUANTITY_SPEC(normal_electric_field_strength,
              mp_units::isq::electric_potential / mp_units::isq::length,
              mp_units::quantity_character::real_scalar);
QUANTITY_SPEC(normal_drift_velocity,
              mp_units::isq::length / mp_units::isq::time,
              mp_units::quantity_character::real_scalar);
QUANTITY_SPEC(electron_mean_energy, mp_units::isq::energy);
QUANTITY_SPEC(electron_energy_density,
              electron_mean_energy / mp_units::isq::volume);
QUANTITY_SPEC(electron_energy_density_rate,
              electron_energy_density / mp_units::isq::time);

}  // namespace pemu::unit::plasma_quantity
