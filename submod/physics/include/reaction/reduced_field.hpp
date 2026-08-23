#pragma once

#include <pemu/field/cell_field.hpp>
#include <pemu/field/utils.hpp>

#include <cmath>
#include <stdexcept>

namespace pemu::physics::reaction {

/** Number of V cm^2 represented by one Townsend. */
inline constexpr double volt_square_centimetres_per_townsend = 1.0e-17;

/** @brief Computes E/N in Townsend from |E| in V/cm and N in 1/cm^3. */
[[nodiscard]] inline double reducedElectricFieldTownsend(
    double electric_field_magnitude_v_per_cm,
    double neutral_number_density_per_cm3) {
  if (!std::isfinite(electric_field_magnitude_v_per_cm) ||
      electric_field_magnitude_v_per_cm < 0.0) {
    throw std::invalid_argument(
        "electric-field magnitude must be finite and non-negative");
  }
  if (!std::isfinite(neutral_number_density_per_cm3) ||
      neutral_number_density_per_cm3 <= 0.0) {
    throw std::invalid_argument(
        "neutral number density must be finite and positive");
  }

  const double value = electric_field_magnitude_v_per_cm /
                       neutral_number_density_per_cm3 /
                       volt_square_centimetres_per_townsend;
  if (!std::isfinite(value)) {
    throw std::overflow_error("reduced electric field overflowed");
  }
  return value;
}

/** @brief Computes a cell E/N field using canonical centimetre-based values. */
inline void computeReducedElectricFieldTownsend(
    const field::CellField<double>& electric_field_magnitude_v_per_cm,
    const field::CellField<double>& neutral_number_density_per_cm3,
    field::CellField<double>& reduced_electric_field_td) {
  field::ensureSameMesh(electric_field_magnitude_v_per_cm,
                        neutral_number_density_per_cm3);
  field::ensureSameMesh(electric_field_magnitude_v_per_cm,
                        reduced_electric_field_td);

  for (mesh::CellId cell = 0;
       cell < electric_field_magnitude_v_per_cm.mesh().numCells(); ++cell) {
    static_cast<void>(
        reducedElectricFieldTownsend(electric_field_magnitude_v_per_cm[cell],
                                     neutral_number_density_per_cm3[cell]));
  }
  for (mesh::CellId cell = 0;
       cell < electric_field_magnitude_v_per_cm.mesh().numCells(); ++cell) {
    reduced_electric_field_td[cell] =
        reducedElectricFieldTownsend(electric_field_magnitude_v_per_cm[cell],
                                     neutral_number_density_per_cm3[cell]);
  }
}

}  // namespace pemu::physics::reaction
