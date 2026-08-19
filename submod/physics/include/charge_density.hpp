#pragma once
#include <pemu/field/cell_field.hpp>
#include <pemu/field/utils.hpp>

namespace pemu::physics {

inline void addSpeciesChargeDensity(
    const field::CellField<double>& number_density, double particle_charge,
    field::CellField<double>& charge_density) {
  field::ensureSameMesh(number_density, charge_density);

  for (mesh::CellId cell = 0; cell < number_density.mesh().numCells(); ++cell) {

    charge_density[cell] += particle_charge * number_density[cell];
  }
}

inline void clearChargeDensity(field::CellField<double>& rho) {
  rho.fill(0.0);
}

};  // namespace pemu::physics