#pragma once

#include <pemu/field/cell_field.hpp>
#include <pemu/field/utils.hpp>

#include <stdexcept>

namespace pemu::physics::reaction {

/**
* @brief 
* 
* @param electron_density 
* @param neutral_density 
* @param rate_coefficient 
* @param reaction_rate output 
*/
inline void electronImpactIonizationRate(
    const field::CellField<double>& electron_density, double neutral_density,
    double rate_coefficient, field::CellField<double>& reaction_rate) {
  field::ensureSameMesh(electron_density, reaction_rate);
  if (neutral_density < 0.0) {
    throw std::invalid_argument("neutral density must be non-negative");
  }
  if (rate_coefficient < 0.0) {
    throw std::invalid_argument(
        "ionization rate coefficient must be non-negative");
  }
  const auto& mesh = electron_density.mesh();
  for (mesh::CellId cell = 0; cell < mesh.numCells(); ++cell) {
    reaction_rate[cell] =
        rate_coefficient * electron_density[cell] * neutral_density;
  }
}

inline void addPairProductionSource(
    const field::CellField<double>& reaction_rate,
    field::CellField<double>& electron_source,
    field::CellField<double>& ion_source) {
  field::ensureSameMesh(reaction_rate, electron_source);
  field::ensureSameMesh(reaction_rate, ion_source);
  const auto& mesh = reaction_rate.mesh();
  for (mesh::CellId cell = 0; cell < mesh.numCells(); ++cell) {
    const double rate = reaction_rate[cell];
    electron_source[cell] += rate;
    ion_source[cell] += rate;
  }
}

}  // namespace pemu::physics::reaction