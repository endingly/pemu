#pragma once

#include <pemu/field/cell_field.hpp>
#include <pemu/physics/reaction/rate_coefficient.hpp>
#include <pemu/physics/reaction/types.hpp>
#include <pemu/physics/species.hpp>
#include <pemu/simulation/plasma_reaction_rate_evaluator.hpp>

namespace pemu::simulation {

/** @brief Selects the physical coordinate used by an electron-impact table. */
enum class ElectronImpactRateCoordinate {
  reduced_electric_field_townsend,
  electron_temperature_ev,
};

/**
 * @brief Evaluates one tabulated binary electron-impact reaction per cell.
 *
 * The target-density field is externally owned and may represent an immobile
 * neutral species or prescribed gas state. Its numerical unit is 1/cm^3.
 * Tables selected by reduced field use a Td coordinate; temperature tables use
 * Te in eV. Coefficients are expected in cm^3/s, producing rates in
 * 1/(cm^3*s).
 */
class TabulatedElectronImpactEvaluator {
 public:
  /** @brief Binds species, reaction, gas density, and immutable lookup data. */
  TabulatedElectronImpactEvaluator(
      physics::SpeciesId electron, physics::reaction::ReactionId reaction,
      const field::CellField<double>& target_number_density,
      ElectronImpactRateCoordinate coordinate,
      physics::reaction::TabulatedRateCoefficient table,
      bool allow_unitless_raw_values = false);

  /** @brief Overwrites the selected reaction field from the current context. */
  void operator()(const PlasmaReactionRateContext& context,
                  physics::reaction::ReactionRateFields& reaction_rates);

 private:
  /** @brief Enforces the canonical centimetre/eV table unit convention. */
  void validateMetadata(const PlasmaReactionRateContext& context,
                        const field::CellField<double>& reaction_rate) const;

  physics::SpeciesId electron_;
  physics::reaction::ReactionId reaction_;
  const field::CellField<double>* target_number_density_;
  ElectronImpactRateCoordinate coordinate_;
  physics::reaction::TabulatedRateCoefficient table_;
  bool allow_unitless_raw_values_;
  field::CellField<double> electric_field_magnitude_;
  field::CellField<double> table_coordinate_;
  field::CellField<double> rate_coefficient_;
};

}  // namespace pemu::simulation
