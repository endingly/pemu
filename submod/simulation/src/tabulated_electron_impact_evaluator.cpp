#include <pemu/simulation/tabulated_electron_impact_evaluator.hpp>

#include <pemu/discretization/operators/cell_vector_reconstruction.hpp>
#include <pemu/field/utils.hpp>
#include <pemu/physics/electron_energy.hpp>
#include <pemu/physics/reaction/mass_action.hpp>
#include <pemu/physics/reaction/reduced_field.hpp>
#include <pemu/unit/quantity_metadata.hpp>

#include <llnl-units/units.hpp>

#include <array>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace pemu::simulation {

/** @copydoc TabulatedElectronImpactEvaluator::TabulatedElectronImpactEvaluator */
TabulatedElectronImpactEvaluator::TabulatedElectronImpactEvaluator(
    physics::SpeciesId electron, physics::reaction::ReactionId reaction,
    const field::CellField<double>& target_number_density,
    ElectronImpactRateCoordinate coordinate,
    physics::reaction::TabulatedRateCoefficient table,
    bool allow_unitless_raw_values)
    : electron_(electron),
      reaction_(reaction),
      target_number_density_(&target_number_density),
      coordinate_(coordinate),
      table_(std::move(table)),
      allow_unitless_raw_values_(allow_unitless_raw_values),
      electric_field_magnitude_(target_number_density.mesh(), 0.0),
      table_coordinate_(target_number_density.mesh(), 0.0),
      rate_coefficient_(
          target_number_density.mesh(), 0.0,
          {.name = "binary electron-impact rate coefficient",
           .physical_quantity = unit::PhysicalQuantityMetadata{
               unit::QuantityKind::reaction_rate_coefficient,
               physics::reaction::canonicalCentimetreRateCoefficientUnit(2)}}) {
  if (coordinate_ !=
          ElectronImpactRateCoordinate::reduced_electric_field_townsend &&
      coordinate_ != ElectronImpactRateCoordinate::electron_temperature_ev) {
    throw std::invalid_argument(
        "electron-impact rate coordinate is not supported");
  }
}

/** @copydoc TabulatedElectronImpactEvaluator::validateMetadata */
void TabulatedElectronImpactEvaluator::validateMetadata(
    const PlasmaReactionRateContext& context,
    const field::CellField<double>& reaction_rate) const {
  const auto& electron_density = context.density[electron_];
  const auto has_quantity = [](const field::FieldMetadata& metadata) {
    return metadata.physical_quantity.has_value();
  };
  const bool any_quantity =
      has_quantity(electron_density.metadata()) ||
      has_quantity(target_number_density_->metadata()) ||
      has_quantity(reaction_rate.metadata()) ||
      (coordinate_ ==
               ElectronImpactRateCoordinate::reduced_electric_field_townsend
           ? has_quantity(context.electric_field_normal.metadata())
           : has_quantity(context.electron_mean_energy.metadata()));

  if (!any_quantity) {
    if (!allow_unitless_raw_values_) {
      throw std::invalid_argument(
          "tabulated electron-impact chemistry requires physical metadata or "
          "an explicit unitless raw-value opt-in");
    }
    return;
  }

  const auto require =
      [](const field::FieldMetadata& metadata, unit::QuantityKind kind,
         units::precise_unit expected_unit, std::string_view name) {
        if (!metadata.physical_quantity.has_value() ||
            metadata.physical_quantity->kind() != kind ||
            metadata.physical_quantity->unit() != expected_unit) {
          throw std::invalid_argument(std::string{name} +
                                      " metadata is not in the canonical M14 "
                                      "unit");
        }
      };
  const auto density_unit = units::precise::one / units::precise::cm.pow(3);
  require(electron_density.metadata(),
          unit::QuantityKind::particle_number_density, density_unit,
          "electron density");
  require(target_number_density_->metadata(),
          unit::QuantityKind::particle_number_density, density_unit,
          "target density");
  require(reaction_rate.metadata(), unit::QuantityKind::reaction_rate_density,
          density_unit / units::precise::s, "reaction rate");

  if (coordinate_ ==
      ElectronImpactRateCoordinate::reduced_electric_field_townsend) {
    require(context.electric_field_normal.metadata(),
            unit::QuantityKind::normal_electric_field_strength,
            units::precise::V / units::precise::cm, "electric field");
  } else {
    require(context.electron_mean_energy.metadata(),
            unit::QuantityKind::electron_mean_energy,
            units::precise::energy::eV, "electron mean energy");
  }
}

/** @copydoc TabulatedElectronImpactEvaluator::operator() */
void TabulatedElectronImpactEvaluator::operator()(
    const PlasmaReactionRateContext& context,
    physics::reaction::ReactionRateFields& reaction_rates) {
  field::ensureSameMesh(context.density[electron_], *target_number_density_);
  field::ensureSameMesh(context.density[electron_], context.potential);
  field::ensureSameMesh(context.density[electron_],
                        context.electric_field_normal);
  field::ensureSameMesh(context.density[electron_],
                        context.electron_mean_energy);
  field::ensureSameMesh(context.density[electron_], reaction_rates[reaction_]);
  validateMetadata(context, reaction_rates[reaction_]);

  switch (coordinate_) {
    case ElectronImpactRateCoordinate::reduced_electric_field_townsend:
      discretization::operators::reconstructCellVectorMagnitudeFromFaceNormal(
          context.electric_field_normal, electric_field_magnitude_);
      physics::reaction::computeReducedElectricFieldTownsend(
          electric_field_magnitude_, *target_number_density_,
          table_coordinate_);
      break;
    case ElectronImpactRateCoordinate::electron_temperature_ev:
      physics::computeElectronTemperatureEv(
          context.electron_mean_energy, table_coordinate_);
      break;
  }

  physics::reaction::evaluateRateCoefficient(table_coordinate_, table_,
                                             rate_coefficient_);
  const std::array reactants{
      physics::reaction::MassActionTerm{context.density[electron_]},
      physics::reaction::MassActionTerm{*target_number_density_}};
  physics::reaction::massActionReactionRate(
      reactants, physics::reaction::RateCoefficientView{rate_coefficient_},
      reaction_rates[reaction_]);
}

}  // namespace pemu::simulation
