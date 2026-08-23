#include <gtest/gtest.h>

#include <pemu/field/cell_field.hpp>
#include <pemu/field/plasma_field_metadata.hpp>
#include <pemu/mesh/moab_mesh.hpp>
#include <pemu/physics/reaction/mass_action.hpp>
#include <pemu/physics/reaction/mass_action_assembler.hpp>
#include <pemu/physics/reaction/network.hpp>
#include <pemu/physics/reaction/rate_coefficient.hpp>
#include <pemu/physics/reaction/reduced_field.hpp>
#include <pemu/physics/species.hpp>

#include <array>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <vector>

namespace pemu::physics::reaction::test {
namespace {

/** @brief Returns the two-cell mesh shared by reaction-kinetics tests. */
[[nodiscard]] std::filesystem::path testMeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
}

class ReactionKineticsTest : public ::testing::Test {
 protected:
  /** @brief Loads two cells for spatial and atomic-failure checks. */
  ReactionKineticsTest() : mesh_(testMeshPath().string()) {}

  mesh::MoabMesh mesh_;
};

}  // namespace

TEST(ReducedElectricFieldTest, UsesTownsendCentimetreConversion) {
  EXPECT_NEAR(reducedElectricFieldTownsend(400.0, 2.5e19), 1.6, 1e-14);
  EXPECT_DOUBLE_EQ(reducedElectricFieldTownsend(0.0, 2.5e19), 0.0);
  EXPECT_THROW(static_cast<void>(reducedElectricFieldTownsend(1.0, 0.0)),
               std::invalid_argument);
  EXPECT_THROW(static_cast<void>(reducedElectricFieldTownsend(-1.0, 1.0)),
               std::invalid_argument);
}

TEST(TabulatedRateCoefficientTest, InterpolatesLinearAndLogLinearValues) {
  const TabulatedRateCoefficient linear({0.0, 10.0}, {2.0, 6.0},
                                        RateInterpolation::linear);
  const TabulatedRateCoefficient logarithmic({0.0, 10.0}, {1.0e-12, 1.0e-10});

  EXPECT_DOUBLE_EQ(linear(0.0), 2.0);
  EXPECT_DOUBLE_EQ(linear(5.0), 4.0);
  EXPECT_DOUBLE_EQ(linear(10.0), 6.0);
  EXPECT_NEAR(logarithmic(5.0), 1.0e-11, 1.0e-25);
}

TEST(TabulatedRateCoefficientTest, AppliesExplicitBoundsPolicy) {
  const TabulatedRateCoefficient strict({1.0, 2.0}, {3.0, 4.0},
                                        RateInterpolation::linear);
  const TabulatedRateCoefficient clamped({1.0, 2.0}, {3.0, 4.0},
                                         RateInterpolation::linear,
                                         RateTableBounds::clamp);

  EXPECT_THROW(static_cast<void>(strict(0.5)), std::out_of_range);
  EXPECT_THROW(static_cast<void>(strict(2.5)), std::out_of_range);
  EXPECT_DOUBLE_EQ(clamped(0.5), 3.0);
  EXPECT_DOUBLE_EQ(clamped(2.5), 4.0);
}

TEST(TabulatedRateCoefficientTest, RejectsMalformedTables) {
  EXPECT_THROW(TabulatedRateCoefficient({1.0}, {2.0}), std::invalid_argument);
  EXPECT_THROW(TabulatedRateCoefficient({1.0, 1.0}, {2.0, 3.0}),
               std::invalid_argument);
  EXPECT_THROW(TabulatedRateCoefficient({1.0, 2.0}, {2.0, -3.0}),
               std::invalid_argument);
  EXPECT_THROW(TabulatedRateCoefficient({1.0, 2.0}, {0.0, 3.0}),
               std::invalid_argument);
  EXPECT_NO_THROW(TabulatedRateCoefficient({1.0, 2.0}, {0.0, 3.0},
                                           RateInterpolation::linear));
}

TEST_F(ReactionKineticsTest, EvaluatesArbitraryThreeBodyMassActionRate) {
  field::CellField<double> electron_density(mesh_, 2.0);
  field::CellField<double> oxygen_density(mesh_, 3.0);
  field::CellField<double> third_body_density(mesh_, 5.0);
  field::CellField<double> rate_coefficient(mesh_, 7.0);
  field::CellField<double> reaction_rate(mesh_, 0.0);
  const std::array reactants{MassActionTerm{electron_density},
                             MassActionTerm{oxygen_density},
                             MassActionTerm{third_body_density}};

  massActionReactionRate(reactants, RateCoefficientView{rate_coefficient},
                         reaction_rate);

  EXPECT_DOUBLE_EQ(totalReactionOrder(reactants), 3.0);
  for (const double value : reaction_rate) {
    EXPECT_DOUBLE_EQ(value, 210.0);
  }
}

TEST_F(ReactionKineticsTest, SupportsFractionalEmpiricalReactionOrder) {
  field::CellField<double> density(mesh_, 9.0);
  field::CellField<double> reaction_rate(mesh_, 0.0);
  const std::array reactants{MassActionTerm{density, 0.5}};

  massActionReactionRate(reactants, RateCoefficientView{2.0}, reaction_rate);

  for (const double value : reaction_rate) {
    EXPECT_DOUBLE_EQ(value, 6.0);
  }
}

TEST_F(ReactionKineticsTest,
       ZeroFactorShortCircuitsAnOtherwiseOverflowingRate) {
  field::CellField<double> enormous_density(mesh_,
                                            std::numeric_limits<double>::max());
  field::CellField<double> zero_density(mesh_, 0.0);
  field::CellField<double> reaction_rate(mesh_, 17.0);
  const std::array reactants{MassActionTerm{zero_density},
                             MassActionTerm{enormous_density, 2.0}};

  massActionReactionRate(
      reactants, RateCoefficientView{std::numeric_limits<double>::max()},
      reaction_rate);

  for (const double value : reaction_rate) {
    EXPECT_DOUBLE_EQ(value, 0.0);
  }
}

TEST_F(ReactionKineticsTest, InvalidLaterCellLeavesOutputUnchanged) {
  field::CellField<double> density(mesh_, 2.0);
  field::CellField<double> coefficient(mesh_, 4.0);
  field::CellField<double> reaction_rate(mesh_, 77.0);
  coefficient[mesh::CellId{1}] = -1.0;
  const std::array reactants{MassActionTerm{density}};

  EXPECT_THROW(massActionReactionRate(
                   reactants, RateCoefficientView{coefficient}, reaction_rate),
               std::invalid_argument);
  EXPECT_DOUBLE_EQ(reaction_rate[mesh::CellId{0}], 77.0);
  EXPECT_DOUBLE_EQ(reaction_rate[mesh::CellId{1}], 77.0);
}

TEST(MassActionUnitTest, DerivesCanonicalUnitsFromMolecularity) {
  EXPECT_EQ(canonicalCentimetreRateCoefficientUnit(1),
            units::precise::one / units::precise::s);
  EXPECT_EQ(canonicalCentimetreRateCoefficientUnit(2),
            units::precise::cm.pow(3) / units::precise::s);
  EXPECT_EQ(canonicalCentimetreRateCoefficientUnit(3),
            units::precise::cm.pow(6) / units::precise::s);
  EXPECT_THROW(static_cast<void>(canonicalCentimetreRateCoefficientUnit(0)),
               std::invalid_argument);
}

TEST_F(ReactionKineticsTest,
       AssemblerEvaluatesThreeBodyLawAndLeavesThirdBodySourceAbsent) {
  physics::SpeciesSet species;
  const auto electron = species.add({.name = "e", .charge = -1.0});
  const auto oxygen = species.add({.name = "O2"});
  const auto third_body = species.add({.name = "M"});
  const auto negative_ion = species.add({.name = "O2-", .charge = -1.0});
  ReactionNetwork network(species);
  const auto attachment = network.addReaction(
      {.name = "three-body attachment",
       .stoichiometry = {{electron, -1.0},
                         {oxygen, -1.0},
                         {negative_ion, +1.0}},
       .kinetic_orders = {{electron, 1.0}, {oxygen, 1.0}, {third_body, 1.0}}});

  const auto plasma_metadata = field::centimetrePlasmaFieldMetadata();
  physics::SpeciesCellFields density(mesh_, species.size(), 1.0,
                                     plasma_metadata.number_density);
  density[electron].fill(2.0);
  density[oxygen].fill(3.0);
  density[third_body].fill(5.0);
  field::FieldMetadata coefficient_metadata{
      .name = "three-body coefficient",
      .physical_quantity = unit::PhysicalQuantityMetadata{
          unit::QuantityKind::reaction_rate_coefficient,
          canonicalCentimetreRateCoefficientUnit(3)}};
  std::vector<field::CellField<double>> coefficients;
  coefficients.emplace_back(mesh_, 7.0, coefficient_metadata);
  ReactionRateFields rates(mesh_, network.size(), 0.0,
                           plasma_metadata.reaction_rate);
  MassActionReactionAssembler assembler(network, density);

  assembler.evaluate(coefficients, rates);

  for (const double value : rates[attachment]) {
    EXPECT_DOUBLE_EQ(value, 210.0);
  }
  physics::SpeciesCellFields source(mesh_, species.size(), 0.0,
                                    plasma_metadata.number_density_source);
  network.accumulateSources(rates.span(), source);
  for (const double value : source[third_body]) {
    EXPECT_DOUBLE_EQ(value, 0.0);
  }
}

TEST_F(ReactionKineticsTest,
       AssemblerRejectsBinaryCoefficientUnitForThreeBodyLaw) {
  physics::SpeciesSet species;
  const auto first = species.add({.name = "a"});
  const auto second = species.add({.name = "b"});
  const auto third = species.add({.name = "M"});
  const auto product = species.add({.name = "c"});
  ReactionNetwork network(species);
  static_cast<void>(network.addReaction(
      {.name = "three body",
       .stoichiometry = {{first, -1.0}, {second, -1.0}, {product, +1.0}},
       .kinetic_orders = {{first, 1.0}, {second, 1.0}, {third, 1.0}}}));
  const auto plasma_metadata = field::centimetrePlasmaFieldMetadata();
  physics::SpeciesCellFields density(mesh_, species.size(), 2.0,
                                     plasma_metadata.number_density);
  field::FieldMetadata wrong_coefficient_metadata{
      .name = "wrong binary coefficient",
      .physical_quantity = unit::PhysicalQuantityMetadata{
          unit::QuantityKind::reaction_rate_coefficient,
          canonicalCentimetreRateCoefficientUnit(2)}};
  std::vector<field::CellField<double>> coefficients;
  coefficients.emplace_back(mesh_, 1.0, wrong_coefficient_metadata);
  ReactionRateFields rates(mesh_, network.size(), 19.0,
                           plasma_metadata.reaction_rate);
  MassActionReactionAssembler assembler(network, density);

  EXPECT_THROW(assembler.evaluate(coefficients, rates), std::invalid_argument);
  for (const double value : rates[ReactionId{0}]) {
    EXPECT_DOUBLE_EQ(value, 19.0);
  }
}

TEST_F(ReactionKineticsTest,
       InvalidLaterReactionLeavesEveryRateFieldUnchanged) {
  physics::SpeciesSet species;
  const auto reactant = species.add({.name = "a"});
  const auto product = species.add({.name = "b"});
  ReactionNetwork network(species);
  static_cast<void>(
      network.addReaction({.name = "first",
                           .stoichiometry = {{reactant, -1.0}, {product, +1.0}},
                           .kinetic_orders = {{reactant, 1.0}}}));
  static_cast<void>(
      network.addReaction({.name = "second",
                           .stoichiometry = {{reactant, -1.0}, {product, +1.0}},
                           .kinetic_orders = {{reactant, 1.0}}}));
  physics::SpeciesCellFields density(mesh_, species.size(), 2.0);
  std::vector<field::CellField<double>> coefficients;
  coefficients.emplace_back(mesh_, 3.0);
  coefficients.emplace_back(mesh_, 4.0);
  coefficients.back()[mesh::CellId{1}] = -1.0;
  ReactionRateFields rates(mesh_, network.size(), 29.0);
  MassActionReactionAssembler assembler(network, density,
                                        /*allow_unitless_raw_values=*/true);

  EXPECT_THROW(assembler.evaluate(coefficients, rates), std::invalid_argument);
  for (const auto& rate : rates) {
    for (const double value : rate) {
      EXPECT_DOUBLE_EQ(value, 29.0);
    }
  }
}

}  // namespace pemu::physics::reaction::test
