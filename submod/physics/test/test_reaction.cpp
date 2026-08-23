#include <gtest/gtest.h>
#include <filesystem>
#include <pemu/mesh/moab_mesh.hpp>
#include <pemu/physics/reaction/mass_action.hpp>
#include <pemu/physics/reaction/network.hpp>

#include <algorithm>
#include <limits>

namespace pemu::physics::test {

namespace {

std::filesystem::path testMeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
}

class ReactionTest : public ::testing::Test {
 protected:
  ReactionTest() : mesh_(testMeshPath().string()) {}

  mesh::MoabMesh mesh_;
};

class ReactionNetworkTest : public ::testing::Test {
 protected:
  ReactionNetworkTest() : mesh_(testMeshPath().string()) {}
  mesh::MoabMesh mesh_;
};

};  // namespace

TEST_F(ReactionTest, ElectronImpactIonizationComputesExpectedRate) {
  field::CellField<double> electron_density(mesh_, 3.0);

  field::CellField<double> rate(mesh_, 0.0);

  constexpr double neutral_density = 4.0;

  constexpr double rate_coefficient = 2.0;

  physics::reaction::binaryReactionRate(
      electron_density, neutral_density, rate_coefficient, rate);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(rate[cell], 24.0, 1e-12);
  }
}

TEST_F(ReactionTest, ZeroElectronDensityProducesNoIonization) {
  field::CellField<double> electron_density(mesh_, 0.0);

  field::CellField<double> rate(mesh_, 123.0);

  physics::reaction::binaryReactionRate(electron_density, 10.0, 2.0, rate);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(rate[cell], 0.0, 1e-12);
  }
}

TEST_F(ReactionNetworkTest, PairIonizationConservesCharge) {
  physics::SpeciesSet species;

  const auto electron = species.add({.name = "e", .charge = -2.0});

  const auto ion = species.add({.name = "ion", .charge = +2.0});

  physics::reaction::ReactionNetwork network(species);

  const auto reaction =
      network.addReaction({.name = "ionization",

                           .stoichiometry = {{electron, +1.0}, {ion, +1.0}}});

  EXPECT_NEAR(network.netChargePerReaction(reaction), 0.0, 1e-12);

  EXPECT_TRUE(network.conservesCharge(reaction));
}

TEST_F(ReactionNetworkTest,
       StoresThirdBodyOrderIndependentlyFromNetStoichiometry) {
  physics::SpeciesSet species;
  const auto electron = species.add({.name = "e", .charge = -1.0});
  const auto oxygen = species.add({.name = "O2"});
  const auto third_body = species.add({.name = "M"});
  const auto negative_ion = species.add({.name = "O2-", .charge = -1.0});
  physics::reaction::ReactionNetwork network(species);

  const auto attachment = network.addReaction(
      {.name = "three-body attachment",
       .stoichiometry = {{electron, -1.0},
                         {oxygen, -1.0},
                         {negative_ion, +1.0}},
       .kinetic_orders = {
           {electron, 1.0}, {oxygen, 1.0}, {third_body, 1.0}}});

  const auto& definition = network.at(attachment);
  ASSERT_EQ(definition.kinetic_orders.size(), 3u);
  EXPECT_EQ(definition.kinetic_orders[2].species, third_body);
  EXPECT_TRUE(std::ranges::none_of(
      definition.stoichiometry,
      [third_body](const auto& term) { return term.species == third_body; }));
  EXPECT_TRUE(network.conservesCharge(attachment));
}

TEST_F(ReactionNetworkTest, RejectsInvalidAndDuplicateKineticOrders) {
  physics::SpeciesSet species;
  const auto electron = species.add({.name = "e", .charge = -1.0});
  physics::reaction::ReactionNetwork network(species);

  EXPECT_THROW(static_cast<void>(network.addReaction(
                   {.name = "duplicate",
                    .stoichiometry = {{electron, -1.0}},
                    .kinetic_orders = {{electron, 1.0}, {electron, 2.0}}})),
               std::invalid_argument);
  EXPECT_THROW(
      static_cast<void>(network.addReaction(
          {.name = "invalid order",
           .stoichiometry = {{electron, -1.0}},
           .kinetic_orders = {
               {electron, std::numeric_limits<double>::infinity()}}})),
      std::invalid_argument);
}

TEST_F(ReactionNetworkTest, DetectsChargeViolatingReaction) {
  physics::SpeciesSet species;

  const auto electron = species.add({.name = "e", .charge = -1.0});

  physics::reaction::ReactionNetwork network(species);

  const auto invalid =
      network.addReaction({.name = "electron from nowhere",

                           .stoichiometry = {{electron, +1.0}}});

  EXPECT_FALSE(network.conservesCharge(invalid));

  EXPECT_NEAR(network.netChargePerReaction(invalid), -1.0, 1e-12);
}

TEST_F(ReactionNetworkTest, AccumulatesStoichiometricSources) {
  physics::SpeciesSet species;

  const auto electron = species.add({.name = "e", .charge = -1.0});

  const auto ion = species.add({.name = "ion", .charge = +1.0});

  physics::reaction::ReactionNetwork network(species);

  [[maybe_unused]] const auto ionization = network.addReaction(
      {.name = "ionization", .stoichiometry = {{electron, +1.0}, {ion, +1.0}}});

  std::vector<field::CellField<double>> rates;

  rates.emplace_back(mesh_, 5.0);

  physics::SpeciesCellFields source(mesh_, species.size(), 0.0);

  network.accumulateSources(rates, source);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(source[electron][cell], 5.0, 1e-12);

    EXPECT_NEAR(source[ion][cell], 5.0, 1e-12);
  }
}

TEST_F(ReactionNetworkTest, MultipleReactionsAccumulateCorrectly) {
  physics::SpeciesSet species;

  const auto electron = species.add({.name = "e", .charge = -1.0});

  const auto ion = species.add({.name = "ion", .charge = +1.0});

  physics::reaction::ReactionNetwork network(species);

  [[maybe_unused]] const auto ionization = network.addReaction(
      {.name = "ionization", .stoichiometry = {{electron, +1.0}, {ion, +1.0}}});

  [[maybe_unused]] const auto recombination =
      network.addReaction({.name = "recombination",
                           .stoichiometry = {{electron, -1.0}, {ion, -1.0}}});

  std::vector<field::CellField<double>> rates;

  rates.emplace_back(mesh_, 5.0);

  rates.emplace_back(mesh_, 2.0);

  physics::SpeciesCellFields source(mesh_, species.size(), 0.0);

  network.accumulateSources(rates, source);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(source[electron][cell], 3.0, 1e-12);

    EXPECT_NEAR(source[ion][cell], 3.0, 1e-12);
  }
}

};  // namespace pemu::physics::test
