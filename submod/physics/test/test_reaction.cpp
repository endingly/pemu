#include <gtest/gtest.h>
#include <filesystem>
#include <pemu/mesh/moab_mesh.hpp>
#include <pemu/physics/ionization_reaction.hpp>
#include <pemu/physics/reaction.hpp>

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

  physics::reaction::electronImpactIonizationRate(
      electron_density, neutral_density, rate_coefficient, rate);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(rate[cell], 24.0, 1e-12);
  }
}

TEST_F(ReactionTest, ZeroElectronDensityProducesNoIonization) {
  field::CellField<double> electron_density(mesh_, 0.0);

  field::CellField<double> rate(mesh_, 123.0);

  physics::reaction::electronImpactIonizationRate(electron_density, 10.0, 2.0,
                                                  rate);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(rate[cell], 0.0, 1e-12);
  }
}

TEST_F(ReactionTest, IonizationCreatesElectronIonPairs) {
  field::CellField<double> rate(mesh_, 5.0);

  field::CellField<double> electron_source(mesh_, 0.0);

  field::CellField<double> ion_source(mesh_, 0.0);

  physics::reaction::addPairProductionSource(rate, electron_source, ion_source);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(electron_source[cell], 5.0, 1e-12);

    EXPECT_NEAR(ion_source[cell], 5.0, 1e-12);
  }
}

TEST_F(ReactionTest, PairProductionAccumulatesIntoExistingSource) {
  field::CellField<double> rate(mesh_, 5.0);

  field::CellField<double> electron_source(mesh_, 2.0);

  field::CellField<double> ion_source(mesh_, 3.0);

  physics::reaction::addPairProductionSource(rate, electron_source, ion_source);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(electron_source[cell], 7.0, 1e-12);

    EXPECT_NEAR(ion_source[cell], 8.0, 1e-12);
  }
}

TEST_F(ReactionTest, PairProductionCreatesNoNetCharge) {
  field::CellField<double> rate(mesh_, 5.0);

  field::CellField<double> electron_source(mesh_, 0.0);

  field::CellField<double> ion_source(mesh_, 0.0);

  physics::reaction::addPairProductionSource(rate, electron_source, ion_source);

  constexpr double electron_charge = -2.0;

  constexpr double ion_charge = +2.0;

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    const double charge_source =
        electron_charge * electron_source[cell] + ion_charge * ion_source[cell];

    EXPECT_NEAR(charge_source, 0.0, 1e-12);
  }
}

TEST_F(ReactionNetworkTest, PairIonizationConservesCharge) {
  physics::SpeciesSet species;

  const auto electron = species.add({.name = "e", .charge = -2.0});

  const auto ion = species.add({.name = "ion", .charge = +2.0});

  physics::ReactionNetwork network(species);

  const auto reaction =
      network.addReaction({.name = "ionization",

                           .stoichiometry = {{electron, +1.0}, {ion, +1.0}}});

  EXPECT_NEAR(network.netChargePerReaction(reaction), 0.0, 1e-12);

  EXPECT_TRUE(network.conservesCharge(reaction));
}

TEST_F(ReactionNetworkTest, DetectsChargeViolatingReaction) {
  physics::SpeciesSet species;

  const auto electron = species.add({.name = "e", .charge = -1.0});

  physics::ReactionNetwork network(species);

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

  physics::ReactionNetwork network(species);

  auto _ = network.addReaction(
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

  physics::ReactionNetwork network(species);

  auto _ = network.addReaction(
      {.name = "ionization", .stoichiometry = {{electron, +1.0}, {ion, +1.0}}});

  auto _ =
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