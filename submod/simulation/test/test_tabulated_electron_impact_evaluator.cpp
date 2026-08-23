#include <gtest/gtest.h>

#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>
#include <pemu/field/plasma_field_metadata.hpp>
#include <pemu/mesh/moab_mesh.hpp>
#include <pemu/physics/reaction/network.hpp>
#include <pemu/physics/species.hpp>
#include <pemu/physics/reaction/rate_coefficient.hpp>
#include <pemu/simulation/plasma_reaction_rate_evaluator.hpp>
#include <pemu/simulation/tabulated_electron_impact_evaluator.hpp>

#include <filesystem>
#include <utility>

namespace pemu::simulation::test {
namespace {

/** @brief Returns the two-cell mesh used by evaluator contract tests. */
[[nodiscard]] std::filesystem::path testMeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
}

class TabulatedElectronImpactEvaluatorTest : public ::testing::Test {
 protected:
  /** @brief Creates synchronized fields for one electron and one reaction. */
  TabulatedElectronImpactEvaluatorTest()
      : mesh_(testMeshPath().string()),
        density_(mesh_, 1, 2.0e6),
        potential_(mesh_, 0.0),
        electric_field_normal_(mesh_, 0.0),
        mean_energy_(mesh_, 3.0),
        target_density_(mesh_, 2.5e19),
        reaction_rates_(mesh_, 1, 0.0) {}

  /** @brief Returns a context bound to the fixture's synchronized state. */
  [[nodiscard]] PlasmaReactionRateContext context() const {
    return {.density = density_,
            .potential = potential_,
            .electric_field_normal = electric_field_normal_,
            .electron_mean_energy = mean_energy_};
  }

  mesh::MoabMesh mesh_;
  physics::SpeciesCellFields density_;
  field::CellField<double> potential_;
  field::FaceField<double> electric_field_normal_;
  field::CellField<double> mean_energy_;
  field::CellField<double> target_density_;
  physics::reaction::ReactionRateFields reaction_rates_;
};

}  // namespace

TEST_F(TabulatedElectronImpactEvaluatorTest,
       EvaluatesReducedFieldDependentReactionThroughTypeErasure) {
  const auto metadata = field::centimetrePlasmaFieldMetadata();
  density_[physics::SpeciesId{0}].setMetadata(metadata.number_density);
  target_density_.setMetadata(metadata.number_density);
  electric_field_normal_.setMetadata(metadata.electric_field);
  reaction_rates_[physics::reaction::ReactionId{0}].setMetadata(metadata.reaction_rate);
  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {
    electric_field_normal_[face] = 400.0 * mesh_.faceNormal(face).x;
  }
  TabulatedElectronImpactEvaluator concrete(
      physics::SpeciesId{0}, physics::reaction::ReactionId{0}, target_density_,
      ElectronImpactRateCoordinate::reduced_electric_field_townsend,
      physics::reaction::TabulatedRateCoefficient(
          {1.0, 2.0}, {1.0e-13, 2.0e-13},
          physics::reaction::RateInterpolation::linear));
  PlasmaReactionRateEvaluator evaluator = std::move(concrete);

  evaluator(context(), reaction_rates_);

  for (const double rate : reaction_rates_[physics::reaction::ReactionId{0}]) {
    EXPECT_NEAR(rate, 8.0e12, 2.0);
  }
}

TEST_F(TabulatedElectronImpactEvaluatorTest,
       EvaluatesElectronTemperatureDependentReaction) {
  TabulatedElectronImpactEvaluator evaluator(
      physics::SpeciesId{0}, physics::reaction::ReactionId{0}, target_density_,
      ElectronImpactRateCoordinate::electron_temperature_ev,
      physics::reaction::TabulatedRateCoefficient(
          {1.0, 3.0}, {1.0e-13, 3.0e-13},
          physics::reaction::RateInterpolation::linear),
      true);

  evaluator(context(), reaction_rates_);

  for (const double rate : reaction_rates_[physics::reaction::ReactionId{0}]) {
    EXPECT_NEAR(rate, 1.0e13, 2.0);
  }
}

TEST_F(TabulatedElectronImpactEvaluatorTest,
       RejectsUnitlessFieldsWithoutExplicitOptIn) {
  TabulatedElectronImpactEvaluator evaluator(
      physics::SpeciesId{0}, physics::reaction::ReactionId{0}, target_density_,
      ElectronImpactRateCoordinate::electron_temperature_ev,
      physics::reaction::TabulatedRateCoefficient(
          {1.0, 3.0}, {1.0e-13, 3.0e-13},
          physics::reaction::RateInterpolation::linear));

  EXPECT_THROW(evaluator(context(), reaction_rates_), std::invalid_argument);
}

}  // namespace pemu::simulation::test
