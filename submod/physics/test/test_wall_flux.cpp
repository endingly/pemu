#include <gtest/gtest.h>

#include <pemu/mesh/moab_mesh.hpp>
#include <pemu/physics/wall/flux_assembler.hpp>

#include <filesystem>
#include <limits>
#include <stdexcept>

namespace pemu::physics::wall::test {
namespace {

/** @brief Returns the two-cell mesh used by wall-flux tests. */
[[nodiscard]] std::filesystem::path testMeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
}

/** @brief Finds the first face belonging to one physical boundary. */
[[nodiscard]] mesh::FaceId findBoundaryFace(const mesh::IMesh& mesh,
                                            mesh::BoundaryId boundary) {
  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {
    if (mesh.isBoundary(face) && mesh.boundaryId(face) == boundary) {
      return face;
    }
  }
  throw std::runtime_error("test boundary face not found");
}

class WallFluxAssemblerTest : public ::testing::Test {
 protected:
  /** @brief Creates electron/ion fields and one secondary-emitting wall. */
  WallFluxAssemblerTest()
      : mesh_(testMeshPath().string()),
        electron_(species_.add({.name = "e", .charge = -1.0})),
        ion_(species_.add({.name = "ion", .charge = +1.0})),
        wall_face_(findBoundaryFace(mesh_, mesh::BoundaryId{1})) {}

  /** @brief Returns a valid wall model used by individual tests. */
  [[nodiscard]] WallBoundarySet makeWallBoundaries() const {
    WallBoundarySet boundaries;
    boundaries.set(mesh::BoundaryId{1},
                   {.losses = {{electron_, 2.0, 5.0}, {ion_, 3.0, 0.0}},
                    .secondary_emissions = {{ion_, electron_, 0.1, 4.0}}});
    return boundaries;
  }

  mesh::MoabMesh mesh_;
  SpeciesSet species_;
  SpeciesId electron_;
  SpeciesId ion_;
  mesh::FaceId wall_face_;
};

}  // namespace

TEST(WallBoundarySetTest, ReplacesDefinitionsByBoundaryId) {
  WallBoundarySet boundaries;
  boundaries.set(mesh::BoundaryId{1}, {.losses = {}});
  boundaries.set(mesh::BoundaryId{1}, {.losses = {{SpeciesId{0}, 2.0, 3.0}}});

  ASSERT_EQ(boundaries.size(), 1u);
  ASSERT_EQ(boundaries.at(mesh::BoundaryId{1}).losses.size(), 1u);
  EXPECT_DOUBLE_EQ(
      boundaries.at(mesh::BoundaryId{1}).losses.front().particle_loss_velocity,
      2.0);
}

TEST_F(WallFluxAssemblerTest,
       AssemblesPrimaryLossAndSecondaryParticleAndEnergyInflux) {
  WallFluxAssembler assembler(mesh_, species_, makeWallBoundaries());
  SpeciesCellFields density(mesh_, species_.size(), 0.0);
  density[electron_].fill(5.0);
  density[ion_].fill(7.0);

  assembler.evaluate(density);

  EXPECT_EQ(assembler.activeFaces()[electron_][wall_face_], std::uint8_t{1});
  EXPECT_EQ(assembler.activeFaces()[ion_][wall_face_], std::uint8_t{1});
  EXPECT_DOUBLE_EQ(assembler.particleLossVelocity()[electron_][wall_face_],
                   2.0);
  EXPECT_DOUBLE_EQ(assembler.particleLossVelocity()[ion_][wall_face_], 3.0);
  EXPECT_DOUBLE_EQ(assembler.particleInwardFlux()[electron_][wall_face_], 2.1);
  EXPECT_DOUBLE_EQ(assembler.energyInwardFlux()[electron_][wall_face_], 8.4);
  EXPECT_DOUBLE_EQ(assembler.particleNormalFlux(density, electron_, wall_face_),
                   7.9);

  field::CellField<double> electron_energy_density(mesh_, 11.0);
  EXPECT_DOUBLE_EQ(assembler.energyNormalFlux(electron_energy_density,
                                              electron_, wall_face_),
                   46.6);
}

TEST_F(WallFluxAssemblerTest, RequiresEvaluationBeforeCompleteFluxQueries) {
  WallFluxAssembler assembler(mesh_, species_, makeWallBoundaries());
  SpeciesCellFields density(mesh_, species_.size(), 1.0);
  field::CellField<double> energy_density(mesh_, 1.0);

  EXPECT_THROW(static_cast<void>(assembler.particleNormalFlux(
                   density, electron_, wall_face_)),
               std::logic_error);
  EXPECT_THROW(static_cast<void>(assembler.energyNormalFlux(
                   energy_density, electron_, wall_face_)),
               std::logic_error);
}

TEST_F(WallFluxAssemblerTest, LeavesFacesOutsideConfiguredBoundaryInactive) {
  WallFluxAssembler assembler(mesh_, species_, makeWallBoundaries());

  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {
    if (mesh_.isBoundary(face) &&
        mesh_.boundaryId(face) != mesh::BoundaryId{1}) {
      EXPECT_EQ(assembler.activeFaces()[electron_][face], std::uint8_t{0});
      EXPECT_DOUBLE_EQ(assembler.particleLossVelocity()[electron_][face], 0.0);
    }
  }
}

TEST_F(WallFluxAssemblerTest, InvalidDensityLeavesPublishedInfluxUnchanged) {
  WallFluxAssembler assembler(mesh_, species_, makeWallBoundaries());
  SpeciesCellFields density(mesh_, species_.size(), 1.0);
  assembler.evaluate(density);
  const double previous_particle =
      assembler.particleInwardFlux()[electron_][wall_face_];
  const double previous_energy =
      assembler.energyInwardFlux()[electron_][wall_face_];
  density[ion_][mesh_.owner(wall_face_)] =
      std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(assembler.evaluate(density), std::invalid_argument);
  EXPECT_DOUBLE_EQ(assembler.particleInwardFlux()[electron_][wall_face_],
                   previous_particle);
  EXPECT_DOUBLE_EQ(assembler.energyInwardFlux()[electron_][wall_face_],
                   previous_energy);
  density[ion_][mesh_.owner(wall_face_)] = 2.0;
  EXPECT_THROW(static_cast<void>(assembler.particleNormalFlux(
                   density, electron_, wall_face_)),
               std::logic_error);
}

TEST_F(WallFluxAssemblerTest, RejectsMissingSpeciesLossForEmissionChannel) {
  WallBoundarySet boundaries;
  boundaries.set(mesh::BoundaryId{1},
                 {.losses = {{ion_, 3.0, 0.0}},
                  .secondary_emissions = {{ion_, electron_, 0.1, 4.0}}});

  EXPECT_THROW(WallFluxAssembler(mesh_, species_, std::move(boundaries)),
               std::invalid_argument);
}

TEST_F(WallFluxAssemblerTest, RejectsUnknownMeshBoundary) {
  WallBoundarySet boundaries;
  boundaries.set(mesh::BoundaryId{999}, {.losses = {{electron_, 2.0, 5.0}}});

  EXPECT_THROW(WallFluxAssembler(mesh_, species_, std::move(boundaries)),
               std::invalid_argument);
}

}  // namespace pemu::physics::wall::test
