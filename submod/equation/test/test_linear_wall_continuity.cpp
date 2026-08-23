#include <gtest/gtest.h>

#include <pemu/equation/adaptive_step_explicit_species_continuity_stepper.hpp>
#include <pemu/mesh/moab_mesh.hpp>

#include <filesystem>
#include <stdexcept>

namespace pemu::equation::test {
namespace {

/** @brief Returns the two-cell continuity test mesh. */
[[nodiscard]] std::filesystem::path testMeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
}

/** @brief Finds one face by physical boundary ID. */
[[nodiscard]] mesh::FaceId findBoundaryFace(const mesh::IMesh& mesh,
                                            mesh::BoundaryId boundary) {
  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {
    if (mesh.isBoundary(face) && mesh.boundaryId(face) == boundary) {
      return face;
    }
  }
  throw std::runtime_error("test boundary face not found");
}

/** @brief Integrates a cell number-density field over the mesh. */
[[nodiscard]] double integratedState(const mesh::IMesh& mesh,
                                     const field::CellField<double>& state) {
  double result = 0.0;
  for (mesh::CellId cell = 0; cell < mesh.numCells(); ++cell) {
    result += state[cell] * mesh.cellVolume(cell);
  }
  return result;
}

class LinearWallContinuityTest : public ::testing::Test {
 protected:
  /** @brief Builds a zero-drift wall transport fixture. */
  LinearWallContinuityTest()
      : mesh_(testMeshPath().string()),
        wall_face_(findBoundaryFace(mesh_, mesh::BoundaryId{1})),
        velocity_(mesh_, 0.0),
        active_(mesh_, std::uint8_t{0}),
        loss_velocity_(mesh_, 0.0),
        inward_flux_(mesh_, 0.0) {
    active_[wall_face_] = std::uint8_t{1};
    conditions_.setDirichlet(mesh::BoundaryId{2}, 3.0);
    conditions_.setDirichlet(mesh::BoundaryId{3}, 3.0);
    conditions_.setDirichlet(mesh::BoundaryId{4}, 3.0);
  }

  mesh::MoabMesh mesh_;
  mesh::FaceId wall_face_;
  field::FaceField<double> velocity_;
  field::FaceField<std::uint8_t> active_;
  field::FaceField<double> loss_velocity_;
  field::FaceField<double> inward_flux_;
  boundary::BoundaryConditionSet conditions_;
};

}  // namespace

TEST_F(LinearWallContinuityTest, WallFluxChangesIntegratedParticlesExactly) {
  loss_velocity_[wall_face_] = 2.0;
  inward_flux_[wall_face_] = 1.0;
  const discretization::operators::LinearBoundaryFluxView wall_flux(
      active_, loss_velocity_, inward_flux_);
  AdaptiveStepExplicitSpeciesContinuityStepper stepper(mesh_, velocity_, 1.0,
                                                       conditions_, wall_flux);
  field::CellField<double> density(mesh_, 3.0);
  field::CellField<double> source(mesh_, 0.0);
  field::CellField<double> increment(mesh_, 0.0);
  constexpr double dt = 0.1;
  const double before = integratedState(mesh_, density);

  stepper.computeIncrement(density, source, dt, increment);
  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {
    density[cell] += increment[cell];
  }

  const double expected_change =
      -dt * (2.0 * 3.0 - 1.0) * mesh_.faceArea(wall_face_);
  EXPECT_NEAR(integratedState(mesh_, density) - before, expected_change, 1e-12);
}

TEST_F(LinearWallContinuityTest, InwardEmissionProducesParticles) {
  loss_velocity_[wall_face_] = 2.0;
  inward_flux_[wall_face_] = 7.0;
  const discretization::operators::LinearBoundaryFluxView wall_flux(
      active_, loss_velocity_, inward_flux_);
  AdaptiveStepExplicitSpeciesContinuityStepper stepper(mesh_, velocity_, 1.0,
                                                       conditions_, wall_flux);
  field::CellField<double> density(mesh_, 3.0);
  field::CellField<double> source(mesh_, 0.0);
  field::CellField<double> increment(mesh_, 0.0);
  const double before = integratedState(mesh_, density);

  stepper.computeIncrement(density, source, 0.1, increment);
  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {
    density[cell] += increment[cell];
  }

  EXPECT_GT(integratedState(mesh_, density), before);
}

}  // namespace pemu::equation::test
