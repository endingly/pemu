#include <gtest/gtest.h>

#include <pemu/discretization/operators/linear_boundary_flux.hpp>
#include <pemu/discretization/operators/scharfetter_gummel_flux.hpp>
#include <pemu/mesh/moab_mesh.hpp>

#include <filesystem>
#include <limits>
#include <stdexcept>

namespace pemu::discretization::operators::test {
namespace {

/** @brief Returns the two-cell operator test mesh. */
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

/** @brief Finds the sole internal face in the two-cell mesh. */
[[nodiscard]] mesh::FaceId findInternalFace(const mesh::IMesh& mesh) {
  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {
    if (!mesh.isBoundary(face)) {
      return face;
    }
  }
  throw std::runtime_error("test internal face not found");
}

class LinearBoundaryFluxTest : public ::testing::Test {
 protected:
  /** @brief Loads a mesh and selects its left boundary. */
  LinearBoundaryFluxTest()
      : mesh_(testMeshPath().string()),
        wall_face_(findBoundaryFace(mesh_, mesh::BoundaryId{1})) {}

  /** @brief Creates density-matching Dirichlet values off the wall. */
  [[nodiscard]] boundary::BoundaryConditionSet remainingBoundaries() const {
    boundary::BoundaryConditionSet conditions;
    conditions.setDirichlet(mesh::BoundaryId{2}, 3.0);
    conditions.setDirichlet(mesh::BoundaryId{3}, 3.0);
    conditions.setDirichlet(mesh::BoundaryId{4}, 3.0);
    return conditions;
  }

  mesh::MoabMesh mesh_;
  mesh::FaceId wall_face_;
};

}  // namespace

TEST_F(LinearBoundaryFluxTest, OverridesOnlySelectedSgBoundaryFace) {
  field::CellField<double> state(mesh_, 3.0);
  field::FaceField<double> velocity(mesh_, 0.0);
  field::FaceField<double> flux(mesh_, 91.0);
  field::FaceField<std::uint8_t> active(mesh_, std::uint8_t{0});
  field::FaceField<double> loss_velocity(mesh_, 0.0);
  field::FaceField<double> inward_flux(mesh_, 0.0);
  active[wall_face_] = std::uint8_t{1};
  loss_velocity[wall_face_] = 2.0;
  inward_flux[wall_face_] = 1.0;
  const LinearBoundaryFluxView wall_flux(active, loss_velocity, inward_flux);
  const auto conditions = remainingBoundaries();

  scharfetterGummelFlux(state, velocity, 1.0, conditions, wall_flux, flux);

  EXPECT_DOUBLE_EQ(flux[wall_face_], 5.0);
  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {
    if (face != wall_face_) {
      EXPECT_NEAR(flux[face], 0.0, 1e-12);
    }
  }
}

TEST_F(LinearBoundaryFluxTest, AddsWallLossToExplicitDiagonal) {
  field::FaceField<std::uint8_t> active(mesh_, std::uint8_t{0});
  field::FaceField<double> loss_velocity(mesh_, 0.0);
  field::FaceField<double> inward_flux(mesh_, 0.0);
  field::CellField<double> loss_rate(mesh_, 0.0);
  active[wall_face_] = std::uint8_t{1};
  loss_velocity[wall_face_] = 2.0;
  const LinearBoundaryFluxView wall_flux(active, loss_velocity, inward_flux);

  computeLinearBoundaryLossRate(wall_flux, loss_rate);

  const auto owner = mesh_.owner(wall_face_);
  const double expected =
      2.0 * mesh_.faceArea(wall_face_) / mesh_.cellVolume(owner);
  EXPECT_DOUBLE_EQ(loss_rate[owner], expected);
}

TEST_F(LinearBoundaryFluxTest, AggregatesSeveralWallFacesOnOneOwnerCell) {
  field::FaceField<std::uint8_t> active(mesh_, std::uint8_t{0});
  field::FaceField<double> loss_velocity(mesh_, 0.0);
  field::FaceField<double> inward_flux(mesh_, 0.0);
  field::CellField<double> loss_rate(mesh_, 91.0);
  const auto owner = mesh_.owner(wall_face_);
  double expected = 0.0;
  std::size_t active_count = 0;
  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {
    if (mesh_.isBoundary(face) && mesh_.owner(face) == owner) {
      active[face] = std::uint8_t{1};
      loss_velocity[face] = 2.0 + static_cast<double>(active_count);
      expected +=
          loss_velocity[face] * mesh_.faceArea(face) / mesh_.cellVolume(owner);
      ++active_count;
    }
  }
  ASSERT_GE(active_count, 2u);
  const LinearBoundaryFluxView wall_flux(active, loss_velocity, inward_flux);

  computeLinearBoundaryLossRate(wall_flux, loss_rate);

  EXPECT_DOUBLE_EQ(loss_rate[owner], expected);
}

TEST_F(LinearBoundaryFluxTest, RejectsAggregateLossRateOverflow) {
  field::FaceField<std::uint8_t> active(mesh_, std::uint8_t{0});
  field::FaceField<double> loss_velocity(mesh_, 0.0);
  field::FaceField<double> inward_flux(mesh_, 0.0);
  field::CellField<double> loss_rate(mesh_, 0.0);
  const auto owner = mesh_.owner(wall_face_);
  const double desired_contribution = 0.75 * std::numeric_limits<double>::max();
  std::size_t active_count = 0;
  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {
    if (!mesh_.isBoundary(face) || mesh_.owner(face) != owner) {
      continue;
    }
    const long double coefficient =
        static_cast<long double>(desired_contribution) *
        static_cast<long double>(mesh_.cellVolume(owner)) /
        static_cast<long double>(mesh_.faceArea(face));
    ASSERT_LE(coefficient,
              static_cast<long double>(std::numeric_limits<double>::max()));
    active[face] = std::uint8_t{1};
    loss_velocity[face] = static_cast<double>(coefficient);
    ++active_count;
    if (active_count == 2) {
      break;
    }
  }
  ASSERT_EQ(active_count, 2u);
  const LinearBoundaryFluxView wall_flux(active, loss_velocity, inward_flux);

  EXPECT_THROW(computeLinearBoundaryLossRate(wall_flux, loss_rate),
               std::overflow_error);
}

TEST_F(LinearBoundaryFluxTest, InvalidInternalMaskLeavesSgOutputUnchanged) {
  field::CellField<double> state(mesh_, 3.0);
  field::FaceField<double> velocity(mesh_, 0.0);
  field::FaceField<double> flux(mesh_, 91.0);
  field::FaceField<std::uint8_t> active(mesh_, std::uint8_t{0});
  field::FaceField<double> loss_velocity(mesh_, 0.0);
  field::FaceField<double> inward_flux(mesh_, 0.0);
  active[findInternalFace(mesh_)] = std::uint8_t{1};
  const LinearBoundaryFluxView wall_flux(active, loss_velocity, inward_flux);
  const auto conditions = remainingBoundaries();

  EXPECT_THROW(
      scharfetterGummelFlux(state, velocity, 1.0, conditions, wall_flux, flux),
      std::invalid_argument);
  for (const double value : flux) {
    EXPECT_DOUBLE_EQ(value, 91.0);
  }
}

}  // namespace pemu::discretization::operators::test
