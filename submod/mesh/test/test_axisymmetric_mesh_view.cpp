#include <gtest/gtest.h>

#include <pemu/mesh/axisymmetric_mesh_view.hpp>
#include <pemu/mesh/moab_mesh.hpp>

#include <cmath>
#include <filesystem>
#include <numbers>
#include <stdexcept>

namespace pemu::mesh::test {

namespace {

constexpr double kTolerance = 1e-12;

// Returns the shared two-cell mesh independently of CTest's working directory.
[[nodiscard]] std::filesystem::path testMeshPath() {
#ifdef PEMU_MESH_TEST_DATA_DIR
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
#else
  return std::filesystem::path{"submod/mesh/test/two_quads.msh"};
#endif
}

// Locates a face by its meridional center without assuming backend face IDs.
[[nodiscard]] FaceId findFace(const IMesh& mesh, const double radius,
                              const double axial_coordinate) {
  for (FaceId face = 0; face < mesh.numFaces(); ++face) {
    const Vec3 center = mesh.faceCenter(face);
    if (std::abs(center.x - radius) < kTolerance &&
        std::abs(center.y - axial_coordinate) < kTolerance) {
      return face;
    }
  }
  throw std::runtime_error("face not found");
}

class AxisymmetricMeshViewTest : public ::testing::Test {
 protected:
  AxisymmetricMeshViewTest()
      : planar_mesh_(testMeshPath().string()), mesh_(planar_mesh_) {}

  MoabMesh planar_mesh_;
  AxisymmetricMeshView mesh_;
};

}  // namespace

TEST_F(AxisymmetricMeshViewTest, PreservesMeridionalTopologyAndCoordinates) {
  EXPECT_EQ(mesh_.dimension(), planar_mesh_.dimension());
  EXPECT_EQ(mesh_.numCells(), planar_mesh_.numCells());
  EXPECT_EQ(mesh_.numFaces(), planar_mesh_.numFaces());
  EXPECT_EQ(mesh_.numVertices(), planar_mesh_.numVertices());

  for (FaceId face = 0; face < mesh_.numFaces(); ++face) {
    EXPECT_EQ(mesh_.owner(face), planar_mesh_.owner(face));
    EXPECT_EQ(mesh_.neighbor(face), planar_mesh_.neighbor(face));
    EXPECT_EQ(mesh_.isBoundary(face), planar_mesh_.isBoundary(face));
    EXPECT_EQ(mesh_.boundaryId(face), planar_mesh_.boundaryId(face));
    EXPECT_DOUBLE_EQ(mesh_.faceCenter(face).x, planar_mesh_.faceCenter(face).x);
    EXPECT_DOUBLE_EQ(mesh_.faceNormal(face).x, planar_mesh_.faceNormal(face).x);
  }
}

TEST_F(AxisymmetricMeshViewTest, RevolvesPlanarCellsIntoAnnularVolumes) {
  double total_volume = 0.0;
  for (CellId cell = 0; cell < mesh_.numCells(); ++cell) {
    const double radius = planar_mesh_.cellCenter(cell).x;
    const double expected =
        2.0 * std::numbers::pi * radius * planar_mesh_.cellVolume(cell);
    EXPECT_NEAR(mesh_.cellVolume(cell), expected, kTolerance);
    total_volume += mesh_.cellVolume(cell);
  }

  // Revolving [0, 2] x [0, 1] produces a radius-two, unit-length cylinder.
  EXPECT_NEAR(total_volume, 4.0 * std::numbers::pi, kTolerance);
}

TEST_F(AxisymmetricMeshViewTest, RevolvesFacesIntoPhysicalSurfaceAreas) {
  const FaceId axis_face = findFace(mesh_, 0.0, 0.5);
  const FaceId internal_face = findFace(mesh_, 1.0, 0.5);
  const FaceId outer_face = findFace(mesh_, 2.0, 0.5);
  const FaceId inner_end_face = findFace(mesh_, 0.5, 0.0);
  const FaceId outer_end_face = findFace(mesh_, 1.5, 1.0);

  EXPECT_DOUBLE_EQ(mesh_.faceArea(axis_face), 0.0);
  EXPECT_NEAR(mesh_.faceArea(internal_face), 2.0 * std::numbers::pi,
              kTolerance);
  EXPECT_NEAR(mesh_.faceArea(outer_face), 4.0 * std::numbers::pi, kTolerance);
  EXPECT_NEAR(mesh_.faceArea(inner_end_face), std::numbers::pi, kTolerance);
  EXPECT_NEAR(mesh_.faceArea(outer_end_face), 3.0 * std::numbers::pi,
              kTolerance);
}

TEST_F(AxisymmetricMeshViewTest, RejectsOutOfRangeMeasureQueries) {
  EXPECT_THROW(
      (void)mesh_.cellVolume(static_cast<CellId>(mesh_.numCells())),
      std::out_of_range);
  EXPECT_THROW((void)mesh_.faceArea(static_cast<FaceId>(mesh_.numFaces())),
               std::out_of_range);
}

}  // namespace pemu::mesh::test
