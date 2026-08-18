#include <gtest/gtest.h>

#include "moab_mesh.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <span>

namespace pemu::mesh::test {

namespace {

constexpr double kTolerance = 1e-12;

double dot(const Vec3& a, const Vec3& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

double norm(const Vec3& v) {
  return std::sqrt(dot(v, v));
}

bool contains(std::span<const FaceId> faces, FaceId target) {
  return std::find(faces.begin(), faces.end(), target) != faces.end();
}

// ------------------------------------------------------------
// Prefer passing test-data directory from CMake:
//
//   PEMU_MESH_TEST_DATA_DIR
//
// This makes CTest independent of current working directory.
// ------------------------------------------------------------

std::filesystem::path testMeshPath() {
#ifdef PEMU_MESH_TEST_DATA_DIR

  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";

#else

  return std::filesystem::path{"submod/mesh/test/two_quads.msh"};

#endif
}

class MoabMeshTest : public ::testing::Test {
 protected:
  MoabMeshTest() : mesh_(testMeshPath().string()) {}

  MoabMesh mesh_;
};

}  // namespace

// ============================================================
// Basic mesh information
// ============================================================

TEST_F(MoabMeshTest, LoadsTwoDimensionalMesh) {
  EXPECT_EQ(mesh_.dimension(), 2);
}

TEST_F(MoabMeshTest, HasExpectedEntityCounts) {
  EXPECT_EQ(mesh_.numCells(), 2u);

  EXPECT_EQ(mesh_.numVertices(), 6u);

  //
  // Two adjacent quads:
  //
  // 4 + 4 - 1 shared edge = 7 unique faces.
  //
  EXPECT_EQ(mesh_.numFaces(), 7u);
}

// ============================================================
// Boundary / internal face topology
// ============================================================

TEST_F(MoabMeshTest, HasSixBoundaryFacesAndOneInternalFace) {
  std::size_t boundary_count = 0;
  std::size_t internal_count = 0;

  for (FaceId f = 0; f < mesh_.numFaces(); ++f) {

    if (mesh_.isBoundary(f)) {

      ++boundary_count;

      EXPECT_NE(mesh_.owner(f), invalid_cell);

      EXPECT_EQ(mesh_.neighbor(f), invalid_cell);

    } else {

      ++internal_count;

      EXPECT_NE(mesh_.owner(f), invalid_cell);

      EXPECT_NE(mesh_.neighbor(f), invalid_cell);

      EXPECT_NE(mesh_.owner(f), mesh_.neighbor(f));
    }
  }

  EXPECT_EQ(boundary_count, 6u);

  EXPECT_EQ(internal_count, 1u);
}

// ============================================================
// Each quad should have exactly four faces.
// ============================================================

TEST_F(MoabMeshTest, EachCellHasFourFaces) {
  ASSERT_EQ(mesh_.numCells(), 2u);

  EXPECT_EQ(mesh_.cellFaces(0).size(), 4u);

  EXPECT_EQ(mesh_.cellFaces(1).size(), 4u);
}

// ============================================================
// Verify:
//
// if face.owner == cell
//
// or
//
// face.neighbor == cell
//
// then:
//
//     face ∈ cellFaces(cell)
//
// ============================================================

TEST_F(MoabMeshTest, CellFaceConnectivityIsConsistent) {
  for (FaceId f = 0; f < mesh_.numFaces(); ++f) {

    const CellId owner = mesh_.owner(f);

    ASSERT_NE(owner, invalid_cell);

    EXPECT_TRUE(contains(mesh_.cellFaces(owner), f));

    if (!mesh_.isBoundary(f)) {

      const CellId neighbor = mesh_.neighbor(f);

      ASSERT_NE(neighbor, invalid_cell);

      EXPECT_TRUE(contains(mesh_.cellFaces(neighbor), f));
    }
  }
}

// ============================================================
// The two cells must share exactly one face.
// ============================================================

TEST_F(MoabMeshTest, CellsShareExactlyOneFace) {
  const auto faces0 = mesh_.cellFaces(0);

  const auto faces1 = mesh_.cellFaces(1);

  std::size_t common_faces = 0;

  for (const FaceId f0 : faces0) {

    if (contains(faces1, f0)) {

      ++common_faces;
    }
  }

  EXPECT_EQ(common_faces, 1u);
}

// ============================================================
// Cell geometry
//
// C0 center = (0.5, 0.5)
// C1 center = (1.5, 0.5)
//
// Note:
//
// We intentionally do not assume MOAB assigned C0 to CellId=0
// based solely on input-file order in future backends.
//
// Therefore inspect centers rather than impose semantic identity.
// ============================================================

TEST_F(MoabMeshTest, ComputesCorrectCellCenters) {
  ASSERT_EQ(mesh_.numCells(), 2u);

  const Vec3 c0 = mesh_.cellCenter(0);

  const Vec3 c1 = mesh_.cellCenter(1);

  //
  // Both cells lie at y = 0.5.
  //
  EXPECT_NEAR(c0.y, 0.5, kTolerance);

  EXPECT_NEAR(c1.y, 0.5, kTolerance);

  EXPECT_NEAR(c0.z, 0.0, kTolerance);

  EXPECT_NEAR(c1.z, 0.0, kTolerance);

  //
  // Their x coordinates must be {0.5, 1.5}.
  //
  const double xmin = std::min(c0.x, c1.x);

  const double xmax = std::max(c0.x, c1.x);

  EXPECT_NEAR(xmin, 0.5, kTolerance);

  EXPECT_NEAR(xmax, 1.5, kTolerance);
}

TEST_F(MoabMeshTest, ComputesCorrectCellAreas) {
  ASSERT_EQ(mesh_.numCells(), 2u);

  for (CellId c = 0; c < mesh_.numCells(); ++c) {

    //
    // cellVolume() in a 2D mesh means geometric area.
    //
    EXPECT_NEAR(mesh_.cellVolume(c), 1.0, kTolerance);
  }
}

// ============================================================
// Face geometry
// ============================================================

TEST_F(MoabMeshTest, ComputesCorrectFaceLengths) {
  for (FaceId f = 0; f < mesh_.numFaces(); ++f) {

    //
    // Every edge in this mesh has unit length.
    //
    EXPECT_NEAR(mesh_.faceArea(f), 1.0, kTolerance);
  }
}

TEST_F(MoabMeshTest, FaceNormalsAreUnitVectors) {
  for (FaceId f = 0; f < mesh_.numFaces(); ++f) {

    const Vec3 n = mesh_.faceNormal(f);

    EXPECT_NEAR(norm(n), 1.0, kTolerance);
  }
}

// ============================================================
// Internal normal orientation invariant:
//
// normal points owner -> neighbor.
//
// Therefore:
//
//     n dot (x_neighbor - x_owner) > 0
//
// This invariant is extremely important for FVM.
// ============================================================

TEST_F(MoabMeshTest, InternalFaceNormalPointsFromOwnerToNeighbor) {
  std::size_t checked = 0;

  for (FaceId f = 0; f < mesh_.numFaces(); ++f) {

    if (mesh_.isBoundary(f)) {
      continue;
    }

    const CellId owner = mesh_.owner(f);

    const CellId neighbor = mesh_.neighbor(f);

    const Vec3 owner_center = mesh_.cellCenter(owner);

    const Vec3 neighbor_center = mesh_.cellCenter(neighbor);

    const Vec3 direction = neighbor_center - owner_center;

    const Vec3 normal = mesh_.faceNormal(f);

    EXPECT_GT(dot(normal, direction), 0.0);

    ++checked;
  }

  EXPECT_EQ(checked, 1u);
}

// ============================================================
// Boundary normal orientation:
//
// normal points from cell center toward boundary/outside.
//
// For this simple convex Cartesian mesh:
//
//     n dot (faceCenter - cellCenter) > 0
//
// ============================================================

TEST_F(MoabMeshTest, BoundaryFaceNormalsPointOutward) {
  std::size_t checked = 0;

  for (FaceId f = 0; f < mesh_.numFaces(); ++f) {

    if (!mesh_.isBoundary(f)) {
      continue;
    }

    const CellId owner = mesh_.owner(f);

    const Vec3 owner_center = mesh_.cellCenter(owner);

    const Vec3 face_center = mesh_.faceCenter(f);

    const Vec3 outward = face_center - owner_center;

    const Vec3 normal = mesh_.faceNormal(f);

    EXPECT_GT(dot(normal, outward), 0.0);

    ++checked;
  }

  EXPECT_EQ(checked, 6u);
}

// ============================================================
// Internal face geometry.
//
// There is exactly one internal face:
//
//     x = 1
//     y = 0.5
//
// ============================================================

TEST_F(MoabMeshTest, InternalFaceHasCorrectGeometry) {
  bool found = false;

  for (FaceId f = 0; f < mesh_.numFaces(); ++f) {

    if (mesh_.isBoundary(f)) {
      continue;
    }

    ASSERT_FALSE(found) << "Expected exactly one internal face";

    found = true;

    const Vec3 center = mesh_.faceCenter(f);

    EXPECT_NEAR(center.x, 1.0, kTolerance);

    EXPECT_NEAR(center.y, 0.5, kTolerance);

    EXPECT_NEAR(center.z, 0.0, kTolerance);

    EXPECT_NEAR(mesh_.faceArea(f), 1.0, kTolerance);
  }

  EXPECT_TRUE(found);
}

// ============================================================
// Boundary metadata.
//
// Current MoabMesh v0.1:
//
// boundary face -> BoundaryId = 0
//
// internal face -> invalid_boundary
//
// ============================================================

TEST_F(MoabMeshTest, AssignsDefaultBoundaryId) {
  for (FaceId f = 0; f < mesh_.numFaces(); ++f) {

    if (mesh_.isBoundary(f)) {

      EXPECT_EQ(mesh_.boundaryId(f), BoundaryId{0});

    } else {

      EXPECT_EQ(mesh_.boundaryId(f), invalid_boundary);
    }
  }
}

// ============================================================
// Dense ID invariants.
//
// All IDs returned by topology must satisfy:
//
// CellId ∈ [0, numCells)
//
// ============================================================

TEST_F(MoabMeshTest, UsesDenseValidCellIds) {
  for (FaceId f = 0; f < mesh_.numFaces(); ++f) {

    const CellId owner = mesh_.owner(f);

    EXPECT_LT(owner, mesh_.numCells());

    if (!mesh_.isBoundary(f)) {

      const CellId neighbor = mesh_.neighbor(f);

      EXPECT_LT(neighbor, mesh_.numCells());
    }
  }
}

// ============================================================
// Invalid index handling.
// ============================================================

TEST_F(MoabMeshTest, RejectsInvalidFaceId) {
  const FaceId invalid = static_cast<FaceId>(mesh_.numFaces());

  EXPECT_THROW(mesh_.owner(invalid), std::out_of_range);

  EXPECT_THROW(mesh_.neighbor(invalid), std::out_of_range);

  EXPECT_THROW(mesh_.isBoundary(invalid), std::out_of_range);

  EXPECT_THROW(mesh_.faceCenter(invalid), std::out_of_range);

  EXPECT_THROW(mesh_.faceArea(invalid), std::out_of_range);

  EXPECT_THROW(mesh_.faceNormal(invalid), std::out_of_range);

  EXPECT_THROW(mesh_.boundaryId(invalid), std::out_of_range);
}

TEST_F(MoabMeshTest, RejectsInvalidCellId) {
  const CellId invalid = static_cast<CellId>(mesh_.numCells());

  EXPECT_THROW(mesh_.cellFaces(invalid), std::out_of_range);

  EXPECT_THROW(mesh_.cellCenter(invalid), std::out_of_range);

  EXPECT_THROW(mesh_.cellVolume(invalid), std::out_of_range);
}

}  // namespace pemu::mesh::test