#include <gtest/gtest.h>

#include <pemu/discretization/operators/cell_vector_reconstruction.hpp>
#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>
#include <pemu/mesh/moab_mesh.hpp>

#include <filesystem>
#include <limits>

namespace pemu::discretization::operators::test {
namespace {

/** @brief Returns the two-cell mesh used by vector-reconstruction tests. */
[[nodiscard]] std::filesystem::path testMeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
}

class CellVectorReconstructionTest : public ::testing::Test {
 protected:
  /** @brief Loads an orthogonal two-dimensional mesh. */
  CellVectorReconstructionTest() : mesh_(testMeshPath().string()) {}

  mesh::MoabMesh mesh_;
};

}  // namespace

TEST_F(CellVectorReconstructionTest, ReconstructsConstantVectorMagnitude) {
  field::FaceField<double> face_normal_component(mesh_, 0.0);
  field::CellField<double> cell_magnitude(mesh_, 0.0);
  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {
    const auto normal = mesh_.faceNormal(face);
    face_normal_component[face] = 3.0 * normal.x + 4.0 * normal.y;
  }

  reconstructCellVectorMagnitudeFromFaceNormal(face_normal_component,
                                               cell_magnitude);

  for (const double magnitude : cell_magnitude) {
    EXPECT_NEAR(magnitude, 5.0, 1e-12);
  }
}

TEST_F(CellVectorReconstructionTest, RejectsNonFiniteFaceDataAtomically) {
  field::FaceField<double> face_normal_component(mesh_, 0.0);
  field::CellField<double> cell_magnitude(mesh_, 12.0);
  face_normal_component[mesh::FaceId{0}] =
      std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(reconstructCellVectorMagnitudeFromFaceNormal(
                   face_normal_component, cell_magnitude),
               std::invalid_argument);
  for (const double magnitude : cell_magnitude) {
    EXPECT_DOUBLE_EQ(magnitude, 12.0);
  }
}

}  // namespace pemu::discretization::operators::test
