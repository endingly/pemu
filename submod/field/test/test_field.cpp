#include <gtest/gtest.h>

#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>
#include <pemu/mesh/moab_mesh.hpp>

#include <filesystem>
#include <stdexcept>

namespace pemu::field::test {

namespace {

std::filesystem::path testMeshPath() {
#ifdef PEMU_MESH_TEST_DATA_DIR
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
#else
  return std::filesystem::path{"two_quads.msh"};
#endif
}

class FieldTest : public ::testing::Test {
 protected:
  FieldTest() : mesh_(testMeshPath().string()) {}

  mesh::MoabMesh mesh_;
};

}  // namespace

// ------------------------------------------------------------
// CellField
// ------------------------------------------------------------

TEST_F(FieldTest, CellFieldSizeMatchesMesh) {
  CellField<double> field(mesh_);

  EXPECT_EQ(field.size(), mesh_.numCells());

  EXPECT_EQ(field.size(), 2u);
}

TEST_F(FieldTest, CellFieldSupportsInitialValue) {
  CellField<double> field(mesh_, 3.5);

  for (mesh::CellId c = 0; c < mesh_.numCells(); ++c) {

    EXPECT_DOUBLE_EQ(field[c], 3.5);
  }
}

TEST_F(FieldTest, CellFieldSupportsDenseIdAccess) {
  CellField<double> field(mesh_);

  field[0] = 10.0;
  field[1] = 20.0;

  EXPECT_DOUBLE_EQ(field[0], 10.0);

  EXPECT_DOUBLE_EQ(field[1], 20.0);
}

TEST_F(FieldTest, CellFieldFillWorks) {
  CellField<double> field(mesh_);

  field.fill(7.0);

  for (const auto value : field) {

    EXPECT_DOUBLE_EQ(value, 7.0);
  }
}

TEST_F(FieldTest, CellFieldAtRejectsInvalidCell) {
  CellField<double> field(mesh_);

  EXPECT_THROW(field.at(static_cast<mesh::CellId>(mesh_.numCells())),
               std::out_of_range);
}

TEST_F(FieldTest, CellFieldKeepsMeshAssociation) {
  CellField<double> field(mesh_);

  EXPECT_EQ(&field.mesh(), &mesh_);
}

// ------------------------------------------------------------
// FaceField
// ------------------------------------------------------------

TEST_F(FieldTest, FaceFieldSizeMatchesMesh) {
  FaceField<double> field(mesh_);

  EXPECT_EQ(field.size(), mesh_.numFaces());

  EXPECT_EQ(field.size(), 7u);
}

TEST_F(FieldTest, FaceFieldSupportsInitialValue) {
  FaceField<double> field(mesh_, 2.0);

  for (mesh::FaceId f = 0; f < mesh_.numFaces(); ++f) {

    EXPECT_DOUBLE_EQ(field[f], 2.0);
  }
}

TEST_F(FieldTest, FaceFieldSupportsDenseIdAccess) {
  FaceField<double> field(mesh_);

  for (mesh::FaceId f = 0; f < mesh_.numFaces(); ++f) {

    field[f] = static_cast<double>(f);
  }

  for (mesh::FaceId f = 0; f < mesh_.numFaces(); ++f) {

    EXPECT_DOUBLE_EQ(field[f], static_cast<double>(f));
  }
}

TEST_F(FieldTest, FaceFieldAtRejectsInvalidFace) {
  FaceField<double> field(mesh_);

  EXPECT_THROW(field.at(static_cast<mesh::FaceId>(mesh_.numFaces())),
               std::out_of_range);
}

TEST_F(FieldTest, FaceFieldKeepsMeshAssociation) {
  FaceField<double> field(mesh_);

  EXPECT_EQ(&field.mesh(), &mesh_);
}

}  // namespace pemu::field::test