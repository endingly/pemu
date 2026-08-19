#include <gtest/gtest.h>
#include <filesystem>
#include <pemu/boundary/boundary_condition_set.hpp>
#include <pemu/discretization/operators/electric_field.hpp>
#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>
#include <pemu/mesh/moab_mesh.hpp>

namespace pemu::discretization::test {

namespace {

mesh::FaceId findBoundaryFace(const mesh::IMesh& mesh, mesh::BoundaryId id) {
  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {

    if (!mesh.isBoundary(face)) {
      continue;
    }

    if (mesh.boundaryId(face) == id) {
      return face;
    }
  }

  throw std::runtime_error("boundary face not found");
}

std::filesystem::path testMeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
}

class ElectricFieldTest : public ::testing::Test {
 protected:
  ElectricFieldTest() : mesh_(testMeshPath().string()) {}

  mesh::MoabMesh mesh_;
};

};  // namespace

TEST_F(ElectricFieldTest, LinearPotentialProducesExactNormalElectricField) {
  constexpr double epsilon = 2.0;

  field::CellField<double> phi(mesh_, 0.0);

  field::FaceField<double> electric_field(mesh_, 0.0);

  // --------------------------------------------------------
  // phi(x,y) = x
  // --------------------------------------------------------

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    phi[cell] = mesh_.cellCenter(cell).x;
  }

  boundary::BoundaryConditionSet bc;

  // --------------------------------------------------------
  // Left / right:
  //
  // Dirichlet phi = x_boundary
  // --------------------------------------------------------

  const auto left_face = findBoundaryFace(mesh_, mesh::BoundaryId{1});

  const auto right_face = findBoundaryFace(mesh_, mesh::BoundaryId{2});

  bc.setDirichlet(mesh::BoundaryId{1}, mesh_.faceCenter(left_face).x);

  bc.setDirichlet(mesh::BoundaryId{2}, mesh_.faceCenter(right_face).x);

  // --------------------------------------------------------
  // Bottom / top:
  //
  // -epsilon grad(phi)·n = 0
  // --------------------------------------------------------

  bc.setNeumann(mesh::BoundaryId{3}, 0.0);

  bc.setNeumann(mesh::BoundaryId{4}, 0.0);

  discretization::operators::electricFieldNormal(phi, epsilon, bc,
                                                 electric_field);

  // --------------------------------------------------------
  // Exact:
  //
  //     E = (-1, 0)
  //
  // therefore:
  //
  //     E_n = -n_x
  // --------------------------------------------------------

  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {

    const auto normal = mesh_.faceNormal(face);

    EXPECT_NEAR(electric_field[face], -normal.x, 1e-12);
  }
}

};  // namespace pemu::discretization::test