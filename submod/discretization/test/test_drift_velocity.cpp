#include <gtest/gtest.h>
#include <filesystem>
#include <pemu/boundary/boundary_condition_set.hpp>
#include <pemu/discretization/operators/drift_velocity.hpp>
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

class ElectrostaticDriftTest : public ::testing::Test {
 protected:
  ElectrostaticDriftTest() : mesh_(testMeshPath().string()) {}

  mesh::MoabMesh mesh_;
};

class DriftVelocityTest : public ::testing::Test {
 protected:
  DriftVelocityTest() : mesh_(testMeshPath().string()) {}

  mesh::MoabMesh mesh_;
};

};  // namespace

TEST_F(DriftVelocityTest, PositiveSpeciesDriftsAlongElectricField) {
  field::FaceField<double> electric_field(mesh_, -2.0);

  field::FaceField<double> velocity(mesh_, 0.0);

  discretization::operators::driftVelocityNormal(
      electric_field, 3.0, physics::ChargePolarity::Positive, velocity);

  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {

    EXPECT_NEAR(velocity[face], -6.0, 1e-12);
  }
}

TEST_F(DriftVelocityTest, NegativeSpeciesDriftsAgainstElectricField) {
  field::FaceField<double> electric_field(mesh_, -2.0);

  field::FaceField<double> velocity(mesh_, 0.0);

  discretization::operators::driftVelocityNormal(
      electric_field, 3.0, physics::ChargePolarity::Negative, velocity);

  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {

    EXPECT_NEAR(velocity[face], 6.0, 1e-12);
  }
}

TEST_F(ElectrostaticDriftTest, PotentialProducesOppositeIonAndElectronDrift) {
  constexpr double epsilon = 1.0;

  constexpr double mobility = 2.0;

  field::CellField<double> phi(mesh_, 0.0);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    phi[cell] = mesh_.cellCenter(cell).x;
  }

  boundary::BoundaryConditionSet bc;

  const auto left = findBoundaryFace(mesh_, mesh::BoundaryId{1});

  const auto right = findBoundaryFace(mesh_, mesh::BoundaryId{2});

  bc.setDirichlet(mesh::BoundaryId{1}, mesh_.faceCenter(left).x);

  bc.setDirichlet(mesh::BoundaryId{2}, mesh_.faceCenter(right).x);

  bc.setNeumann(mesh::BoundaryId{3}, 0.0);

  bc.setNeumann(mesh::BoundaryId{4}, 0.0);

  field::FaceField<double> electric_field(mesh_, 0.0);

  field::FaceField<double> ion_velocity(mesh_, 0.0);

  field::FaceField<double> electron_velocity(mesh_, 0.0);

  discretization::operators::electricFieldNormal(phi, epsilon, bc,
                                                 electric_field);

  discretization::operators::driftVelocityNormal(
      electric_field, mobility, physics::ChargePolarity::Positive,
      ion_velocity);

  discretization::operators::driftVelocityNormal(
      electric_field, mobility, physics::ChargePolarity::Negative,
      electron_velocity);

  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {

    const auto normal = mesh_.faceNormal(face);

    EXPECT_NEAR(electric_field[face], -normal.x, 1e-12);

    EXPECT_NEAR(ion_velocity[face], -2.0 * normal.x, 1e-12);

    EXPECT_NEAR(electron_velocity[face], +2.0 * normal.x, 1e-12);
  }
}

};  // namespace pemu::discretization::test