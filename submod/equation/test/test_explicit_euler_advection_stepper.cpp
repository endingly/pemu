#include <gtest/gtest.h>
#include <filesystem>
#include <pemu/equation/fixed_step_explicit_euler_advection_stepper.hpp>
#include <pemu/mesh/moab_mesh.hpp>

namespace pemu::equation::test {

namespace {

std::filesystem::path testMeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
}

class FixedStepExplicitEulerAdvectionStepperTest : public ::testing::Test {
 protected:
  FixedStepExplicitEulerAdvectionStepperTest()
      : mesh_(testMeshPath().string()) {}

  mesh::MoabMesh mesh_;
};

mesh::FaceId findInternalFace(const mesh::IMesh& mesh) {
  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {

    if (!mesh.isBoundary(face)) {
      return face;
    }
  }
  throw std::runtime_error("internal face not found");
}

};  // namespace

TEST_F(FixedStepExplicitEulerAdvectionStepperTest, ComputesExpectedCfl) {
  field::FaceField<double> velocity(mesh_, 0.0);

  const auto face = findInternalFace(mesh_);

  velocity[face] = 2.0;

  boundary::BoundaryConditionSet bc;

  equation::FixedStepExplicitEulerAdvectionStepper stepper(mesh_, velocity,
                                                           0.25, bc);

  EXPECT_NEAR(stepper.maxCfl(), 0.5, 1e-12);
}

TEST_F(FixedStepExplicitEulerAdvectionStepperTest, RejectsCflGreaterThanOne) {
  field::CellField<double> u(mesh_, 1.0);

  field::FaceField<double> velocity(mesh_, 0.0);

  const auto face = findInternalFace(mesh_);

  velocity[face] = 2.0;

  boundary::BoundaryConditionSet bc;

  equation::FixedStepExplicitEulerAdvectionStepper stepper(mesh_, velocity,
                                                           0.75, bc);

  EXPECT_GT(stepper.maxCfl(), 1.0);

  EXPECT_THROW(stepper.step(u), std::runtime_error);
}

TEST_F(FixedStepExplicitEulerAdvectionStepperTest,
       CflOneMovesStateAcrossInternalFace) {
  field::CellField<double> u(mesh_, 0.0);

  field::FaceField<double> velocity(mesh_, 0.0);

  const auto face = findInternalFace(mesh_);

  const auto owner = mesh_.owner(face);

  const auto neighbor = mesh_.neighbor(face);

  u[owner] = 1.0;

  u[neighbor] = 0.0;

  velocity[face] = 1.0;

  boundary::BoundaryConditionSet bc;

  equation::FixedStepExplicitEulerAdvectionStepper stepper(mesh_, velocity, 1.0,
                                                           bc);

  ASSERT_NEAR(stepper.maxCfl(), 1.0, 1e-12);

  stepper.step(u);

  EXPECT_NEAR(u[owner], 0.0, 1e-12);

  EXPECT_NEAR(u[neighbor], 1.0, 1e-12);

  const double total = u[owner] * mesh_.cellVolume(owner) +
                       u[neighbor] * mesh_.cellVolume(neighbor);

  EXPECT_NEAR(total, 1.0, 1e-12);
}

};  // namespace pemu::equation::test
