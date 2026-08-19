#include <gtest/gtest.h>
#include <filesystem>
#include <pemu/boundary/boundary_condition_set.hpp>
#include <pemu/discretization/operators/diffusion_flux.hpp>
#include <pemu/discretization/operators/divergence.hpp>
#include <pemu/mesh/moab_mesh.hpp>

#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>

namespace pemu::discretization::test {

namespace {
mesh::FaceId findInternalFace(const mesh::IMesh& mesh) {
  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {

    if (!mesh.isBoundary(face)) {
      return face;
    }
  }
  throw std::runtime_error("internal face not found");
}

std::filesystem::path testMeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
}

class OperatorTest : public ::testing::Test {
 protected:
  OperatorTest() : mesh_(testMeshPath().string()) {}

  mesh::MoabMesh mesh_;
};

};  // namespace

TEST_F(OperatorTest, DivergenceConservesInternalFlux) {
  field::FaceField<double> flux(mesh_, 0.0);

  field::CellField<double> div(mesh_, 0.0);

  const auto internal_face = findInternalFace(mesh_);

  flux[internal_face] = 2.5;

  discretization::operators::divergence(flux, div);

  double integral = 0.0;

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    integral += div[cell] * mesh_.cellVolume(cell);
  }

  EXPECT_NEAR(integral, 0.0, 1e-12);
}

TEST_F(OperatorTest, DivergenceUsesOwnerNeighborOrientation) {
  field::FaceField<double> flux(mesh_, 0.0);

  field::CellField<double> div(mesh_, 0.0);

  const auto face = findInternalFace(mesh_);

  const auto owner = mesh_.owner(face);

  const auto neighbor = mesh_.neighbor(face);

  flux[face] = 2.0;

  discretization::operators::divergence(flux, div);

  EXPECT_NEAR(div[owner], +2.0, 1e-12);

  EXPECT_NEAR(div[neighbor], -2.0, 1e-12);
}

TEST_F(OperatorTest, DivergenceSatisfiesDiscreteGaussTheorem) {
  field::FaceField<double> flux(mesh_, 0.0);

  field::CellField<double> div(mesh_, 0.0);

  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {

    if (mesh_.isBoundary(face)) {
      flux[face] = 1.0;
    }
  }

  discretization::operators::divergence(flux, div);

  double volume_integral = 0.0;

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    volume_integral += div[cell] * mesh_.cellVolume(cell);
  }

  double boundary_integral = 0.0;

  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {

    if (!mesh_.isBoundary(face)) {
      continue;
    }

    boundary_integral += flux[face] * mesh_.faceArea(face);
  }

  EXPECT_NEAR(volume_integral, boundary_integral, 1e-12);

  EXPECT_NEAR(volume_integral, 6.0, 1e-12);
}

TEST_F(OperatorTest, DiffusionFluxIsExactForLinearField) {
  field::CellField<double> u(mesh_);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    const auto center = mesh_.cellCenter(cell);

    u[cell] = center.x;
  }

  boundary::BoundaryConditionSet bc;

  // u = x

  bc.setDirichlet(1, 0.0);  // left
  bc.setDirichlet(2, 2.0);  // right

  //
  // y-boundaries:
  //
  // grad(u) dot n = 0
  //
  bc.setNeumann(3, 0.0);
  bc.setNeumann(4, 0.0);

  field::FaceField<double> flux(mesh_, 0.0);

  discretization::operators::diffusionFlux(u, 2.0, bc, flux);

  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {

    const auto n = mesh_.faceNormal(face);

    const double expected = -2.0 * n.x;

    EXPECT_NEAR(flux[face], expected, 1e-12);
  }
}

TEST_F(OperatorTest, DivergenceOfLinearDiffusionFluxIsZero) {
  field::CellField<double> u(mesh_);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    const auto center = mesh_.cellCenter(cell);

    u[cell] = center.x;
  }

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(1, 0.0);
  bc.setDirichlet(2, 2.0);

  bc.setNeumann(3, 0.0);
  bc.setNeumann(4, 0.0);

  field::FaceField<double> flux(mesh_, 0.0);

  field::CellField<double> div(mesh_, 0.0);

  discretization::operators::diffusionFlux(u, 1.0, bc, flux);

  discretization::operators::divergence(flux, div);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(div[cell], 0.0, 1e-12);
  }
}

TEST_F(OperatorTest, DivergenceRejectsDifferentMeshes) {
  mesh::MoabMesh another_mesh(testMeshPath().string());

  field::FaceField<double> flux(mesh_, 0.0);

  field::CellField<double> div(another_mesh, 0.0);

  EXPECT_THROW(discretization::operators::divergence(flux, div),
               std::invalid_argument);
}

};  // namespace pemu::discretization::test
