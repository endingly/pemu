#include <gtest/gtest.h>

#include <pemu/discretization/operators/cell_vector_reconstruction.hpp>
#include <pemu/discretization/operators/divergence.hpp>
#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>
#include <pemu/mesh/axisymmetric_mesh_view.hpp>
#include <pemu/mesh/moab_mesh.hpp>

#include <filesystem>

namespace pemu::discretization::test {

namespace {

// Returns the shared radial two-cell mesh independently of CTest's directory.
[[nodiscard]] std::filesystem::path radialMeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
}

}  // namespace

TEST(AxisymmetricDivergenceTest, RecoversRadialLinearFluxDivergence) {
  mesh::MoabMesh planar_mesh(radialMeshPath().string());
  mesh::AxisymmetricMeshView mesh(planar_mesh);
  field::FaceField<double> radial_flux(mesh, 0.0);
  field::CellField<double> divergence(mesh, 0.0);

  // For Gamma = r e_r, Gamma.n = r n_r and
  // div_axi(Gamma) = (1/r) d(r^2)/dr = 2.
  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {
    radial_flux[face] = mesh.faceCenter(face).x * mesh.faceNormal(face).x;
  }

  operators::divergence(radial_flux, divergence);

  for (mesh::CellId cell = 0; cell < mesh.numCells(); ++cell) {
    EXPECT_NEAR(divergence[cell], 2.0, 1e-12);
  }
}

TEST(AxisymmetricDivergenceTest,
     ReconstructsVectorMagnitudeAcrossZeroAreaAxisFace) {
  mesh::MoabMesh planar_mesh(radialMeshPath().string());
  mesh::AxisymmetricMeshView mesh(planar_mesh);
  field::FaceField<double> normal_component(mesh, 0.0);
  field::CellField<double> magnitude(mesh, 0.0);

  // A constant meridional vector (1, 2) remains reconstructible in the axis
  // cell because its positive-area internal and axial faces span both axes.
  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {
    const auto normal = mesh.faceNormal(face);
    normal_component[face] = normal.x + 2.0 * normal.y;
  }

  operators::reconstructCellVectorMagnitudeFromFaceNormal(normal_component,
                                                          magnitude);

  for (const double value : magnitude) {
    EXPECT_NEAR(value, std::sqrt(5.0), 1e-12);
  }
}

}  // namespace pemu::discretization::test
