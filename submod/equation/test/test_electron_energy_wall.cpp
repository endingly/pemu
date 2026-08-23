#include <gtest/gtest.h>

#include <pemu/equation/explicit_electron_energy_stepper.hpp>
#include <pemu/mesh/moab_mesh.hpp>

#include <filesystem>
#include <stdexcept>
#include <utility>

namespace pemu::equation::test {
namespace {

/** @brief Returns the two-cell energy-wall test mesh. */
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

/** @brief Integrates energy density over physical cell volumes. */
[[nodiscard]] double integratedEnergy(
    const mesh::IMesh& mesh, const field::CellField<double>& energy_density) {
  double result = 0.0;
  for (mesh::CellId cell = 0; cell < mesh.numCells(); ++cell) {
    result += energy_density[cell] * mesh.cellVolume(cell);
  }
  return result;
}

}  // namespace

TEST(ElectronEnergyWallTest,
     AppliesIndependentEnergyLossAndSecondaryEmissionFlux) {
  mesh::MoabMesh mesh(testMeshPath().string());
  const auto wall_face = findBoundaryFace(mesh, mesh::BoundaryId{1});
  field::FaceField<double> electron_velocity(mesh, 0.0);
  field::FaceField<std::uint8_t> active(mesh, std::uint8_t{0});
  field::FaceField<double> energy_loss_velocity(mesh, 0.0);
  field::FaceField<double> energy_inward_flux(mesh, 0.0);
  active[wall_face] = std::uint8_t{1};
  energy_loss_velocity[wall_face] = 5.0;
  energy_inward_flux[wall_face] = 8.4;
  const discretization::operators::LinearBoundaryFluxView wall_flux(
      active, energy_loss_velocity, energy_inward_flux);
  boundary::BoundaryConditionSet conditions;
  conditions.setDirichlet(mesh::BoundaryId{2}, 11.0);
  conditions.setDirichlet(mesh::BoundaryId{3}, 11.0);
  conditions.setDirichlet(mesh::BoundaryId{4}, 11.0);
  ExplicitElectronEnergyStepper stepper(mesh, electron_velocity, 1.0,
                                        std::move(conditions), wall_flux);
  field::CellField<double> energy_density(mesh, 11.0);
  field::CellField<double> source(mesh, 0.0);
  field::CellField<double> increment(mesh, 0.0);
  constexpr double dt = 0.01;
  const double before = integratedEnergy(mesh, energy_density);

  stepper.computeIncrement(energy_density, source, dt, increment);
  for (mesh::CellId cell = 0; cell < mesh.numCells(); ++cell) {
    energy_density[cell] += increment[cell];
  }

  const double expected_flux = 5.0 * 11.0 - 8.4;
  EXPECT_NEAR(integratedEnergy(mesh, energy_density) - before,
              -dt * expected_flux * mesh.faceArea(wall_face), 1e-12);
}

TEST(ElectronEnergyWallTest, WallLossContributesToPositivityLimit) {
  mesh::MoabMesh mesh(testMeshPath().string());
  const auto wall_face = findBoundaryFace(mesh, mesh::BoundaryId{1});
  field::FaceField<double> electron_velocity(mesh, 0.0);
  field::FaceField<std::uint8_t> active(mesh, std::uint8_t{0});
  field::FaceField<double> energy_loss_velocity(mesh, 0.0);
  field::FaceField<double> energy_inward_flux(mesh, 0.0);
  active[wall_face] = std::uint8_t{1};
  energy_loss_velocity[wall_face] = 1000.0;
  const discretization::operators::LinearBoundaryFluxView wall_flux(
      active, energy_loss_velocity, energy_inward_flux);
  boundary::BoundaryConditionSet conditions;
  conditions.setDirichlet(mesh::BoundaryId{2}, 1.0);
  conditions.setDirichlet(mesh::BoundaryId{3}, 1.0);
  conditions.setDirichlet(mesh::BoundaryId{4}, 1.0);
  ExplicitElectronEnergyStepper stepper(mesh, electron_velocity, 1.0,
                                        std::move(conditions), wall_flux);
  field::CellField<double> energy_density(mesh, 1.0);
  field::CellField<double> source(mesh, 0.0);

  EXPECT_LT(stepper.maxPositiveTimeStep(energy_density, source), 0.01);
}

}  // namespace pemu::equation::test
