#include <gtest/gtest.h>

#include <pemu/equation/adaptive_step_multi_species_drift_diffusion_stepper.hpp>
#include <pemu/equation/fixed_step_multi_species_drift_diffusion_stepper.hpp>
#include <pemu/linalg/cholmod_solver.hpp>
#include <pemu/mesh/moab_mesh.hpp>
#include <pemu/physics/wall/types.hpp>

#include <filesystem>
#include <memory>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace pemu::equation::test {
namespace {

static_assert(
    !std::is_move_constructible_v<FixedStepMultiSpeciesDriftDiffusionStepper>);
static_assert(
    !std::is_move_assignable_v<FixedStepMultiSpeciesDriftDiffusionStepper>);

/** @brief Returns the two-cell coupled-wall test mesh. */
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

/** @brief Creates grounded electrostatic boundaries. */
[[nodiscard]] boundary::BoundaryConditionSet potentialBoundaries() {
  boundary::BoundaryConditionSet conditions;
  for (mesh::BoundaryId boundary = 1; boundary <= 4; ++boundary) {
    conditions.setDirichlet(boundary, 0.0);
  }
  return conditions;
}

/** @brief Creates density-matching SG boundaries outside the left wall. */
[[nodiscard]] std::vector<boundary::BoundaryConditionSet> speciesBoundaries(
    std::size_t count, double density) {
  std::vector<boundary::BoundaryConditionSet> conditions(count);
  for (auto& species_conditions : conditions) {
    species_conditions.setDirichlet(mesh::BoundaryId{2}, density);
    species_conditions.setDirichlet(mesh::BoundaryId{3}, density);
    species_conditions.setDirichlet(mesh::BoundaryId{4}, density);
  }
  return conditions;
}

/** @brief Integrates one species density over physical cell volumes. */
[[nodiscard]] double integratedDensity(
    const mesh::IMesh& mesh, const field::CellField<double>& density) {
  double result = 0.0;
  for (mesh::CellId cell = 0; cell < mesh.numCells(); ++cell) {
    result += density[cell] * mesh.cellVolume(cell);
  }
  return result;
}

}  // namespace

TEST(AdaptiveMultiSpeciesWallTransportTest,
     CouplesIonLossToSecondaryElectronFluxAndConservesFaceBalance) {
  mesh::MoabMesh mesh(testMeshPath().string());
  physics::SpeciesSet species;
  const auto electron = species.add(
      {.name = "e",
       .charge = -1.0,
       .mobility = 1.0,
       .diffusivity = 0.1,
       .transport_model = physics::SpeciesTransportModel::DriftDiffusion});
  const auto ion = species.add(
      {.name = "ion",
       .charge = +1.0,
       .mobility = 0.5,
       .diffusivity = 0.1,
       .transport_model = physics::SpeciesTransportModel::DriftDiffusion});
  constexpr double initial_density = 7.0;
  physics::SpeciesCellFields density(mesh, species.size(), initial_density);
  physics::SpeciesCellFields source(mesh, species.size(), 0.0);
  physics::wall::WallBoundarySet wall_boundaries;
  wall_boundaries.set(mesh::BoundaryId{1},
                      {.losses = {{electron, 2.0, 5.0}, {ion, 3.0, 0.0}},
                       .secondary_emissions = {{ion, electron, 0.1, 4.0}}});
  AdaptiveStepMultiSpeciesDriftDiffusionStepper stepper(
      mesh, species, 1.0, potentialBoundaries(),
      speciesBoundaries(species.size(), initial_density),
      std::make_unique<linalg::CholmodSolver>(),
      {.safety = 1.0, .min_dt = 1e-8, .max_dt = 1e-3, .max_growth = 2.0}, {},
      std::nullopt, std::move(wall_boundaries));
  const auto wall_face = findBoundaryFace(mesh, mesh::BoundaryId{1});

  ASSERT_TRUE(stepper.prepareElectrostatics(density).success());
  field::FaceField<double> electron_flux(mesh, 0.0);
  field::FaceField<double> ion_flux(mesh, 0.0);
  stepper.computeParticleFluxNormal(density, electron, electron_flux);
  stepper.computeParticleFluxNormal(density, ion, ion_flux);
  EXPECT_DOUBLE_EQ(ion_flux[wall_face], 21.0);
  EXPECT_DOUBLE_EQ(electron_flux[wall_face], 11.9);
  EXPECT_DOUBLE_EQ(
      stepper.wallFluxAssembler().energyInwardFlux()[electron][wall_face], 8.4);

  const double electron_before = integratedDensity(mesh, density[electron]);
  const double ion_before = integratedDensity(mesh, density[ion]);
  const auto proposal = stepper.proposeTimeStep(density, source, 1e-3);
  stepper.advancePrepared(density, source, proposal);
  const double boundary_area = mesh.faceArea(wall_face);

  EXPECT_NEAR(integratedDensity(mesh, density[electron]) - electron_before,
              -proposal.dt * 11.9 * boundary_area, 1e-12);
  EXPECT_NEAR(integratedDensity(mesh, density[ion]) - ion_before,
              -proposal.dt * 21.0 * boundary_area, 1e-12);
}

TEST(FixedMultiSpeciesWallTransportTest,
     UsesTheSameSecondaryEmissionFluxContract) {
  mesh::MoabMesh mesh(testMeshPath().string());
  physics::SpeciesSet species;
  const auto electron = species.add(
      {.name = "e",
       .charge = -1.0,
       .mobility = 1.0,
       .diffusivity = 0.1,
       .transport_model = physics::SpeciesTransportModel::DriftDiffusion});
  const auto ion = species.add(
      {.name = "ion",
       .charge = +1.0,
       .mobility = 0.5,
       .diffusivity = 0.1,
       .transport_model = physics::SpeciesTransportModel::DriftDiffusion});
  constexpr double initial_density = 7.0;
  physics::SpeciesCellFields density(mesh, species.size(), initial_density);
  physics::SpeciesCellFields source(mesh, species.size(), 0.0);
  physics::wall::WallBoundarySet wall_boundaries;
  wall_boundaries.set(mesh::BoundaryId{1},
                      {.losses = {{electron, 2.0, 5.0}, {ion, 3.0, 0.0}},
                       .secondary_emissions = {{ion, electron, 0.1, 4.0}}});
  FixedStepMultiSpeciesDriftDiffusionStepper stepper(
      mesh, species, 1.0, 1e-3, potentialBoundaries(),
      speciesBoundaries(species.size(), initial_density),
      std::make_unique<linalg::CholmodSolver>(), {}, std::nullopt,
      std::move(wall_boundaries));
  const auto wall_face = findBoundaryFace(mesh, mesh::BoundaryId{1});

  ASSERT_TRUE(stepper.updateElectrostatics(density).success());
  field::FaceField<double> electron_flux(mesh, 0.0);
  stepper.computeParticleFluxNormal(density, electron, electron_flux);
  EXPECT_DOUBLE_EQ(electron_flux[wall_face], 11.9);
  EXPECT_DOUBLE_EQ(
      stepper.wallFluxAssembler().energyInwardFlux()[electron][wall_face], 8.4);
  const double before = integratedDensity(mesh, density[electron]);

  stepper.advanceTransport(density, source);

  EXPECT_NEAR(integratedDensity(mesh, density[electron]) - before,
              -1e-3 * 11.9 * mesh.faceArea(wall_face), 1e-12);
}

}  // namespace pemu::equation::test
