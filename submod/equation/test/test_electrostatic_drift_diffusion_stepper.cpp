#include <gtest/gtest.h>

#include <pemu/boundary/boundary_condition_set.hpp>

#include <pemu/equation/electrostatic_drift_diffusion_stepper.hpp>

#include <pemu/field/cell_field.hpp>

#include <pemu/linalg/cholmod_solver.hpp>

#include <pemu/mesh/moab_mesh.hpp>

#include <pemu/physics/charged_species_transport.hpp>

#include <cmath>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <vector>

namespace {

using namespace pemu;

// ============================================================
// Test data path
// ============================================================

std::filesystem::path testMeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
}

// ============================================================
// Find one face belonging to a physical boundary.
// ============================================================

mesh::FaceId findBoundaryFace(const mesh::IMesh& mesh,
                              mesh::BoundaryId boundary_id) {
  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {

    if (!mesh.isBoundary(face)) {
      continue;
    }

    if (mesh.boundaryId(face) == boundary_id) {

      return face;
    }
  }

  throw std::runtime_error("boundary face not found");
}

// ============================================================
// Test fixture
// ============================================================

class ElectrostaticDriftDiffusionTest : public ::testing::Test {
 protected:
  ElectrostaticDriftDiffusionTest() : mesh_(testMeshPath().string()) {}

  // --------------------------------------------------------
  // phi = 0 on every boundary.
  // --------------------------------------------------------

  boundary::BoundaryConditionSet makeZeroPotentialBoundaryConditions() {
    boundary::BoundaryConditionSet bc;

    bc.setDirichlet(mesh::BoundaryId{1}, 0.0);

    bc.setDirichlet(mesh::BoundaryId{2}, 0.0);

    bc.setDirichlet(mesh::BoundaryId{3}, 0.0);

    bc.setDirichlet(mesh::BoundaryId{4}, 0.0);

    return bc;
  }

  // --------------------------------------------------------
  // Constant density Dirichlet BC.
  // --------------------------------------------------------

  boundary::BoundaryConditionSet makeConstantSpeciesBoundaryConditions(
      double value) {
    boundary::BoundaryConditionSet bc;

    bc.setDirichlet(mesh::BoundaryId{1}, value);

    bc.setDirichlet(mesh::BoundaryId{2}, value);

    bc.setDirichlet(mesh::BoundaryId{3}, value);

    bc.setDirichlet(mesh::BoundaryId{4}, value);

    return bc;
  }

  // --------------------------------------------------------
  // Potential:
  //
  //     phi(x,y) = x
  //
  // Left/right use Dirichlet.
  // Top/bottom use zero normal derivative.
  // --------------------------------------------------------

  boundary::BoundaryConditionSet makeLinearXPotentialBoundaryConditions() {
    boundary::BoundaryConditionSet bc;

    const auto left_face = findBoundaryFace(mesh_, mesh::BoundaryId{1});

    const auto right_face = findBoundaryFace(mesh_, mesh::BoundaryId{2});

    bc.setDirichlet(mesh::BoundaryId{1}, mesh_.faceCenter(left_face).x);

    bc.setDirichlet(mesh::BoundaryId{2}, mesh_.faceCenter(right_face).x);

    bc.setNeumann(mesh::BoundaryId{3}, 0.0);

    bc.setNeumann(mesh::BoundaryId{4}, 0.0);

    return bc;
  }

 protected:
  mesh::MoabMesh mesh_;
};

// ============================================================
// 1. Uniform neutral plasma
//
//     n_e = n_i
//
// therefore:
//
//     rho = 0
//     phi = 0
//     E = 0
//     v_e = v_i = 0
//
// With constant species density and zero sources,
// the complete coupled timestep must leave everything
// stationary.
// ============================================================

TEST_F(ElectrostaticDriftDiffusionTest, UniformNeutralPlasmaRemainsStationary) {
  constexpr double density_value = 3.0;

  constexpr double epsilon = 1.0;

  constexpr double dt = 0.01;

  field::CellField<double> electron_density(mesh_, density_value);

  field::CellField<double> ion_density(mesh_, density_value);

  field::CellField<double> electron_source(mesh_, 0.0);

  field::CellField<double> ion_source(mesh_, 0.0);

  auto potential_bc = makeZeroPotentialBoundaryConditions();

  auto electron_bc = makeConstantSpeciesBoundaryConditions(density_value);

  auto ion_bc = makeConstantSpeciesBoundaryConditions(density_value);

  physics::ChargedSpeciesTransport electron{
      .charge = -1.0, .mobility = 1.0, .diffusivity = 0.1};

  physics::ChargedSpeciesTransport ion{
      .charge = +1.0, .mobility = 0.5, .diffusivity = 0.1};

  equation::ElectrostaticDriftDiffusionStepper stepper(
      mesh_, epsilon, dt, electron, ion, potential_bc, electron_bc, ion_bc,
      std::make_unique<linalg::CholmodSolver>());

  const auto result =
      stepper.step(electron_density, electron_source, ion_density, ion_source);

  ASSERT_TRUE(result.success());

  // --------------------------------------------------------
  // Cell quantities
  // --------------------------------------------------------

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(stepper.chargeDensity()[cell], 0.0, 1e-12);

    EXPECT_NEAR(stepper.potential()[cell], 0.0, 1e-12);

    EXPECT_NEAR(electron_density[cell], density_value, 1e-12);

    EXPECT_NEAR(ion_density[cell], density_value, 1e-12);
  }

  // --------------------------------------------------------
  // Face quantities
  // --------------------------------------------------------

  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {

    EXPECT_NEAR(stepper.electricFieldNormal()[face], 0.0, 1e-12);

    EXPECT_NEAR(stepper.electronDriftVelocityNormal()[face], 0.0, 1e-12);

    EXPECT_NEAR(stepper.ionDriftVelocityNormal()[face], 0.0, 1e-12);
  }
}

// ============================================================
// 2. Applied linear potential
//
//     phi = x
//
// therefore:
//
//     E = (-1, 0)
//
// electron:
//     v_e = -mu_e E = (+mu_e, 0)
//
// positive ion:
//     v_i = +mu_i E = (-mu_i, 0)
//
// This test only updates electrostatics.
// ============================================================

TEST_F(ElectrostaticDriftDiffusionTest,
       AppliedPotentialProducesOppositeSpeciesDrift) {
  constexpr double density_value = 1.0;

  constexpr double epsilon = 1.0;

  constexpr double dt = 0.01;

  field::CellField<double> electron_density(mesh_, density_value);

  field::CellField<double> ion_density(mesh_, density_value);

  auto potential_bc = makeLinearXPotentialBoundaryConditions();

  auto electron_bc = makeConstantSpeciesBoundaryConditions(density_value);

  auto ion_bc = makeConstantSpeciesBoundaryConditions(density_value);

  physics::ChargedSpeciesTransport electron{
      .charge = -1.0, .mobility = 2.0, .diffusivity = 0.1};

  physics::ChargedSpeciesTransport ion{
      .charge = +1.0, .mobility = 1.0, .diffusivity = 0.1};

  equation::ElectrostaticDriftDiffusionStepper stepper(
      mesh_, epsilon, dt, electron, ion, potential_bc, electron_bc, ion_bc,
      std::make_unique<linalg::CholmodSolver>());

  const auto result =
      stepper.updateElectrostatics(electron_density, ion_density);

  ASSERT_TRUE(result.success());

  // Neutral density:
  //
  //     rho = -n_e + n_i = 0
  //
  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(stepper.chargeDensity()[cell], 0.0, 1e-12);
  }

  // phi=x gives:
  //
  //     E_n = (-1,0) dot n
  //         = -n_x
  //
  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {

    const auto normal = mesh_.faceNormal(face);

    EXPECT_NEAR(stepper.electricFieldNormal()[face], -normal.x, 1e-12);

    // electron:
    //
    // v_e = -2 E
    //
    //     = (+2, 0)
    //
    EXPECT_NEAR(stepper.electronDriftVelocityNormal()[face], +2.0 * normal.x,
                1e-12);

    // ion:
    //
    // v_i = +1 E
    //
    //     = (-1, 0)
    //
    EXPECT_NEAR(stepper.ionDriftVelocityNormal()[face], -1.0 * normal.x, 1e-12);
  }
}

// ============================================================
// 3. CFL rejection must not partially advance species.
//
// We deliberately choose:
//
//     phi = x
//
// so E != 0, and choose a very large dt.
//
// Both density fields are snapshotted before step().
//
// If transport CFL rejects the timestep, neither density
// field may change.
// ============================================================

TEST_F(ElectrostaticDriftDiffusionTest,
       CflFailureDoesNotPartiallyAdvanceSpecies) {
  constexpr double epsilon = 1.0;

  //
  // two_quads consists of unit cells.
  //
  // This is deliberately enormous for explicit transport.
  //
  constexpr double dt = 100.0;

  // --------------------------------------------------------
  // Use a nonuniform state as well, so accidental updates
  // are easy to detect.
  // --------------------------------------------------------

  field::CellField<double> electron_density(mesh_, 0.0);

  field::CellField<double> ion_density(mesh_, 0.0);

  ASSERT_EQ(mesh_.numCells(), 2u);

  electron_density[0] = 1.0;

  electron_density[1] = 2.0;

  ion_density[0] = 3.0;

  ion_density[1] = 4.0;

  field::CellField<double> electron_source(mesh_, 0.0);

  field::CellField<double> ion_source(mesh_, 0.0);

  // --------------------------------------------------------
  // Actual snapshots.
  //
  // Do NOT use span here because span is non-owning.
  // --------------------------------------------------------

  const std::vector<double> electron_snapshot(electron_density.begin(),
                                              electron_density.end());

  const std::vector<double> ion_snapshot(ion_density.begin(),
                                         ion_density.end());

  // --------------------------------------------------------
  // Applied electric field:
  //
  // phi = x
  // --------------------------------------------------------

  auto potential_bc = makeLinearXPotentialBoundaryConditions();

  //
  // SG currently requires Dirichlet species BC.
  //
  auto electron_bc = makeConstantSpeciesBoundaryConditions(1.0);

  auto ion_bc = makeConstantSpeciesBoundaryConditions(1.0);

  physics::ChargedSpeciesTransport electron{
      .charge = -1.0, .mobility = 2.0, .diffusivity = 0.1};

  physics::ChargedSpeciesTransport ion{
      .charge = +1.0, .mobility = 1.0, .diffusivity = 0.1};

  // ========================================================
  // THIS was the missing object in the previous answer.
  // ========================================================

  equation::ElectrostaticDriftDiffusionStepper stepper(
      mesh_, epsilon, dt, electron, ion, potential_bc, electron_bc, ion_bc,
      std::make_unique<linalg::CholmodSolver>());

  // --------------------------------------------------------
  // First explicitly update electrostatics so that the CFL
  // diagnostic is meaningful and testable before step().
  // --------------------------------------------------------

  const auto electrostatic_result =
      stepper.updateElectrostatics(electron_density, ion_density);

  ASSERT_TRUE(electrostatic_result.success());

  // At least one species must violate the explicit
  // transport timestep restriction.
  EXPECT_TRUE(stepper.electronTransportCfl() > 1.0 ||
              stepper.ionTransportCfl() > 1.0);

  // --------------------------------------------------------
  // Coupled step must reject before either species advances.
  // --------------------------------------------------------

  EXPECT_THROW(auto r = stepper.step(electron_density, electron_source,
                                     ion_density, ion_source),
               std::runtime_error);

  // --------------------------------------------------------
  // No half-step modification is allowed.
  // --------------------------------------------------------

  ASSERT_EQ(electron_snapshot.size(), electron_density.size());

  ASSERT_EQ(ion_snapshot.size(), ion_density.size());

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_DOUBLE_EQ(electron_density[cell], electron_snapshot[cell]);

    EXPECT_DOUBLE_EQ(ion_density[cell], ion_snapshot[cell]);
  }
}

// ============================================================
// 4. Invalid species sign should be rejected.
//
// This locks the semantic convention:
//
//     electron charge < 0
//     ion charge      > 0
// ============================================================

TEST_F(ElectrostaticDriftDiffusionTest, RejectsInvalidSpeciesChargePolarity) {
  constexpr double epsilon = 1.0;

  constexpr double dt = 0.01;

  auto potential_bc = makeZeroPotentialBoundaryConditions();

  auto electron_bc = makeConstantSpeciesBoundaryConditions(1.0);

  auto ion_bc = makeConstantSpeciesBoundaryConditions(1.0);

  physics::ChargedSpeciesTransport invalid_electron{
      .charge = +1.0, .mobility = 1.0, .diffusivity = 0.1};

  physics::ChargedSpeciesTransport ion{
      .charge = +1.0, .mobility = 1.0, .diffusivity = 0.1};

  EXPECT_THROW(
      equation::ElectrostaticDriftDiffusionStepper(
          mesh_, epsilon, dt, invalid_electron, ion, potential_bc, electron_bc,
          ion_bc, std::make_unique<linalg::CholmodSolver>()),
      std::invalid_argument);
}

}  // namespace