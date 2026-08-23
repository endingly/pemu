#include <gtest/gtest.h>
#include <pemu/boundary/boundary_condition_set.hpp>
#include <pemu/equation/fixed_step_electrostatic_drift_diffusion_stepper.hpp>
#include <pemu/field/cell_field.hpp>
#include <pemu/linalg/cholmod_solver.hpp>
#include <pemu/mesh/moab_mesh.hpp>
#include <pemu/physics/charged_species_transport.hpp>
#include <pemu/physics/reaction/mass_action.hpp>
#include <pemu/physics/reaction/network.hpp>

#include <cmath>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <vector>

namespace pemu::equation::test {

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

class FixedStepElectrostaticDriftDiffusionStepperTest : public ::testing::Test {
 protected:
  FixedStepElectrostaticDriftDiffusionStepperTest()
      : mesh_(testMeshPath().string()) {}

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

class FixedStepElectrostaticDriftDiffusionReactionTest
    : public ::testing::Test {
 protected:
  FixedStepElectrostaticDriftDiffusionReactionTest()
      : mesh_(testMeshPath().string()) {}
  static constexpr double kTolerance = 1e-12;
  mesh::MoabMesh mesh_;
};

// ============================================================
// Integral of number density:
//
//     N = sum_P n_P V_P
//
// This is not necessarily the literal integer particle count;
// it is the finite-volume integral of number density.
// ============================================================

double totalParticles(const mesh::IMesh& mesh,
                      const field::CellField<double>& density) {
  double total = 0.0;
  for (mesh::CellId cell = 0; cell < mesh.numCells(); ++cell) {
    total += density[cell] * mesh.cellVolume(cell);
  }
  return total;
}

};  // namespace

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

TEST_F(FixedStepElectrostaticDriftDiffusionStepperTest,
       UniformNeutralPlasmaRemainsStationary) {
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

  equation::FixedStepElectrostaticDriftDiffusionStepper stepper(
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

TEST_F(FixedStepElectrostaticDriftDiffusionStepperTest,
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

  equation::FixedStepElectrostaticDriftDiffusionStepper stepper(
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

TEST_F(FixedStepElectrostaticDriftDiffusionStepperTest,
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

  equation::FixedStepElectrostaticDriftDiffusionStepper stepper(
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

TEST_F(FixedStepElectrostaticDriftDiffusionStepperTest,
       RejectsInvalidSpeciesChargePolarity) {
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
      equation::FixedStepElectrostaticDriftDiffusionStepper(
          mesh_, epsilon, dt, invalid_electron, ion, potential_bc, electron_bc,
          ion_bc, std::make_unique<linalg::CholmodSolver>()),
      std::invalid_argument);
}

// ============================================================
// Electron-impact ionization:
//
//     e + N -> 2e + N+
//
// Net production:
//
//     S_e = R
//     S_i = R
//
// with:
//
//     R = k_ion n_e n_N
//
// Test configuration:
//
//     n_e^0 = n_i^0 = 1
//
//     n_N   = 4
//     k_ion = 0.5
//
// therefore:
//
//     R = 0.5 * 1 * 4 = 2
//
// dt = 0.1
//
// therefore:
//
//     n_e^1 = 1 + 0.1 * 2 = 1.2
//     n_i^1 = 1 + 0.1 * 2 = 1.2
//
// Because electron and ion are produced in equal numbers:
//
//     rho = q_e n_e + q_i n_i = 0
//
// both before and after the reaction step.
// ============================================================

TEST_F(FixedStepElectrostaticDriftDiffusionReactionTest,
       IonizationProducesNeutralElectronIonPairs) {
  // ========================================================
  // Physical / numerical parameters
  // ========================================================

  constexpr double permittivity = 1.0;

  constexpr double dt = 0.1;

  constexpr double initial_density = 1.0;

  constexpr double neutral_density = 4.0;

  constexpr double ionization_rate_coefficient = 0.5;

  // ========================================================
  // Species densities
  //
  // Initially quasi-neutral:
  //
  //     n_e = n_i = 1
  // ========================================================

  field::CellField<double> electron_density(mesh_, initial_density);

  field::CellField<double> ion_density(mesh_, initial_density);

  // ========================================================
  // Reaction rate and species sources
  // ========================================================

  physics::SpeciesSet reaction_species;
  const auto electron_id = reaction_species.add({.name = "e", .charge = -1.0});
  const auto ion_id = reaction_species.add({.name = "ion", .charge = +1.0});
  physics::reaction::ReactionNetwork reaction_network(reaction_species);
  const auto ionization = reaction_network.addReaction(
      {.name = "ionization",
       .stoichiometry = {{electron_id, +1.0}, {ion_id, +1.0}},
       .kinetic_orders = {{electron_id, 1.0}}});
  physics::reaction::ReactionRateFields reaction_rates(
      mesh_, reaction_network.size(), 0.0);
  physics::SpeciesCellFields species_source(mesh_, reaction_species.size(),
                                            0.0);
  auto& reaction_rate = reaction_rates[ionization];
  auto& electron_source = species_source[electron_id];
  auto& ion_source = species_source[ion_id];

  // ========================================================
  // Calculate:
  //
  //     R = k_ion n_e n_N
  // ========================================================

  physics::reaction::binaryReactionRate(
      electron_density, neutral_density, ionization_rate_coefficient,
      reaction_rate);

  // Expected rate:
  //
  //     0.5 * 1 * 4 = 2
  //
  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(reaction_rate[cell], 2.0, kTolerance);
  }

  // ========================================================
  // Convert reaction rate into species sources:
  //
  //     S_e += R
  //     S_i += R
  // ========================================================

  species_source.fill(0.0);
  reaction_network.accumulateSources(reaction_rates.span(), species_source);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(electron_source[cell], 2.0, kTolerance);

    EXPECT_NEAR(ion_source[cell], 2.0, kTolerance);
  }

  // ========================================================
  // Check reaction charge source before stepping:
  //
  //     q_e S_e + q_i S_i = 0
  // ========================================================

  constexpr double electron_charge = -1.0;

  constexpr double ion_charge = +1.0;

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    const double charge_source =
        electron_charge * electron_source[cell] + ion_charge * ion_source[cell];

    EXPECT_NEAR(charge_source, 0.0, kTolerance);
  }

  // ========================================================
  // Potential boundary conditions
  //
  // Homogeneous Dirichlet:
  //
  //     phi = 0
  //
  // Since plasma starts neutral:
  //
  //     rho = 0
  //
  // therefore:
  //
  //     phi = 0
  //     E   = 0
  // ========================================================

  boundary::BoundaryConditionSet potential_bc;

  potential_bc.setDirichlet(mesh::BoundaryId{1}, 0.0);

  potential_bc.setDirichlet(mesh::BoundaryId{2}, 0.0);

  potential_bc.setDirichlet(mesh::BoundaryId{3}, 0.0);

  potential_bc.setDirichlet(mesh::BoundaryId{4}, 0.0);

  // ========================================================
  // Species boundary conditions
  //
  // At the beginning of this step:
  //
  //     n_e = n_i = 1
  //
  // Therefore setting boundary state to 1 makes transport
  // flux zero initially.
  //
  // This deliberately isolates the reaction source term.
  // ========================================================

  boundary::BoundaryConditionSet electron_bc;

  boundary::BoundaryConditionSet ion_bc;

  for (const mesh::BoundaryId id : {mesh::BoundaryId{1}, mesh::BoundaryId{2},
                                    mesh::BoundaryId{3}, mesh::BoundaryId{4}}) {

    electron_bc.setDirichlet(id, initial_density);

    ion_bc.setDirichlet(id, initial_density);
  }

  // ========================================================
  // Species transport parameters
  //
  // Mobility does not matter in this particular test because
  // E = 0.
  //
  // Diffusivity also produces zero flux because initial
  // density is spatially constant.
  // ========================================================

  physics::ChargedSpeciesTransport electron{.charge = electron_charge,

                                            .mobility = 1.0,

                                            .diffusivity = 0.1};

  physics::ChargedSpeciesTransport ion{.charge = ion_charge,

                                       .mobility = 0.5,

                                       .diffusivity = 0.1};

  // ========================================================
  // Construct coupled electrostatic drift-diffusion stepper
  // ========================================================

  equation::FixedStepElectrostaticDriftDiffusionStepper stepper(
      mesh_, permittivity, dt, electron, ion, potential_bc, electron_bc, ion_bc,
      std::make_unique<linalg::CholmodSolver>());

  // ========================================================
  // Particle inventory BEFORE timestep
  // ========================================================

  const double electron_before = totalParticles(mesh_, electron_density);

  const double ion_before = totalParticles(mesh_, ion_density);

  // ========================================================
  // Integrated reaction rate:
  //
  //     integral R dV
  // ========================================================

  double integrated_reaction_rate = 0.0;

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    integrated_reaction_rate += reaction_rate[cell] * mesh_.cellVolume(cell);
  }

  // ========================================================
  // Advance one coupled timestep
  //
  // Internally:
  //
  // n_e^k, n_i^k
  //       ↓
  //      rho^k
  //       ↓
  //   Poisson solve
  //       ↓
  //      phi^k
  //       ↓
  //       E^k
  //       ↓
  //   drift velocity
  //       ↓
  //      SG flux
  //       ↓
  // continuity + source
  //       ↓
  // n_e^(k+1), n_i^(k+1)
  // ========================================================

  const auto result =
      stepper.step(electron_density, electron_source, ion_density, ion_source);

  ASSERT_TRUE(result.success());

  // ========================================================
  // Local density verification
  //
  //     n^(k+1)
  //
  //       = 1 + dt * R
  //
  //       = 1 + 0.1 * 2
  //
  //       = 1.2
  // ========================================================

  constexpr double expected_density = 1.2;

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(electron_density[cell], expected_density, kTolerance);

    EXPECT_NEAR(ion_density[cell], expected_density, kTolerance);
  }

  // ========================================================
  // Global particle balance
  // ========================================================

  const double electron_after = totalParticles(mesh_, electron_density);

  const double ion_after = totalParticles(mesh_, ion_density);

  const double expected_particle_increase = dt * integrated_reaction_rate;

  EXPECT_NEAR(electron_after - electron_before, expected_particle_increase,
              kTolerance);

  EXPECT_NEAR(ion_after - ion_before, expected_particle_increase, kTolerance);

  // ========================================================
  // Electron and ion production must be identical.
  // ========================================================

  EXPECT_NEAR(electron_after - electron_before, ion_after - ion_before,
              kTolerance);

  // ========================================================
  // IMPORTANT:
  //
  // step() used rho^k / phi^k / E^k to advance species.
  //
  // electron_density / ion_density are now k+1, but the
  // stepper's electrostatic fields still correspond to k.
  //
  // Recompute electrostatics so that:
  //
  //     rho, phi, E
  //
  // correspond to the newly updated densities.
  // ========================================================

  const auto electrostatic_result =
      stepper.updateElectrostatics(electron_density, ion_density);

  ASSERT_TRUE(electrostatic_result.success());

  // ========================================================
  // Pair production must not create net charge:
  //
  //     rho^(k+1)
  //
  //       = -n_e^(k+1)
  //         +n_i^(k+1)
  //
  //       = 0
  // ========================================================

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(stepper.chargeDensity()[cell], 0.0, kTolerance);

    EXPECT_NEAR(stepper.potential()[cell], 0.0, kTolerance);
  }

  // ========================================================
  // Consequently electric field and drift velocity remain
  // zero after electrostatics is synchronized to k+1.
  // ========================================================

  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {
    EXPECT_NEAR(stepper.electricFieldNormal()[face], 0.0, kTolerance);
    EXPECT_NEAR(stepper.electronDriftVelocityNormal()[face], 0.0, kTolerance);
    EXPECT_NEAR(stepper.ionDriftVelocityNormal()[face], 0.0, kTolerance);
  }
}

}  // namespace pemu::equation::test
