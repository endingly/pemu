#include <gtest/gtest.h>
#include <pemu/boundary/boundary_condition_set.hpp>
#include <pemu/equation/adaptive_step_multi_species_drift_diffusion_stepper.hpp>
#include <pemu/equation/fixed_step_electrostatic_drift_diffusion_stepper.hpp>
#include <pemu/equation/fixed_step_multi_species_drift_diffusion_stepper.hpp>
#include <pemu/field/cell_field.hpp>
#include <pemu/linalg/cholmod_solver.hpp>
#include <pemu/mesh/moab_mesh.hpp>
#include <pemu/physics/charged_species_transport.hpp>
#include <pemu/physics/species.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <vector>

namespace pemu::equation::test {

namespace {

// ============================================================
// Mesh helpers
// ============================================================

[[nodiscard]]
std::filesystem::path twoQuadsMeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
}

[[nodiscard]]
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
// Boundary-condition helpers
// ============================================================

[[nodiscard]]
boundary::BoundaryConditionSet makeZeroPotentialBoundaryConditions() {
  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(mesh::BoundaryId{1}, 0.0);

  bc.setDirichlet(mesh::BoundaryId{2}, 0.0);

  bc.setDirichlet(mesh::BoundaryId{3}, 0.0);

  bc.setDirichlet(mesh::BoundaryId{4}, 0.0);

  return bc;
}

// ------------------------------------------------------------
// Construct BC corresponding to:
//
//     phi(x,y) = x
//
// Therefore:
//
//     E = -grad(phi) = (-1, 0)
//
// Left/right:
//     Dirichlet phi = x
//
// Bottom/top:
//     Neumann
//
//     -epsilon grad(phi) · n = 0
// ------------------------------------------------------------

[[nodiscard]]
boundary::BoundaryConditionSet makeLinearXPotentialBoundaryConditions(
    const mesh::IMesh& mesh) {
  boundary::BoundaryConditionSet bc;

  const auto left_face = findBoundaryFace(mesh, mesh::BoundaryId{1});

  const auto right_face = findBoundaryFace(mesh, mesh::BoundaryId{2});

  bc.setDirichlet(mesh::BoundaryId{1}, mesh.faceCenter(left_face).x);

  bc.setDirichlet(mesh::BoundaryId{2}, mesh.faceCenter(right_face).x);

  bc.setNeumann(mesh::BoundaryId{3}, 0.0);

  bc.setNeumann(mesh::BoundaryId{4}, 0.0);

  return bc;
}

// ------------------------------------------------------------
// Create one BoundaryConditionSet for each species.
//
// This is deliberately indexed exactly like SpeciesId.
// ------------------------------------------------------------

[[nodiscard]]
std::vector<boundary::BoundaryConditionSet>
makeConstantSpeciesBoundaryConditions(std::size_t species_count, double value) {
  std::vector<boundary::BoundaryConditionSet> result(species_count);

  for (auto& bc : result) {

    bc.setDirichlet(mesh::BoundaryId{1}, value);

    bc.setDirichlet(mesh::BoundaryId{2}, value);

    bc.setDirichlet(mesh::BoundaryId{3}, value);

    bc.setDirichlet(mesh::BoundaryId{4}, value);
  }

  return result;
}

// ============================================================
// Field helpers
// ============================================================

[[nodiscard]]
std::vector<double> snapshot(const field::CellField<double>& field) {
  return {field.begin(), field.end()};
}

void expectFieldEqualsSnapshot(const field::CellField<double>& field,
                               const std::vector<double>& expected,
                               double tolerance = 1e-12) {
  ASSERT_EQ(field.size(), expected.size());

  for (std::size_t i = 0; i < expected.size(); ++i) {

    EXPECT_NEAR(field[static_cast<mesh::CellId>(i)], expected[i], tolerance);
  }
}

// ============================================================
// Common two-species definition
// ============================================================

struct ElectronIonIds {
  physics::SpeciesId electron;
  physics::SpeciesId ion;
};

[[nodiscard]]
ElectronIonIds addElectronAndIon(physics::SpeciesSet& species,
                                 double electron_mobility = 1.0,
                                 double ion_mobility = 0.5,
                                 double electron_diffusivity = 0.1,
                                 double ion_diffusivity = 0.1) {
  const auto electron = species.add(
      {.name = "e",
       .charge = -1.0,
       .mobility = electron_mobility,
       .diffusivity = electron_diffusivity,
       .transport_model = physics::SpeciesTransportModel::DriftDiffusion});

  const auto ion = species.add(
      {.name = "Ar+",
       .charge = +1.0,
       .mobility = ion_mobility,
       .diffusivity = ion_diffusivity,
       .transport_model = physics::SpeciesTransportModel::DriftDiffusion});

  return {.electron = electron, .ion = ion};
}

}  // namespace

// ============================================================
// Fixture
// ============================================================

class FixedStepMultiSpeciesDriftDiffusionStepperTest : public ::testing::Test {
 protected:
  FixedStepMultiSpeciesDriftDiffusionStepperTest()
      : mesh_(twoQuadsMeshPath()) {}

  mesh::MoabMesh mesh_;
};

class AdaptiveStepMultiSpeciesDriftDiffusionStepperTest
    : public ::testing::Test {
 protected:
  AdaptiveStepMultiSpeciesDriftDiffusionStepperTest()
      : mesh_(twoQuadsMeshPath()) {}

  mesh::MoabMesh mesh_;
};

// ============================================================
// 1. Regression:
//    uniform neutral electron-ion plasma remains stationary.
//
//    n_e = n_i
//
//    rho = -n_e + n_i = 0
//
//    phi = 0
//    E   = 0
//    v_e = 0
//    v_i = 0
//
//    constant density -> zero SG flux.
// ============================================================

TEST_F(FixedStepMultiSpeciesDriftDiffusionStepperTest,
       UniformNeutralPlasmaRemainsStationary) {
  constexpr double density_value = 3.0;

  // --------------------------------------------------------
  // Species definition
  // --------------------------------------------------------

  physics::SpeciesSet species;

  const auto ids = addElectronAndIon(species);

  // --------------------------------------------------------
  // n_s
  // --------------------------------------------------------

  physics::SpeciesCellFields density(mesh_, species.size(), density_value);

  // --------------------------------------------------------
  // S_s = 0
  // --------------------------------------------------------

  physics::SpeciesCellFields source(mesh_, species.size(), 0.0);

  // --------------------------------------------------------
  // Boundary conditions
  // --------------------------------------------------------

  auto potential_bc = makeZeroPotentialBoundaryConditions();

  auto species_bc =
      makeConstantSpeciesBoundaryConditions(species.size(), density_value);

  // --------------------------------------------------------
  // Solver
  // --------------------------------------------------------

  FixedStepMultiSpeciesDriftDiffusionStepper stepper(
      mesh_, species,

      1.0,   // epsilon
      0.01,  // dt

      std::move(potential_bc),

      std::move(species_bc),

      std::make_unique<linalg::CholmodSolver>());

  // --------------------------------------------------------
  // Advance
  // --------------------------------------------------------

  const auto result = stepper.step(density, source);

  ASSERT_TRUE(result.success());

  // --------------------------------------------------------
  // rho = 0
  // phi = 0
  // densities unchanged
  // --------------------------------------------------------

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(stepper.chargeDensity()[cell], 0.0, 1e-12);

    EXPECT_NEAR(stepper.potential()[cell], 0.0, 1e-12);

    EXPECT_NEAR(density[ids.electron][cell], density_value, 1e-12);

    EXPECT_NEAR(density[ids.ion][cell], density_value, 1e-12);
  }

  // --------------------------------------------------------
  // E = 0
  // v_e = v_i = 0
  // --------------------------------------------------------

  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {

    EXPECT_NEAR(stepper.electricFieldNormal()[face], 0.0, 1e-12);

    EXPECT_NEAR(stepper.driftVelocityNormal(ids.electron)[face], 0.0, 1e-12);

    EXPECT_NEAR(stepper.driftVelocityNormal(ids.ion)[face], 0.0, 1e-12);
  }
}

TEST_F(FixedStepMultiSpeciesDriftDiffusionStepperTest,
       PureNeumannGaugeIsAvailableFromCoupledStepper) {
  physics::SpeciesSet species;
  [[maybe_unused]] const auto ids = addElectronAndIon(species);
  physics::SpeciesCellFields density(mesh_, species.size(), 3.0);

  boundary::BoundaryConditionSet potential_bc;
  potential_bc.setNeumann(mesh::BoundaryId{1}, 0.0);
  potential_bc.setNeumann(mesh::BoundaryId{2}, 0.0);
  potential_bc.setNeumann(mesh::BoundaryId{3}, 0.0);
  potential_bc.setNeumann(mesh::BoundaryId{4}, 0.0);

  auto species_bc = makeConstantSpeciesBoundaryConditions(species.size(), 3.0);
  FixedStepMultiSpeciesDriftDiffusionStepper stepper(
      mesh_, species, 1.0, 0.01, std::move(potential_bc), std::move(species_bc),
      std::make_unique<linalg::CholmodSolver>(), field::PlasmaFieldMetadata{},
      PureNeumannOptions{.gauge = PinCellGauge{}});

  const auto result = stepper.updateElectrostatics(density);

  ASSERT_TRUE(result.success());
  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {
    EXPECT_NEAR(stepper.potential()[cell], 0.0, 1e-12);
  }
  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {
    EXPECT_NEAR(stepper.electricFieldNormal()[face], 0.0, 1e-12);
  }
}

// ============================================================
// 2. The electrostatic solve must support more than two
//    species.
//
//        e      q = -1, n = 3
//        Ar+    q = +1, n = 1
//        Ar2+   q = +2, n = 1
//
//    Therefore:
//
//        rho = -3 + 1 + 2 = 0
//
//    All three are marked Immobile because this test only
//    verifies:
//
//        SpeciesSet -> rho -> Poisson
//
//    and deliberately does not exercise transport.
// ============================================================

TEST_F(FixedStepMultiSpeciesDriftDiffusionStepperTest,
       ChargeDensitySupportsMoreThanTwoSpecies) {
  physics::SpeciesSet species;

  const auto electron = species.add(
      {.name = "e",
       .charge = -1.0,
       .mobility = 0.0,
       .diffusivity = 0.0,
       .transport_model = physics::SpeciesTransportModel::Immobile});

  const auto ion1 = species.add(
      {.name = "Ar+",
       .charge = +1.0,
       .mobility = 0.0,
       .diffusivity = 0.0,
       .transport_model = physics::SpeciesTransportModel::Immobile});

  const auto ion2 = species.add(
      {.name = "Ar2+",
       .charge = +2.0,
       .mobility = 0.0,
       .diffusivity = 0.0,
       .transport_model = physics::SpeciesTransportModel::Immobile});

  physics::SpeciesCellFields density(mesh_, species.size(), 0.0);

  density[electron].fill(3.0);

  density[ion1].fill(1.0);

  density[ion2].fill(1.0);

  auto potential_bc = makeZeroPotentialBoundaryConditions();

  std::vector<boundary::BoundaryConditionSet> species_bc(species.size());

  FixedStepMultiSpeciesDriftDiffusionStepper stepper(
      mesh_, species,

      1.0, 0.01,

      std::move(potential_bc),

      std::move(species_bc),

      std::make_unique<linalg::CholmodSolver>());

  const auto result = stepper.updateElectrostatics(density);

  ASSERT_TRUE(result.success());

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(stepper.chargeDensity()[cell], 0.0, 1e-12);

    EXPECT_NEAR(stepper.potential()[cell], 0.0, 1e-12);
  }
}

// ============================================================
// 3. Immobile species must not be evolved by this solver.
//
//    Current semantic:
//
//        Immobile
//
//    means:
//
//        FixedStepMultiSpeciesDriftDiffusionStepper does not update
//        the species at all.
//
//    Even a non-zero source passed for that species is ignored
//    by this solver.
//
//    Chemistry-only evolution can be introduced separately.
// ============================================================

TEST_F(FixedStepMultiSpeciesDriftDiffusionStepperTest,
       ImmobileSpeciesIsNotUpdated) {
  physics::SpeciesSet species;

  const auto ids = addElectronAndIon(species);

  const auto neutral = species.add(
      {.name = "Ar",
       .charge = 0.0,
       .mobility = 0.0,
       .diffusivity = 0.0,
       .transport_model = physics::SpeciesTransportModel::Immobile});

  physics::SpeciesCellFields density(mesh_, species.size(), 0.0);

  density[ids.electron].fill(3.0);

  density[ids.ion].fill(3.0);

  density[neutral].fill(100.0);

  physics::SpeciesCellFields source(mesh_, species.size(), 0.0);

  //
  // Deliberately non-zero.
  //
  source[neutral].fill(50.0);

  auto potential_bc = makeZeroPotentialBoundaryConditions();

  auto species_bc = makeConstantSpeciesBoundaryConditions(species.size(), 3.0);

  FixedStepMultiSpeciesDriftDiffusionStepper stepper(
      mesh_, species,

      1.0, 0.01,

      std::move(potential_bc),

      std::move(species_bc),

      std::make_unique<linalg::CholmodSolver>());

  const auto before = snapshot(density[neutral]);

  const auto result = stepper.step(density, source);

  ASSERT_TRUE(result.success());

  expectFieldEqualsSnapshot(density[neutral], before);
}

// ============================================================
// 4. Applied potential:
//
//        phi = x
//
//    gives:
//
//        E = (-1, 0)
//
//    Electron:
//
//        q < 0
//        mu = 2
//
//        v_e = -mu E
//            = (+2, 0)
//
//    Positive ion:
//
//        q > 0
//        mu = 1
//
//        v_i = +mu E
//            = (-1, 0)
//
//    On a face:
//
//        E_n = -n_x
//        v_e,n = +2 n_x
//        v_i,n = -1 n_x
// ============================================================

TEST_F(FixedStepMultiSpeciesDriftDiffusionStepperTest,
       AppliedPotentialProducesCorrectDriftForAllChargedSpecies) {
  physics::SpeciesSet species;

  const auto ids = addElectronAndIon(species,

                                     2.0,  // electron mobility
                                     1.0,  // ion mobility

                                     0.1, 0.1);

  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);

  auto potential_bc = makeLinearXPotentialBoundaryConditions(mesh_);

  auto species_bc = makeConstantSpeciesBoundaryConditions(species.size(), 1.0);

  FixedStepMultiSpeciesDriftDiffusionStepper stepper(
      mesh_, species,

      1.0, 0.001,

      std::move(potential_bc),

      std::move(species_bc),

      std::make_unique<linalg::CholmodSolver>());

  const auto result = stepper.updateElectrostatics(density);

  ASSERT_TRUE(result.success());

  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {

    const auto normal = mesh_.faceNormal(face);

    // ----------------------------------------------------
    // phi = x
    //
    // E = (-1, 0)
    //
    // E_n = -n_x
    // ----------------------------------------------------

    EXPECT_NEAR(stepper.electricFieldNormal()[face],

                -normal.x,

                1e-12);

    // ----------------------------------------------------
    // Electron:
    //
    // v = (+2, 0)
    // ----------------------------------------------------

    EXPECT_NEAR(stepper.driftVelocityNormal(ids.electron)[face],

                +2.0 * normal.x,

                1e-12);

    // ----------------------------------------------------
    // Positive ion:
    //
    // v = (-1, 0)
    // ----------------------------------------------------

    EXPECT_NEAR(stepper.driftVelocityNormal(ids.ion)[face],

                -1.0 * normal.x,

                1e-12);
  }
}

// ============================================================
// 5. Generalization test with THREE transported charged species.
//
//    This verifies that drift computation really iterates over
//    SpeciesSet instead of containing electron/ion special
//    cases.
//
//    phi = x -> E = (-1, 0)
//
//      e:    q < 0, mu=2    -> v=(+2,0)
//      Ar+:  q > 0, mu=1    -> v=(-1,0)
//      Ar2+: q > 0, mu=.25  -> v=(-.25,0)
//
//    Note:
//
//    charge magnitude does NOT multiply mobility here.
//
//        v_s = sign(q_s) mu_s E
//
//    The species-specific mobility already contains the
//    transport response.
// ============================================================

TEST_F(FixedStepMultiSpeciesDriftDiffusionStepperTest,
       ComputesDriftForThreeTransportedChargedSpecies) {
  physics::SpeciesSet species;

  const auto electron = species.add(
      {.name = "e",
       .charge = -1.0,
       .mobility = 2.0,
       .diffusivity = 0.1,
       .transport_model = physics::SpeciesTransportModel::DriftDiffusion});

  const auto ion1 = species.add(
      {.name = "Ar+",
       .charge = +1.0,
       .mobility = 1.0,
       .diffusivity = 0.1,
       .transport_model = physics::SpeciesTransportModel::DriftDiffusion});

  const auto ion2 = species.add(
      {.name = "Ar2+",
       .charge = +2.0,
       .mobility = 0.25,
       .diffusivity = 0.1,
       .transport_model = physics::SpeciesTransportModel::DriftDiffusion});

  physics::SpeciesCellFields density(mesh_, species.size(), 0.0);

  //
  // Neutral total charge:
  //
  // -3 + 1 + 2 = 0
  //
  density[electron].fill(3.0);

  density[ion1].fill(1.0);

  density[ion2].fill(1.0);

  auto potential_bc = makeLinearXPotentialBoundaryConditions(mesh_);

  auto species_bc = makeConstantSpeciesBoundaryConditions(species.size(), 1.0);

  FixedStepMultiSpeciesDriftDiffusionStepper stepper(
      mesh_, species,

      1.0, 0.001,

      std::move(potential_bc),

      std::move(species_bc),

      std::make_unique<linalg::CholmodSolver>());

  const auto result = stepper.updateElectrostatics(density);

  ASSERT_TRUE(result.success());

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(stepper.chargeDensity()[cell], 0.0, 1e-12);
  }

  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {

    const auto nx = mesh_.faceNormal(face).x;

    EXPECT_NEAR(stepper.driftVelocityNormal(electron)[face], +2.0 * nx, 1e-12);

    EXPECT_NEAR(stepper.driftVelocityNormal(ion1)[face], -1.0 * nx, 1e-12);

    EXPECT_NEAR(stepper.driftVelocityNormal(ion2)[face], -0.25 * nx, 1e-12);
  }
}

// ============================================================
// 6. CFL failure must be atomic.
//
//    Critical invariant:
//
//        validate/update field
//        -> calculate all CFL values
//        -> only then modify densities
//
//    If any species violates its explicit transport CFL,
//    NO species density may be changed.
//
//    A deliberately huge dt is used here.
// ============================================================

TEST_F(FixedStepMultiSpeciesDriftDiffusionStepperTest,
       CflFailureDoesNotPartiallyUpdateSpecies) {
  physics::SpeciesSet species;

  const auto ids = addElectronAndIon(species,

                                     2.0, 1.0,

                                     0.1, 0.1);

  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);

  // Make fields non-identical so accidental writes are easier
  // to detect.
  density[ids.electron][0] = 2.0;

  density[ids.ion][1] = 2.0;

  physics::SpeciesCellFields source(mesh_, species.size(), 0.0);

  auto potential_bc = makeLinearXPotentialBoundaryConditions(mesh_);

  auto species_bc = makeConstantSpeciesBoundaryConditions(species.size(), 1.0);

  //
  // Deliberately enormous for two_quads.msh.
  //
  constexpr double dt = 100.0;

  FixedStepMultiSpeciesDriftDiffusionStepper stepper(
      mesh_, species,

      1.0, dt,

      std::move(potential_bc),

      std::move(species_bc),

      std::make_unique<linalg::CholmodSolver>());

  const auto electron_before = snapshot(density[ids.electron]);

  const auto ion_before = snapshot(density[ids.ion]);

  EXPECT_THROW(auto _ = stepper.step(density, source), std::runtime_error);

  // --------------------------------------------------------
  // Neither species may have been modified.
  // --------------------------------------------------------

  expectFieldEqualsSnapshot(density[ids.electron], electron_before);

  expectFieldEqualsSnapshot(density[ids.ion], ion_before);
}

// ============================================================
// 7. SpeciesCellFields with the wrong number of species must
//    be rejected.
//
//    SpeciesSet is part of the solver's structural definition;
//    fields indexed by another species layout are not legal.
// ============================================================

TEST_F(FixedStepMultiSpeciesDriftDiffusionStepperTest,
       RejectsWrongSpeciesFieldCount) {
  physics::SpeciesSet species;

  auto _ = addElectronAndIon(species);

  auto potential_bc = makeZeroPotentialBoundaryConditions();

  auto species_bc = makeConstantSpeciesBoundaryConditions(species.size(), 1.0);

  FixedStepMultiSpeciesDriftDiffusionStepper stepper(
      mesh_, species,

      1.0, 0.01,

      std::move(potential_bc),

      std::move(species_bc),

      std::make_unique<linalg::CholmodSolver>());

  //
  // Solver expects two fields.
  //
  // Deliberately provide three.
  //
  physics::SpeciesCellFields wrong_density(mesh_, 3, 1.0);

  EXPECT_THROW(auto _ = stepper.updateElectrostatics(wrong_density),
               std::invalid_argument);
}

// ============================================================
// 8. Fields from another mesh object must be rejected.
//
//    Even if the two meshes contain identical topology,
//    our runtime invariant is object identity:
//
//        &fields.mesh() == solver.mesh
// ============================================================

TEST_F(FixedStepMultiSpeciesDriftDiffusionStepperTest,
       RejectsFieldsFromDifferentMesh) {
  physics::SpeciesSet species;

  auto _ = addElectronAndIon(species);

  auto potential_bc = makeZeroPotentialBoundaryConditions();

  auto species_bc = makeConstantSpeciesBoundaryConditions(species.size(), 1.0);

  FixedStepMultiSpeciesDriftDiffusionStepper stepper(
      mesh_, species,

      1.0, 0.01,

      std::move(potential_bc),

      std::move(species_bc),

      std::make_unique<linalg::CholmodSolver>());

  mesh::MoabMesh other_mesh(twoQuadsMeshPath());

  physics::SpeciesCellFields density(other_mesh, species.size(), 1.0);

  EXPECT_THROW(auto _ = stepper.updateElectrostatics(density),
               std::invalid_argument);
}

// ============================================================
// 9. Optional but strongly recommended:
//
//    New multi-species solver must reproduce the old verified
//    two-species implementation.
//
//    Once this test passes, the old solver can eventually be
//    removed.
//
//    This test assumes the old interfaces are still:
//
//      ChargedSpeciesTransport
//      FixedStepElectrostaticDriftDiffusionStepper
//
//    as implemented in the previous stage.
// ============================================================

TEST_F(FixedStepMultiSpeciesDriftDiffusionStepperTest,
       MatchesLegacyTwoSpeciesSolver) {
  constexpr double dt = 1e-3;

  constexpr double epsilon = 1.0;

  // ========================================================
  // Shared initial state
  // ========================================================

  field::CellField<double> old_electron_density(mesh_, 0.0);

  field::CellField<double> old_ion_density(mesh_, 0.0);

  //
  // Deliberately nonuniform.
  //
  old_electron_density[0] = 1.0;

  old_electron_density[1] = 2.0;

  old_ion_density[0] = 1.0;

  old_ion_density[1] = 2.0;

  field::CellField<double> old_electron_source(mesh_, 0.0);

  field::CellField<double> old_ion_source(mesh_, 0.0);

  // ========================================================
  // Shared BC
  // ========================================================

  auto old_potential_bc = makeLinearXPotentialBoundaryConditions(mesh_);

  boundary::BoundaryConditionSet old_electron_bc;

  boundary::BoundaryConditionSet old_ion_bc;

  for (const auto id : {mesh::BoundaryId{1}, mesh::BoundaryId{2},
                        mesh::BoundaryId{3}, mesh::BoundaryId{4}}) {

    old_electron_bc.setDirichlet(id, 1.0);

    old_ion_bc.setDirichlet(id, 1.0);
  }

  // ========================================================
  // Legacy solver
  // ========================================================

  physics::ChargedSpeciesTransport old_electron{
      .charge = -1.0, .mobility = 2.0, .diffusivity = 0.1};

  physics::ChargedSpeciesTransport old_ion{
      .charge = +1.0, .mobility = 1.0, .diffusivity = 0.1};

  FixedStepElectrostaticDriftDiffusionStepper old_solver(
      mesh_, epsilon, dt,

      old_electron, old_ion,

      old_potential_bc, old_electron_bc, old_ion_bc,

      std::make_unique<linalg::CholmodSolver>());

  // ========================================================
  // New multi-species solver
  // ========================================================

  physics::SpeciesSet species;

  const auto ids = addElectronAndIon(species,

                                     2.0, 1.0,

                                     0.1, 0.1);

  physics::SpeciesCellFields new_density(mesh_, species.size(), 0.0);

  new_density[ids.electron][0] = 1.0;

  new_density[ids.electron][1] = 2.0;

  new_density[ids.ion][0] = 1.0;

  new_density[ids.ion][1] = 2.0;

  physics::SpeciesCellFields new_source(mesh_, species.size(), 0.0);

  auto new_potential_bc = makeLinearXPotentialBoundaryConditions(mesh_);

  auto new_species_bc =
      makeConstantSpeciesBoundaryConditions(species.size(), 1.0);

  FixedStepMultiSpeciesDriftDiffusionStepper new_solver(
      mesh_, species,

      epsilon, dt,

      std::move(new_potential_bc),

      std::move(new_species_bc),

      std::make_unique<linalg::CholmodSolver>());

  // ========================================================
  // Advance both implementations exactly one timestep.
  // ========================================================

  const auto old_result =
      old_solver.step(old_electron_density, old_electron_source,

                      old_ion_density, old_ion_source);

  const auto new_result = new_solver.step(new_density, new_source);

  ASSERT_TRUE(old_result.success());

  ASSERT_TRUE(new_result.success());

  // ========================================================
  // Compare cell quantities.
  // ========================================================

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(old_electron_density[cell],

                new_density[ids.electron][cell],

                1e-12);

    EXPECT_NEAR(old_ion_density[cell],

                new_density[ids.ion][cell],

                1e-12);

    EXPECT_NEAR(old_solver.chargeDensity()[cell],

                new_solver.chargeDensity()[cell],

                1e-12);

    EXPECT_NEAR(old_solver.potential()[cell],

                new_solver.potential()[cell],

                1e-12);
  }

  // ========================================================
  // Compare face quantities.
  // ========================================================

  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {

    EXPECT_NEAR(old_solver.electricFieldNormal()[face],

                new_solver.electricFieldNormal()[face],

                1e-12);

    EXPECT_NEAR(old_solver.electronDriftVelocityNormal()[face],

                new_solver.driftVelocityNormal(ids.electron)[face],

                1e-12);

    EXPECT_NEAR(old_solver.ionDriftVelocityNormal()[face],

                new_solver.driftVelocityNormal(ids.ion)[face],

                1e-12);
  }
}

TEST_F(AdaptiveStepMultiSpeciesDriftDiffusionStepperTest,
       UsesConfiguredMaximumTimeStep) {
  physics::SpeciesSet species;
  auto _ = addElectronAndIon(species);
  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);
  physics::SpeciesCellFields source(mesh_, species.size(), 0.0);

  AdaptiveStepMultiSpeciesDriftDiffusionStepper stepper(
      mesh_, species, 1.0, makeZeroPotentialBoundaryConditions(),
      makeConstantSpeciesBoundaryConditions(species.size(), 1.0),
      std::make_unique<linalg::CholmodSolver>(),
      {.safety = 0.9, .min_dt = 1e-8, .max_dt = 0.02, .max_growth = 2.0});

  ASSERT_TRUE(stepper.prepareElectrostatics(density).success());
  const auto proposal = stepper.advancePrepared(density, source, 0.1);

  EXPECT_NEAR(proposal.dt, 0.02, 1e-14);
  EXPECT_NEAR(stepper.previousTimeStep(), 0.02, 1e-14);
  EXPECT_TRUE(stepper.hasLastTimeStepProposal());
}

TEST_F(AdaptiveStepMultiSpeciesDriftDiffusionStepperTest,
       FinalTimeStepEqualsRemainingTime) {
  physics::SpeciesSet species;
  auto _ = addElectronAndIon(species);
  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);
  physics::SpeciesCellFields source(mesh_, species.size(), 0.0);

  AdaptiveStepMultiSpeciesDriftDiffusionStepper stepper(
      mesh_, species, 1.0, makeZeroPotentialBoundaryConditions(),
      makeConstantSpeciesBoundaryConditions(species.size(), 1.0),
      std::make_unique<linalg::CholmodSolver>(),
      {.safety = 0.9, .min_dt = 1e-8, .max_dt = 0.02, .max_growth = 2.0});

  ASSERT_TRUE(stepper.prepareElectrostatics(density).success());
  const auto proposal = stepper.advancePrepared(density, source, 0.005);

  EXPECT_NEAR(proposal.dt, 0.005, 1e-14);
}

}  // namespace pemu::equation::test
