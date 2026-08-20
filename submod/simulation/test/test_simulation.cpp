#include <gtest/gtest.h>
#include <pemu/boundary/boundary_condition_set.hpp>
#include <pemu/equation/adaptive_step_multi_species_drift_diffusion_stepper.hpp>
#include <pemu/equation/fixed_step_multi_species_drift_diffusion_stepper.hpp>
#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>
#include <pemu/linalg/cholmod_solver.hpp>
#include <pemu/linalg/i_solver.hpp>
#include <pemu/mesh/moab_mesh.hpp>
#include <pemu/physics/reaction.hpp>
#include <pemu/physics/species.hpp>
#include <pemu/simulation/adaptive_step_plasma_simulation.hpp>
#include <pemu/simulation/adaptive_time_clock.hpp>
#include <pemu/simulation/fixed_step_clock.hpp>
#include <pemu/simulation/fixed_step_plasma_simulation.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pemu::simulation::test {

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
// Exact potential:
//
//     phi(x,y) = x
//
// Therefore:
//
//     E = -grad(phi)
//       = (-1,0)
//
// left/right:
//     Dirichlet
//
// bottom/top:
//     -epsilon grad(phi) . n = 0
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
// Species helpers
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

// ============================================================
// Reaction evaluator used by the normal ionization tests.
//
//     R_ion
//
//       = k_ion * n_e * n_N
// ============================================================

struct ElectronImpactIonizationEvaluator {
  physics::SpeciesId electron;

  physics::ReactionId ionization;

  double neutral_density{};
  double rate_coefficient{};

  void operator()(const physics::SpeciesCellFields& density,

                  const field::CellField<double>&,

                  const field::FaceField<double>&,

                  physics::ReactionRateFields& reaction_rates) const {
    physics::reaction::electronImpactIonizationRate(
        density[electron], neutral_density, rate_coefficient,
        reaction_rates[ionization]);
  }
};

// ============================================================
// Evaluator used to prove:
//
//     reaction rate is recomputed every timestep.
//
// It records n_e[0] observed at each call.
// ============================================================

struct RecordingIonizationEvaluator {
  physics::SpeciesId electron;

  physics::ReactionId ionization;

  double neutral_density{};
  double rate_coefficient{};

  std::vector<double>* observed_electron_density{};

  void operator()(const physics::SpeciesCellFields& density,

                  const field::CellField<double>&,

                  const field::FaceField<double>&,

                  physics::ReactionRateFields& reaction_rates) const {
    if (observed_electron_density == nullptr) {

      throw std::logic_error("observation buffer is null");
    }

    observed_electron_density->push_back(density[electron][0]);

    physics::reaction::electronImpactIonizationRate(
        density[electron], neutral_density, rate_coefficient,
        reaction_rates[ionization]);
  }
};

// ============================================================
// Evaluator used to prove ordering:
//
//       updateElectrostatics()
//               ↓
//       reaction evaluation
//
// It simply fills the reaction rate with:
//
//       max_f |E_n,f|
//
// For phi=x:
//
//       E=(-1,0)
//
// so the maximum should be exactly 1.
// ============================================================

struct ElectricFieldAwareEvaluator {
  physics::ReactionId reaction;

  void operator()(const physics::SpeciesCellFields&,

                  const field::CellField<double>&,

                  const field::FaceField<double>& electric_field,

                  physics::ReactionRateFields& reaction_rates) const {
    double maximum = 0.0;

    for (const double value : electric_field) {

      maximum = std::max(maximum, std::abs(value));
    }

    reaction_rates[reaction].fill(maximum);
  }
};

// ============================================================
// Empty evaluator.
//
// Useful when testing the time loop or Poisson call count
// without involving chemistry.
// ============================================================

struct NoReactionEvaluator {
  void operator()(const physics::SpeciesCellFields&,

                  const field::CellField<double>&,

                  const field::FaceField<double>&,

                  physics::ReactionRateFields&) const noexcept {}
};

// ============================================================
// Counting linear-solver backend.
//
// This is intentionally minimal.
//
// For the tests using it:
//
//     rho = 0
//     phi boundary = 0
//
// therefore:
//
//     b = 0
//
// and x=b is already the correct Poisson solution.
//
// The point of this backend is NOT numerical correctness;
// it is verifying solver lifecycle:
//
//     analyze once
//     factorize once
//     solve once per timestep
// ============================================================

class CountingSolver final : public linalg::ISolver {
 public:
  int analyze_count{};
  int factorize_count{};
  int solve_count{};
  int reset_count{};

  linalg::SolverStatus analyzePattern(const linalg::SparseMatrix&) override {
    ++analyze_count;

    analyzed_ = true;

    return linalg::SolverStatus::Success;
  }

  linalg::SolverStatus factorize(const linalg::SparseMatrix&) override {
    ++factorize_count;

    factorized_ = true;

    return linalg::SolverStatus::Success;
  }

  linalg::SolverResult solve(linalg::ConstVectorRef b,
                             linalg::VectorRef x) override {
    ++solve_count;

    x = b;

    return {.status = linalg::SolverStatus::Success};
  }

  void reset() override {
    ++reset_count;

    analyzed_ = false;

    factorized_ = false;
  }

  [[nodiscard]]
  bool isAnalyzed() const noexcept override {
    return analyzed_;
  }

  [[nodiscard]]
  bool isFactorized() const noexcept override {
    return factorized_;
  }

 private:
  bool analyzed_{false};
  bool factorized_{false};
};

}  // namespace

// ============================================================
// Fixture
// ============================================================

class FixedStepPlasmaSimulationTest : public ::testing::Test {
 protected:
  FixedStepPlasmaSimulationTest() : mesh_(twoQuadsMeshPath()) {}

  mesh::MoabMesh mesh_;
};

class AdaptiveStepPlasmaSimulationTest : public ::testing::Test {
 protected:
  AdaptiveStepPlasmaSimulationTest() : mesh_(twoQuadsMeshPath()) {}

  mesh::MoabMesh mesh_;
};

// ============================================================
// FixedStepClock
// ============================================================

TEST(FixedStepClockTest, AdvancesTimeFromStepIndex) {
  FixedStepClock clock(0.1, 3);

  EXPECT_EQ(clock.step(), 0u);

  EXPECT_NEAR(clock.time(), 0.0, 1e-15);

  EXPECT_NEAR(clock.endTime(), 0.3, 1e-15);

  EXPECT_FALSE(clock.finished());

  clock.advance();

  EXPECT_EQ(clock.step(), 1u);

  EXPECT_NEAR(clock.time(), 0.1, 1e-15);

  EXPECT_FALSE(clock.finished());

  clock.advance();
  clock.advance();

  EXPECT_EQ(clock.step(), 3u);

  EXPECT_NEAR(clock.time(), 0.3, 1e-15);

  EXPECT_TRUE(clock.finished());
}

TEST(FixedStepClockTest, RejectsNonPositiveTimeStep) {
  EXPECT_THROW(FixedStepClock(0.0, 10), std::invalid_argument);

  EXPECT_THROW(FixedStepClock(-0.1, 10), std::invalid_argument);
}

TEST(FixedStepClockTest, RejectsAdvanceAfterCompletion) {
  FixedStepClock clock(0.1, 1);

  clock.advance();

  ASSERT_TRUE(clock.finished());

  EXPECT_THROW(clock.advance(), std::out_of_range);
}

// ============================================================
// 1. Complete single timestep:
//
//     n^k
//       ↓
//     electrostatics
//       ↓
//     R_ion
//       ↓
//     S = nu R
//       ↓
//     transport
//       ↓
//     n^(k+1)
//       ↓
//     clock.advance()
//
// Initial:
//
//     n_e = n_i = 1
//
// therefore rho=0 and E=0.
//
// Reaction:
//
//     k = 0.5
//     n_N = 4
//
//     R = k n_e n_N
//       = 0.5 * 1 * 4
//       = 2
//
// dt = 0.1
//
// therefore:
//
//     n_e^(1) = 1 + 0.1*2 = 1.2
//     n_i^(1) = 1 + 0.1*2 = 1.2
// ============================================================

TEST_F(FixedStepPlasmaSimulationTest,
       AdvanceOneStepEvaluatesReactionAndUpdatesSpecies) {
  constexpr double dt = 0.1;

  // --------------------------------------------------------
  // Species
  // --------------------------------------------------------

  physics::SpeciesSet species;

  const auto ids = addElectronAndIon(species);

  // --------------------------------------------------------
  // State
  // --------------------------------------------------------

  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);

  // --------------------------------------------------------
  // Reaction network
  // --------------------------------------------------------

  physics::ReactionNetwork reactions(species);

  const auto ionization = reactions.addReaction(
      {.name = "electron impact ionization",

       .stoichiometry = {{ids.electron, +1.0}, {ids.ion, +1.0}}});

  // --------------------------------------------------------
  // Transport solver
  // --------------------------------------------------------

  auto potential_bc = makeZeroPotentialBoundaryConditions();

  auto species_bc = makeConstantSpeciesBoundaryConditions(species.size(), 1.0);

  equation::FixedStepMultiSpeciesDriftDiffusionStepper transport(
      mesh_, species,

      1.0,  // epsilon
      dt,

      std::move(potential_bc),

      std::move(species_bc),

      std::make_unique<linalg::CholmodSolver>());

  // --------------------------------------------------------
  // Rate model
  // --------------------------------------------------------

  ElectronImpactIonizationEvaluator evaluator{.electron = ids.electron,

                                              .ionization = ionization,

                                              .neutral_density = 4.0,

                                              .rate_coefficient = 0.5};

  // --------------------------------------------------------
  // Simulation
  // --------------------------------------------------------

  FixedStepPlasmaSimulation simulation(density, reactions, transport, evaluator,

                                       FixedStepClock(dt, 1));

  // --------------------------------------------------------
  // Advance
  // --------------------------------------------------------

  const auto result = simulation.advanceOneStep();

  ASSERT_TRUE(result.success());

  // --------------------------------------------------------
  // Clock
  // --------------------------------------------------------

  EXPECT_EQ(simulation.step(), 1u);

  EXPECT_NEAR(simulation.time(), 0.1, 1e-14);

  EXPECT_TRUE(simulation.finished());

  // --------------------------------------------------------
  // Reaction + source + density
  // --------------------------------------------------------

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(simulation.reactionRates()[ionization][cell], 2.0, 1e-12);

    EXPECT_NEAR(simulation.source()[ids.electron][cell], 2.0, 1e-12);

    EXPECT_NEAR(simulation.source()[ids.ion][cell], 2.0, 1e-12);

    EXPECT_NEAR(density[ids.electron][cell], 1.2, 1e-12);

    EXPECT_NEAR(density[ids.ion][cell], 1.2, 1e-12);
  }
}

// ============================================================
// 2. Reaction rates MUST be reevaluated every timestep.
//
// We record electron density seen by the evaluator.
//
// First call:
//
//     n_e = 1
//
// after first timestep:
//
//     n_e = 1.2
//
// Therefore second evaluator invocation must observe the
// updated density.
//
// We deliberately do NOT assert the density after the second
// timestep, because fixed Dirichlet species BC may contribute
// transport once n differs from the original boundary value.
// ============================================================

TEST_F(FixedStepPlasmaSimulationTest,
       ReactionRateIsReevaluatedFromUpdatedStateEveryStep) {
  constexpr double dt = 0.1;

  physics::SpeciesSet species;

  const auto ids = addElectronAndIon(species);

  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);

  physics::ReactionNetwork reactions(species);

  const auto ionization = reactions.addReaction(
      {.name = "ionization",

       .stoichiometry = {{ids.electron, +1.0}, {ids.ion, +1.0}}});

  auto potential_bc = makeZeroPotentialBoundaryConditions();

  auto species_bc = makeConstantSpeciesBoundaryConditions(species.size(), 1.0);

  equation::FixedStepMultiSpeciesDriftDiffusionStepper transport(
      mesh_, species,

      1.0, dt,

      std::move(potential_bc),

      std::move(species_bc),

      std::make_unique<linalg::CholmodSolver>());

  std::vector<double> observations;

  RecordingIonizationEvaluator evaluator{
      .electron = ids.electron,

      .ionization = ionization,

      .neutral_density = 4.0,

      .rate_coefficient = 0.5,

      .observed_electron_density = &observations};

  FixedStepPlasmaSimulation simulation(density, reactions, transport, evaluator,

                                       FixedStepClock(dt, 2));

  ASSERT_TRUE(simulation.advanceOneStep().success());

  ASSERT_TRUE(simulation.advanceOneStep().success());

  ASSERT_EQ(observations.size(), 2u);

  EXPECT_NEAR(observations[0], 1.0, 1e-12);

  //
  // The first timestep produced:
  //
  //     n_e = 1.2
  //
  // Reaction evaluation in timestep 2 must see
  // that UPDATED state.
  //
  EXPECT_NEAR(observations[1], 1.2, 1e-12);
}

// ============================================================
// 3. Critical ordering test:
//
//     electrostatics
//          ↓
//     reaction evaluation
//
// must hold.
//
// Applied:
//
//     phi = x
//
// Therefore:
//
//     E = (-1,0)
//
// and:
//
//     max |E_n| = 1.
//
// ElectricFieldAwareEvaluator simply uses max|E_n| as the
// reaction rate.
//
// If reaction evaluation happened BEFORE Poisson/E update,
// this test would see the default zero electric field instead.
// ============================================================

TEST_F(FixedStepPlasmaSimulationTest,
       ReactionEvaluatorSeesCurrentElectricField) {
  constexpr double dt = 1e-3;

  physics::SpeciesSet species;

  const auto ids = addElectronAndIon(species);

  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);

  physics::ReactionNetwork reactions(species);

  const auto reaction = reactions.addReaction(
      {.name = "field dependent test reaction",

       .stoichiometry = {{ids.electron, +1.0}, {ids.ion, +1.0}}});

  auto potential_bc = makeLinearXPotentialBoundaryConditions(mesh_);

  auto species_bc = makeConstantSpeciesBoundaryConditions(species.size(), 1.0);

  equation::FixedStepMultiSpeciesDriftDiffusionStepper transport(
      mesh_, species,

      1.0, dt,

      std::move(potential_bc),

      std::move(species_bc),

      std::make_unique<linalg::CholmodSolver>());

  ElectricFieldAwareEvaluator evaluator{.reaction = reaction};

  FixedStepPlasmaSimulation simulation(density, reactions, transport, evaluator,

                                       FixedStepClock(dt, 1));

  const auto result = simulation.advanceOneStep();

  ASSERT_TRUE(result.success());

  // --------------------------------------------------------
  // phi=x -> E=(-1,0)
  //
  // Therefore the max absolute face-normal field is 1.
  // --------------------------------------------------------

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(simulation.reactionRates()[reaction][cell], 1.0, 1e-12);
  }

  // Also independently confirm the field.
  double maximum = 0.0;

  for (const double value : transport.electricFieldNormal()) {

    maximum = std::max(maximum, std::abs(value));
  }

  EXPECT_NEAR(maximum, 1.0, 1e-12);
}

// ============================================================
// 4. run() must advance until the configured number of steps
//    is reached.
//
// No reaction network is required for this lifecycle test.
// Uniform neutral plasma + zero potential remains stationary.
// ============================================================

TEST_F(FixedStepPlasmaSimulationTest, RunAdvancesUntilClockIsFinished) {
  constexpr double dt = 0.01;

  constexpr std::size_t total_steps = 5;

  physics::SpeciesSet species;

  auto _ = addElectronAndIon(species);

  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);

  physics::ReactionNetwork reactions(species);

  auto potential_bc = makeZeroPotentialBoundaryConditions();

  auto species_bc = makeConstantSpeciesBoundaryConditions(species.size(), 1.0);

  equation::FixedStepMultiSpeciesDriftDiffusionStepper transport(
      mesh_, species,

      1.0, dt,

      std::move(potential_bc),

      std::move(species_bc),

      std::make_unique<linalg::CholmodSolver>());

  FixedStepPlasmaSimulation simulation(density, reactions, transport,

                                       NoReactionEvaluator{},

                                       FixedStepClock(dt, total_steps));

  simulation.run();

  EXPECT_TRUE(simulation.finished());

  EXPECT_EQ(simulation.step(), total_steps);

  EXPECT_NEAR(simulation.time(), dt * static_cast<double>(total_steps), 1e-14);
}

// ============================================================
// 5. Exactly one Poisson solve per timestep.
//
// Furthermore:
//
//     analyzePattern = once
//     factorize      = once
//     solve          = once per timestep
//
// because Poisson matrix stays fixed while rho changes.
//
// This is one of the main reasons we separated:
//
//     updateElectrostatics()
//
// from:
//
//     advanceTransport()
// ============================================================

TEST_F(FixedStepPlasmaSimulationTest,
       SolvesPoissonExactlyOncePerTimeStepAndReusesFactorization) {
  constexpr double dt = 0.001;

  constexpr std::size_t total_steps = 3;

  physics::SpeciesSet species;

  auto _ = addElectronAndIon(species);

  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);

  physics::ReactionNetwork reactions(species);

  auto potential_bc = makeZeroPotentialBoundaryConditions();

  auto species_bc = makeConstantSpeciesBoundaryConditions(species.size(), 1.0);

  auto backend = std::make_unique<CountingSolver>();

  auto* counting_solver = backend.get();

  equation::FixedStepMultiSpeciesDriftDiffusionStepper transport(
      mesh_, species,

      1.0, dt,

      std::move(potential_bc),

      std::move(species_bc),

      std::move(backend));

  FixedStepPlasmaSimulation simulation(density, reactions, transport,

                                       NoReactionEvaluator{},

                                       FixedStepClock(dt, total_steps));

  simulation.run();

  EXPECT_EQ(counting_solver->analyze_count, 1);

  EXPECT_EQ(counting_solver->factorize_count, 1);

  EXPECT_EQ(counting_solver->solve_count, static_cast<int>(total_steps));
}

// ============================================================
// 6. Simulation clock and transport solver MUST agree on dt.
//
// Otherwise:
//
//     clock says dt_1
//
// while:
//
//     continuity update actually uses dt_2.
//
// Such a simulation would have no well-defined physical time.
// ============================================================

TEST_F(FixedStepPlasmaSimulationTest,
       RejectsClockTimeStepDifferentFromTransportTimeStep) {
  constexpr double transport_dt = 0.01;

  constexpr double clock_dt = 0.02;

  physics::SpeciesSet species;

  auto _ = addElectronAndIon(species);

  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);

  physics::ReactionNetwork reactions(species);

  auto potential_bc = makeZeroPotentialBoundaryConditions();

  auto species_bc = makeConstantSpeciesBoundaryConditions(species.size(), 1.0);

  equation::FixedStepMultiSpeciesDriftDiffusionStepper transport(
      mesh_, species,

      1.0, transport_dt,

      std::move(potential_bc),

      std::move(species_bc),

      std::make_unique<linalg::CholmodSolver>());

  EXPECT_THROW(FixedStepPlasmaSimulation(density, reactions, transport,

                                         NoReactionEvaluator{},

                                         FixedStepClock(clock_dt, 10)),
               std::invalid_argument);
}

// ============================================================
// 7. Once simulation is finished, another timestep must be
//    rejected.
//
// This prevents accidentally advancing:
//
//     n^(N)
//
// while the clock still claims the configured interval has
// already ended.
// ============================================================

TEST_F(FixedStepPlasmaSimulationTest, RejectsAdvanceAfterSimulationFinished) {
  constexpr double dt = 0.01;

  physics::SpeciesSet species;

  auto _ = addElectronAndIon(species);

  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);

  physics::ReactionNetwork reactions(species);

  auto potential_bc = makeZeroPotentialBoundaryConditions();

  auto species_bc = makeConstantSpeciesBoundaryConditions(species.size(), 1.0);

  equation::FixedStepMultiSpeciesDriftDiffusionStepper transport(
      mesh_, species,

      1.0, dt,

      std::move(potential_bc),

      std::move(species_bc),

      std::make_unique<linalg::CholmodSolver>());

  FixedStepPlasmaSimulation simulation(density, reactions, transport,

                                       NoReactionEvaluator{},

                                       FixedStepClock(dt, 1));

  ASSERT_TRUE(simulation.advanceOneStep().success());

  ASSERT_TRUE(simulation.finished());

  EXPECT_THROW((void)simulation.advanceOneStep(), std::out_of_range);
}

TEST(AdaptiveTimeClockTest, LandsExactlyOnEndTime) {
  AdaptiveTimeClock clock(0.1);

  clock.advance(0.04);
  clock.advance(0.04);
  clock.advance(0.02);

  EXPECT_TRUE(clock.finished());
  EXPECT_EQ(clock.step(), 3u);
  EXPECT_NEAR(clock.time(), 0.1, 1e-15);
  EXPECT_NEAR(clock.remainingTime(), 0.0, 1e-15);
}

TEST_F(AdaptiveStepPlasmaSimulationTest,
       RunSelectsVariableStepsAndReachesEndTime) {
  physics::SpeciesSet species;
  const auto ids = addElectronAndIon(species);
  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);
  physics::ReactionNetwork reactions(species);

  equation::AdaptiveStepMultiSpeciesDriftDiffusionStepper transport(
      mesh_, species, 1.0, makeZeroPotentialBoundaryConditions(),
      makeConstantSpeciesBoundaryConditions(species.size(), 1.0),
      std::make_unique<linalg::CholmodSolver>(),
      {.safety = 0.9, .min_dt = 1e-8, .max_dt = 0.04, .max_growth = 2.0});

  AdaptiveStepPlasmaSimulation simulation(density, reactions, transport,
                                          NoReactionEvaluator{},
                                          AdaptiveTimeClock(0.1));

  simulation.run();

  EXPECT_TRUE(simulation.finished());
  EXPECT_EQ(simulation.step(), 3u);
  EXPECT_NEAR(simulation.time(), 0.1, 1e-15);
  EXPECT_NEAR(simulation.lastTimeStep(), 0.02, 1e-14);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {
    EXPECT_NEAR(density[ids.electron][cell], 1.0, 1e-12);
    EXPECT_NEAR(density[ids.ion][cell], 1.0, 1e-12);
  }
}

}  // namespace pemu::simulation::test
