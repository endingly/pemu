#include <gtest/gtest.h>
#include <pemu/boundary/boundary_condition_set.hpp>
#include <pemu/equation/adaptive_step_multi_species_drift_diffusion_stepper.hpp>
#include <pemu/equation/fixed_step_multi_species_drift_diffusion_stepper.hpp>
#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>
#include <pemu/field/plasma_field_metadata.hpp>
#include <pemu/field/quantity_io.hpp>
#include <pemu/linalg/cholmod_solver.hpp>
#include <pemu/linalg/i_solver.hpp>
#include <pemu/mesh/moab_mesh.hpp>
#include <pemu/physics/reaction.hpp>
#include <pemu/physics/species.hpp>
#include <pemu/simulation/adaptive_step_plasma_simulation.hpp>
#include <pemu/simulation/adaptive_time_clock.hpp>
#include <pemu/simulation/fixed_step_clock.hpp>
#include <pemu/simulation/fixed_step_plasma_simulation.hpp>

#include <mp-units/systems/si.h>

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
std::filesystem::path poisson64x64MeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "poisson_64x64.msh";
}

[[nodiscard]]
double integratedDensity(const mesh::IMesh& mesh,
                         const field::CellField<double>& density) {
  double total = 0.0;

  for (mesh::CellId cell = 0; cell < mesh.numCells(); ++cell) {
    total += density[cell] * mesh.cellVolume(cell);
  }

  return total;
}

[[nodiscard]]
double densityCentroidX(const mesh::IMesh& mesh,
                        const field::CellField<double>& density) {
  double total_density = 0.0;
  double first_moment = 0.0;

  for (mesh::CellId cell = 0; cell < mesh.numCells(); ++cell) {
    const double weighted_density = density[cell] * mesh.cellVolume(cell);

    total_density += weighted_density;

    first_moment += weighted_density * mesh.cellCenter(cell).x;
  }

  if (total_density <= 0.0) {
    throw std::runtime_error(
        "density centroid requires positive total density");
  }

  return first_moment / total_density;
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
boundary::BoundaryConditionSet makeParallelPlatePotentialBoundaryConditions(
    double left_voltage, double right_voltage) {
  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(mesh::BoundaryId{1}, left_voltage);

  bc.setDirichlet(mesh::BoundaryId{2}, right_voltage);

  // The remaining two sides are electrically insulating.
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

class AdaptiveStepPlasmaSimulation64x64Test : public ::testing::Test {
 protected:
  AdaptiveStepPlasmaSimulation64x64Test() : mesh_(poisson64x64MeshPath()) {}

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

TEST_F(AdaptiveStepPlasmaSimulationTest,
       AdvanceOneStepEvaluatesReactionAndUpdatesSpecies) {
  physics::SpeciesSet species;
  const auto ids = addElectronAndIon(species);
  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);

  physics::ReactionNetwork reactions(species);
  const auto ionization = reactions.addReaction(
      {.name = "electron impact ionization",
       .stoichiometry = {{ids.electron, +1.0}, {ids.ion, +1.0}}});

  equation::AdaptiveStepMultiSpeciesDriftDiffusionStepper transport(
      mesh_, species, 1.0, makeZeroPotentialBoundaryConditions(),
      makeConstantSpeciesBoundaryConditions(species.size(), 1.0),
      std::make_unique<linalg::CholmodSolver>(),
      {.safety = 0.9, .min_dt = 1e-8, .max_dt = 0.04, .max_growth = 2.0});

  ElectronImpactIonizationEvaluator evaluator{.electron = ids.electron,
                                              .ionization = ionization,
                                              .neutral_density = 4.0,
                                              .rate_coefficient = 0.5};

  AdaptiveStepPlasmaSimulation simulation(density, reactions, transport,
                                          evaluator, AdaptiveTimeClock(0.1));

  ASSERT_TRUE(simulation.advanceOneStep().success());
  ASSERT_TRUE(simulation.hasLastTimeStepProposal());

  EXPECT_NEAR(simulation.time(), 0.04, 1e-14);
  EXPECT_NEAR(simulation.lastTimeStep(), 0.04, 1e-14);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {
    EXPECT_NEAR(simulation.reactionRates()[ionization][cell], 2.0, 1e-12);
    EXPECT_NEAR(simulation.source()[ids.electron][cell], 2.0, 1e-12);
    EXPECT_NEAR(simulation.source()[ids.ion][cell], 2.0, 1e-12);
    EXPECT_NEAR(density[ids.electron][cell], 1.08, 1e-12);
    EXPECT_NEAR(density[ids.ion][cell], 1.08, 1e-12);
  }
}

TEST_F(AdaptiveStepPlasmaSimulationTest,
       ReactionRateIsReevaluatedFromUpdatedStateEveryStep) {
  physics::SpeciesSet species;
  const auto ids = addElectronAndIon(species);
  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);

  physics::ReactionNetwork reactions(species);
  const auto ionization = reactions.addReaction(
      {.name = "ionization",
       .stoichiometry = {{ids.electron, +1.0}, {ids.ion, +1.0}}});

  equation::AdaptiveStepMultiSpeciesDriftDiffusionStepper transport(
      mesh_, species, 1.0, makeZeroPotentialBoundaryConditions(),
      makeConstantSpeciesBoundaryConditions(species.size(), 1.0),
      std::make_unique<linalg::CholmodSolver>(),
      {.safety = 0.9, .min_dt = 1e-8, .max_dt = 0.04, .max_growth = 2.0});

  std::vector<double> observations;
  RecordingIonizationEvaluator evaluator{
      .electron = ids.electron,
      .ionization = ionization,
      .neutral_density = 4.0,
      .rate_coefficient = 0.5,
      .observed_electron_density = &observations};

  AdaptiveStepPlasmaSimulation simulation(density, reactions, transport,
                                          evaluator, AdaptiveTimeClock(0.08));

  simulation.run();

  ASSERT_EQ(observations.size(), 2u);
  EXPECT_NEAR(observations[0], 1.0, 1e-12);
  EXPECT_NEAR(observations[1], 1.08, 1e-12);
}

TEST_F(AdaptiveStepPlasmaSimulationTest,
       ReactionEvaluatorSeesCurrentElectricField) {
  physics::SpeciesSet species;
  const auto ids = addElectronAndIon(species);
  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);

  physics::ReactionNetwork reactions(species);
  const auto reaction = reactions.addReaction(
      {.name = "field dependent test reaction",
       .stoichiometry = {{ids.electron, +1.0}, {ids.ion, +1.0}}});

  equation::AdaptiveStepMultiSpeciesDriftDiffusionStepper transport(
      mesh_, species, 1.0, makeLinearXPotentialBoundaryConditions(mesh_),
      makeConstantSpeciesBoundaryConditions(species.size(), 1.0),
      std::make_unique<linalg::CholmodSolver>(),
      {.safety = 0.9, .min_dt = 1e-8, .max_dt = 1e-3, .max_growth = 2.0});

  AdaptiveStepPlasmaSimulation simulation(
      density, reactions, transport,
      ElectricFieldAwareEvaluator{.reaction = reaction},
      AdaptiveTimeClock(1e-3));

  ASSERT_TRUE(simulation.advanceOneStep().success());
  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {
    EXPECT_NEAR(simulation.reactionRates()[reaction][cell], 1.0, 1e-12);
  }
}

TEST_F(AdaptiveStepPlasmaSimulationTest,
       ReactionSinkLimitsTimeStepAndPreservesNonNegativeDensity) {
  physics::SpeciesSet species;
  const auto ids = addElectronAndIon(species);
  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);

  physics::ReactionNetwork reactions(species);
  const auto loss = reactions.addReaction(
      {.name = "pair loss",
       .stoichiometry = {{ids.electron, -1.0}, {ids.ion, -1.0}}});

  equation::AdaptiveStepMultiSpeciesDriftDiffusionStepper transport(
      mesh_, species, 1.0, makeZeroPotentialBoundaryConditions(),
      makeConstantSpeciesBoundaryConditions(species.size(), 1.0),
      std::make_unique<linalg::CholmodSolver>(),
      {.safety = 0.9, .min_dt = 1e-8, .max_dt = 0.5, .max_growth = 2.0});

  ElectronImpactIonizationEvaluator evaluator{.electron = ids.electron,
                                              .ionization = loss,
                                              .neutral_density = 1.0,
                                              .rate_coefficient = 10.0};

  AdaptiveStepPlasmaSimulation simulation(density, reactions, transport,
                                          evaluator, AdaptiveTimeClock(1.0));

  ASSERT_TRUE(simulation.advanceOneStep().success());
  const auto& proposal = simulation.lastTimeStepProposal();

  EXPECT_LT(proposal.dt, 0.1);
  EXPECT_GT(proposal.dt, 0.0);
  EXPECT_LT(proposal.positivity_limit, 0.1);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {
    EXPECT_NEAR(density[ids.electron][cell], 1.0 - 10.0 * proposal.dt, 1e-12);
    EXPECT_NEAR(density[ids.ion][cell], 1.0 - 10.0 * proposal.dt, 1e-12);
    EXPECT_GE(density[ids.electron][cell], 0.0);
    EXPECT_GE(density[ids.ion][cell], 0.0);
  }
}

TEST_F(AdaptiveStepPlasmaSimulationTest,
       SolvesPoissonOncePerAdaptiveStepAndReusesFactorization) {
  physics::SpeciesSet species;
  auto _ = addElectronAndIon(species);
  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);
  physics::ReactionNetwork reactions(species);

  auto backend = std::make_unique<CountingSolver>();
  auto* counting_solver = backend.get();
  equation::AdaptiveStepMultiSpeciesDriftDiffusionStepper transport(
      mesh_, species, 1.0, makeZeroPotentialBoundaryConditions(),
      makeConstantSpeciesBoundaryConditions(species.size(), 1.0),
      std::move(backend),
      {.safety = 0.9, .min_dt = 1e-8, .max_dt = 0.04, .max_growth = 2.0});

  AdaptiveStepPlasmaSimulation simulation(density, reactions, transport,
                                          NoReactionEvaluator{},
                                          AdaptiveTimeClock(0.1));
  simulation.run();

  EXPECT_EQ(simulation.step(), 3u);
  EXPECT_EQ(counting_solver->analyze_count, 1);
  EXPECT_EQ(counting_solver->factorize_count, 1);
  EXPECT_EQ(counting_solver->solve_count, 3);
}

TEST_F(AdaptiveStepPlasmaSimulationTest,
       RejectsAdvanceAfterSimulationFinished) {
  physics::SpeciesSet species;
  auto _ = addElectronAndIon(species);
  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);
  physics::ReactionNetwork reactions(species);

  equation::AdaptiveStepMultiSpeciesDriftDiffusionStepper transport(
      mesh_, species, 1.0, makeZeroPotentialBoundaryConditions(),
      makeConstantSpeciesBoundaryConditions(species.size(), 1.0),
      std::make_unique<linalg::CholmodSolver>(),
      {.safety = 0.9, .min_dt = 1e-8, .max_dt = 0.04, .max_growth = 2.0});

  AdaptiveStepPlasmaSimulation simulation(density, reactions, transport,
                                          NoReactionEvaluator{},
                                          AdaptiveTimeClock(0.04));
  ASSERT_TRUE(simulation.advanceOneStep().success());
  ASSERT_TRUE(simulation.finished());
  EXPECT_THROW((void)simulation.advanceOneStep(), std::out_of_range);
}

TEST_F(AdaptiveStepPlasmaSimulation64x64Test,
       SolvesUniformElectronImpactIonizationAcrossMultipleSteps) {
  constexpr double initial_density = 1.0;
  constexpr double neutral_density = 4.0;
  constexpr double rate_coefficient = 0.5;
  constexpr double reaction_rate = neutral_density * rate_coefficient;
  constexpr double time_step = 0.01;
  constexpr std::size_t total_steps = 5;
  constexpr double end_time = time_step * static_cast<double>(total_steps);

  ASSERT_EQ(mesh_.numCells(), 64u * 64u);

  physics::SpeciesSet species;
  // Keep a strictly positive diffusivity, as required by the SG transport
  // operator, while making its effect negligible for this uniform reaction
  // regression. This leaves a multi-step discrete reference solution.
  const auto ids = addElectronAndIon(species, 1.0, 0.5, 1e-16, 1e-16);
  physics::SpeciesCellFields density(mesh_, species.size(), initial_density);

  physics::ReactionNetwork reactions(species);
  const auto ionization = reactions.addReaction(
      {.name = "uniform electron-impact ionization",
       .stoichiometry = {{ids.electron, +1.0}, {ids.ion, +1.0}}});

  equation::AdaptiveStepMultiSpeciesDriftDiffusionStepper transport(
      mesh_, species, 1.0, makeZeroPotentialBoundaryConditions(),
      makeConstantSpeciesBoundaryConditions(species.size(), initial_density),
      std::make_unique<linalg::CholmodSolver>(),
      {.safety = 0.9, .min_dt = 1e-10, .max_dt = time_step, .max_growth = 2.0});

  ElectronImpactIonizationEvaluator evaluator{
      .electron = ids.electron,
      .ionization = ionization,
      .neutral_density = neutral_density,
      .rate_coefficient = rate_coefficient};

  AdaptiveStepPlasmaSimulation simulation(
      density, reactions, transport, evaluator, AdaptiveTimeClock(end_time));

  simulation.run();
  ASSERT_TRUE(simulation.finished());
  ASSERT_TRUE(simulation.hasLastTimeStepProposal());
  EXPECT_EQ(simulation.step(), total_steps);
  EXPECT_NEAR(simulation.time(), end_time, 1e-15);
  EXPECT_NEAR(simulation.lastTimeStep(), time_step, 1e-15);

  const double growth_per_step = 1.0 + time_step * reaction_rate;
  const double expected_density =
      initial_density * std::pow(growth_per_step, total_steps);
  const double final_reaction_rate = reaction_rate * initial_density *
                                     std::pow(growth_per_step, total_steps - 1);

  EXPECT_NEAR(integratedDensity(mesh_, density[ids.electron]), expected_density,
              1e-12);
  EXPECT_NEAR(integratedDensity(mesh_, density[ids.ion]), expected_density,
              1e-12);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {
    EXPECT_TRUE(std::isfinite(density[ids.electron][cell]));
    EXPECT_TRUE(std::isfinite(density[ids.ion][cell]));
    EXPECT_GE(density[ids.electron][cell], 0.0);
    EXPECT_GE(density[ids.ion][cell], 0.0);
    EXPECT_NEAR(simulation.reactionRates()[ionization][cell],
                final_reaction_rate, 1e-12);
    EXPECT_NEAR(simulation.source()[ids.electron][cell], final_reaction_rate,
                1e-12);
    EXPECT_NEAR(simulation.source()[ids.ion][cell], final_reaction_rate, 1e-12);
    EXPECT_NEAR(density[ids.electron][cell], expected_density, 1e-12);
    EXPECT_NEAR(density[ids.ion][cell], expected_density, 1e-12);
    EXPECT_NEAR(simulation.chargeDensity()[cell], 0.0, 1e-12);
    EXPECT_NEAR(simulation.potential()[cell], 0.0, 1e-12);
  }

  for (const double electric_field : simulation.electricFieldNormal()) {
    EXPECT_TRUE(std::isfinite(electric_field));
    EXPECT_NEAR(electric_field, 0.0, 1e-12);
  }
}

TEST_F(AdaptiveStepPlasmaSimulation64x64Test,
       ParallelPlate400VDrivesOppositeDriftAndIonizationInCentimeterMesh) {
  using namespace mp_units;
  using namespace mp_units::si::unit_symbols;

  constexpr auto number_density_unit = one / cubic(cm);
  constexpr auto mobility_unit = square(cm) / (V * s);
  constexpr auto diffusivity_unit = square(cm) / s;
  constexpr auto ionization_coefficient_unit = cubic(cm) / s;

  // Strong quantities are converted to raw doubles only at solver boundaries.
  constexpr auto domain_length = 1.0 * cm;
  constexpr auto left_potential = 0.0 * V;
  constexpr auto right_potential = 400.0 * V;
  constexpr auto expected_electric_field_x =
      -(right_potential - left_potential) / domain_length;
  constexpr auto elementary_charge = 1.602176634e-19 * C;
  constexpr auto vacuum_permittivity = 8.8541878128e-14 * F / cm;
  constexpr auto initial_density = 1.0e6 * number_density_unit;
  constexpr auto electron_mobility = 1.0e3 * mobility_unit;
  constexpr auto ion_mobility = 1.5 * mobility_unit;
  constexpr auto electron_diffusivity = 1.0e2 * diffusivity_unit;
  constexpr auto ion_diffusivity = 4.0e-2 * diffusivity_unit;
  constexpr auto neutral_density = 2.5e19 * number_density_unit;
  constexpr auto ionization_rate_coefficient =
      1.0e-13 * ionization_coefficient_unit;
  constexpr auto end_time = 2.0e-7 * s;
  constexpr auto maximum_time_step = 1.0e-6 * s;

  constexpr double domain_length_cm = domain_length.numerical_value_in(cm);
  constexpr double left_voltage = left_potential.numerical_value_in(V);
  constexpr double right_voltage = right_potential.numerical_value_in(V);
  constexpr double elementary_charge_coulomb =
      elementary_charge.numerical_value_in(C);
  constexpr double initial_density_cm3 =
      initial_density.numerical_value_in(number_density_unit);
  constexpr double end_time_s = end_time.numerical_value_in(s);
  constexpr double max_time_step_s = maximum_time_step.numerical_value_in(s);

  const auto field_metadata = field::centimetrePlasmaFieldMetadata();

  ASSERT_EQ(mesh_.numCells(), 64u * 64u);

  physics::SpeciesSet species;
  const auto electron = species.add(
      {.name = "e",
       .charge = -elementary_charge_coulomb,
       .mobility = electron_mobility.numerical_value_in(mobility_unit),
       .diffusivity = electron_diffusivity.numerical_value_in(diffusivity_unit),
       .transport_model = physics::SpeciesTransportModel::DriftDiffusion});
  const auto ion = species.add(
      {.name = "Ar+",
       .charge = elementary_charge_coulomb,
       .mobility = ion_mobility.numerical_value_in(mobility_unit),
       .diffusivity = ion_diffusivity.numerical_value_in(diffusivity_unit),
       .transport_model = physics::SpeciesTransportModel::DriftDiffusion});
  physics::SpeciesCellFields density(mesh_, species.size(), initial_density_cm3,
                                     field_metadata.number_density);

  physics::ReactionNetwork reactions(species);
  const auto ionization =
      reactions.addReaction({.name = "electron-impact ionization",
                             .stoichiometry = {{electron, +1.0}, {ion, +1.0}}});

  equation::AdaptiveStepMultiSpeciesDriftDiffusionStepper transport(
      mesh_, species, vacuum_permittivity.numerical_value_in(F / cm),
      makeParallelPlatePotentialBoundaryConditions(left_voltage, right_voltage),
      makeConstantSpeciesBoundaryConditions(species.size(),
                                            initial_density_cm3),
      std::make_unique<linalg::CholmodSolver>(),
      {.safety = 0.8,
       .min_dt = 1.0e-12,
       .max_dt = max_time_step_s,
       .max_growth = 1.5},
      field_metadata);

  // Initial neutrality makes the parallel-plate potential an exact solution.
  ASSERT_TRUE(transport.prepareElectrostatics(density).success());
  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {
    const double expected_potential =
        left_voltage + (right_voltage - left_voltage) *
                           mesh_.cellCenter(cell).x / domain_length_cm;
    EXPECT_NEAR(transport.potential()[cell], expected_potential, 1.0e-8);
  }

  const auto left_face = findBoundaryFace(mesh_, mesh::BoundaryId{1});
  const auto right_face = findBoundaryFace(mesh_, mesh::BoundaryId{2});
  EXPECT_NEAR(transport.electricFieldNormal()[left_face],
              -expected_electric_field_x.numerical_value_in(V / cm), 1.0e-8);
  EXPECT_NEAR(transport.electricFieldNormal()[right_face],
              expected_electric_field_x.numerical_value_in(V / cm), 1.0e-8);

  ASSERT_TRUE(density[electron].metadata().hasPhysicalQuantity());
  ASSERT_TRUE(transport.potential().metadata().hasPhysicalQuantity());
  ASSERT_TRUE(transport.electricFieldNormal().metadata().hasPhysicalQuantity());
  EXPECT_TRUE(density[electron].metadata().physical_quantity->represents(
      pemu::unit::plasma_quantity::particle_number_density,
      number_density_unit));
  EXPECT_TRUE(transport.potential().metadata().physical_quantity->represents(
      isq::electric_potential, V));
  EXPECT_TRUE(
      transport.electricFieldNormal().metadata().physical_quantity->represents(
          pemu::unit::plasma_quantity::normal_electric_field_strength, V / cm));

  const double initial_integrated_density =
      integratedDensity(mesh_, density[electron]);
  const double initial_centroid_x = densityCentroidX(mesh_, density[electron]);

  ElectronImpactIonizationEvaluator evaluator{
      .electron = electron,
      .ionization = ionization,
      .neutral_density =
          neutral_density.numerical_value_in(number_density_unit),
      .rate_coefficient = ionization_rate_coefficient.numerical_value_in(
          ionization_coefficient_unit)};
  AdaptiveStepPlasmaSimulation simulation(
      density, reactions, transport, evaluator, AdaptiveTimeClock(end_time_s));

  simulation.run();

  ASSERT_TRUE(simulation.finished());
  ASSERT_TRUE(simulation.hasLastTimeStepProposal());
  ASSERT_TRUE(
      simulation.reactionRates()[ionization].metadata().hasPhysicalQuantity());
  ASSERT_TRUE(simulation.source()[electron].metadata().hasPhysicalQuantity());
  EXPECT_TRUE(simulation.reactionRates()[ionization]
                  .metadata()
                  .physical_quantity->represents(
                      pemu::unit::plasma_quantity::reaction_rate_density,
                      one / (cubic(cm) * s)));
  EXPECT_TRUE(
      simulation.source()[electron].metadata().physical_quantity->represents(
          pemu::unit::plasma_quantity::particle_number_density_rate,
          one / (cubic(cm) * s)));
  EXPECT_NEAR(simulation.time(), end_time_s, 1.0e-18);
  EXPECT_GT(simulation.step(), 1u);
  EXPECT_LT(simulation.lastTimeStep(), max_time_step_s);

  EXPECT_GT(integratedDensity(mesh_, density[electron]),
            initial_integrated_density);
  EXPECT_GT(integratedDensity(mesh_, density[ion]), initial_integrated_density);
  const double electron_centroid_x = densityCentroidX(mesh_, density[electron]);
  const double ion_centroid_x = densityCentroidX(mesh_, density[ion]);
  EXPECT_GT(electron_centroid_x, initial_centroid_x);
  // Ionization is stronger where the electron density has drifted. Therefore
  // the ion centroid need not move left of its initial value, but the faster
  // electron population must still lie to its right.
  EXPECT_GT(electron_centroid_x, ion_centroid_x);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {
    EXPECT_TRUE(std::isfinite(density[electron][cell]));
    EXPECT_TRUE(std::isfinite(density[ion][cell]));
    EXPECT_GE(density[electron][cell], 0.0);
    EXPECT_GE(density[ion][cell], 0.0);
    EXPECT_GT(simulation.reactionRates()[ionization][cell], 0.0);
    EXPECT_GT(simulation.source()[electron][cell], 0.0);
    EXPECT_GT(simulation.source()[ion][cell], 0.0);
  }
}

}  // namespace pemu::simulation::test
