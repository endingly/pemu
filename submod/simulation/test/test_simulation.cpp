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
#include <pemu/output/checkpoint/vtkhdf.hpp>
#include <pemu/output/dump/i_writer.hpp>
#include <pemu/output/dump/vtkhdf_writer.hpp>
#include <pemu/physics/reaction.hpp>
#include <pemu/physics/species.hpp>
#include <pemu/simulation/adaptive_time_clock.hpp>
#include <pemu/simulation/fixed_step_clock.hpp>
#include <pemu/simulation/workflow.hpp>
#include <pemu/trace/ostream_trace_sink.hpp>
#include <pemu/trace/split_trace_sink.hpp>
#include <pemu/trace/trace.hpp>

#include <mp-units/systems/si.h>

#include <vtkFieldData.h>
#include <vtkHDFReader.h>
#include <vtkNew.h>
#include <vtkUnstructuredGrid.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
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
  double* observed_maximum{};

  void operator()(const physics::SpeciesCellFields&,

                  const field::CellField<double>&,

                  const field::FaceField<double>& electric_field,

                  physics::ReactionRateFields& reaction_rates) const {
    double maximum = 0.0;

    for (const double value : electric_field) {

      maximum = std::max(maximum, std::abs(value));
    }

    if (observed_maximum != nullptr) {
      *observed_maximum = maximum;
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

struct RecordingTraceSink {
  std::vector<std::string>* names{};
  std::vector<std::string>* categories{};
  std::vector<pemu::trace::DiagDomain>* domains{};
  std::vector<pemu::trace::EventKind>* kinds{};
  std::vector<double>* selected_time_steps{};
  std::vector<double>* completed_times{};

  void operator()(const pemu::trace::TraceEvent& event) const noexcept {
    names->emplace_back(event.name);
    categories->emplace_back(event.category);
    domains->push_back(event.domain);
    kinds->push_back(event.kind);

    for (const auto& attribute : event.attributes) {
      const auto* value = std::get_if<double>(&attribute.value);

      if (value == nullptr) {
        continue;
      }

      if (event.name == "timestep.selected" && attribute.name == "dt") {
        selected_time_steps->push_back(*value);
      }

      if (event.name == "step.completed" && attribute.name == "time") {
        completed_times->push_back(*value);
      }
    }
  }

  void flush() const noexcept {}
};

static_assert(pemu::trace::TraceSink<RecordingTraceSink>);

struct FlushCountingTraceSink {
  std::size_t* flush_count{};

  void operator()(const pemu::trace::TraceEvent&) const noexcept {}
  void flush() const noexcept { ++*flush_count; }
};

static_assert(pemu::trace::TraceSink<FlushCountingTraceSink>);

struct CapturedFieldOutput {
  output::OutputStamp stamp;
  std::filesystem::path path;
  std::vector<std::string> cell_names;
  std::vector<std::string> face_names;
};

class CapturingFieldOutputSeries final : public output::dump::ISeries {
 public:
  /** @brief Creates a capture-only series backed by caller-owned records. */
  CapturingFieldOutputSeries(std::vector<CapturedFieldOutput>& requests,
                             std::filesystem::path path,
                             std::size_t& finish_count)
      : requests_(&requests),
        path_(std::move(path)),
        finish_count_(&finish_count) {}

  /** @copydoc output::dump::ISeries::append */
  [[nodiscard]]
  output::OutputRecord append(const output::dump::Request& request) override {
    CapturedFieldOutput captured{.stamp = request.stamp, .path = request.path};
    for (const auto& selection : request.cell_field_selections) {
      if (selection.field == nullptr) {
        throw std::invalid_argument("captured cell field must not be null");
      }
      captured.cell_names.push_back(selection.name.empty()
                                        ? selection.field->metadata().name
                                        : selection.name);
    }
    for (const auto& selection : request.face_fields) {
      if (selection.field == nullptr) {
        throw std::invalid_argument("captured face field must not be null");
      }
      captured.face_names.push_back(selection.name.empty()
                                        ? selection.field->metadata().name
                                        : selection.name);
    }
    requests_->push_back(std::move(captured));
    return {.path = request.path, .stamp = request.stamp};
  }

  /** @copydoc output::dump::ISeries::finish */
  [[nodiscard]] output::OutputRecord finish() override {
    if (requests_->empty()) {
      throw std::logic_error("cannot finish an empty captured series");
    }
    ++*finish_count_;
    return {.path = path_, .stamp = requests_->back().stamp};
  }

 private:
  std::vector<CapturedFieldOutput>* requests_{};
  std::filesystem::path path_;
  std::size_t* finish_count_{};
};

class CapturingFieldOutputWriter final : public output::dump::IWriter {
 public:
  /** @copydoc output::dump::IWriter::openSeries */
  [[nodiscard]] std::unique_ptr<output::dump::ISeries> openSeries(
      const output::dump::SeriesRequest& request) const override {
    return std::make_unique<CapturingFieldOutputSeries>(requests, request.path,
                                                        finish_count);
  }

  mutable std::vector<CapturedFieldOutput> requests;
  mutable std::size_t finish_count{};
};

class FailingFinishSeries final : public output::dump::ISeries {
 public:
  [[nodiscard]] output::OutputRecord append(
      const output::dump::Request& request) override {
    last_record_ = {.path = request.path, .stamp = request.stamp};
    return last_record_;
  }

  [[nodiscard]] output::OutputRecord finish() override {
    throw std::runtime_error("injected field output finish failure");
  }

 private:
  output::OutputRecord last_record_;
};

class FailingFinishWriter final : public output::dump::IWriter {
 public:
  [[nodiscard]] std::unique_ptr<output::dump::ISeries> openSeries(
      const output::dump::SeriesRequest&) const override {
    return std::make_unique<FailingFinishSeries>();
  }
};

class CountingCheckpointReader final : public output::checkpoint::IReader {
 public:
  [[nodiscard]] output::checkpoint::Manifest inspect(
      const std::filesystem::path& path) const override {
    ++inspect_count;
    return reader_.inspect(path);
  }

  [[nodiscard]] output::OutputRecord restore(
      const std::filesystem::path& path,
      const output::checkpoint::RestoreRequest& request) const override {
    ++restore_count;
    return reader_.restore(path, request);
  }

  mutable std::size_t inspect_count{};
  mutable std::size_t restore_count{};

 private:
  output::checkpoint::VtkHdfReader reader_;
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

  linalg::SolverResult analyzePattern(const linalg::SparseMatrix&) override {
    ++analyze_count;

    analyzed_ = true;

    return {.status = linalg::SolverStatus::Success};
  }

  linalg::SolverResult factorize(const linalg::SparseMatrix&) override {
    ++factorize_count;

    factorized_ = true;

    return {.status = linalg::SolverStatus::Success};
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

class DiagnosticFailingSolver final : public linalg::ISolver {
 public:
  linalg::SolverResult analyzePattern(const linalg::SparseMatrix&) override {
    analyzed_ = true;
    return {.status = linalg::SolverStatus::Success};
  }

  linalg::SolverResult factorize(const linalg::SparseMatrix&) override {
    factorized_ = true;
    return {.status = linalg::SolverStatus::Success};
  }

  linalg::SolverResult solve(linalg::ConstVectorRef,
                             linalg::VectorRef) override {
    return {
        .status = linalg::SolverStatus::SolveFailed,
        .residual_norm = 2.5,
        .relative_residual = 0.25,
        .diagnostic = trace::makeDiagnosticEvent(
            trace::DiagDomain::linalg, "test_backend", "solve.failed",
            "injected linear solve failure"),
    };
  }

  void reset() override {
    analyzed_ = false;
    factorized_ = false;
  }

  bool isAnalyzed() const noexcept override { return analyzed_; }

  bool isFactorized() const noexcept override { return factorized_; }

 private:
  bool analyzed_{false};
  bool factorized_{false};
};

// Runs a real CHOLMOD solve before injecting a failure.  This makes the
// 64x64 integration test exercise one complete physical timestep before it
// verifies that the lower-level linalg diagnostic survives simulation's
// failure handling.
class FailAfterSuccessfulCholmodSolves final : public linalg::ISolver {
 public:
  explicit FailAfterSuccessfulCholmodSolves(
      std::size_t successful_solves_before_failure)
      : successful_solves_before_failure_(successful_solves_before_failure) {}

  linalg::SolverResult analyzePattern(
      const linalg::SparseMatrix& matrix) override {
    return backend_.analyzePattern(matrix);
  }

  linalg::SolverResult factorize(const linalg::SparseMatrix& matrix) override {
    return backend_.factorize(matrix);
  }

  linalg::SolverResult solve(linalg::ConstVectorRef right_hand_side,
                             linalg::VectorRef solution) override {
    if (successful_solve_count_ >= successful_solves_before_failure_) {
      return {
          .status = linalg::SolverStatus::SolveFailed,
          .residual_norm = 2.5,
          .relative_residual = 0.25,
          .diagnostic = trace::makeDiagnosticEvent(
              trace::DiagDomain::linalg, "test_backend", "solve.failed",
              "injected failure after successful CHOLMOD solve"),
      };
    }

    auto result = backend_.solve(right_hand_side, solution);
    if (result.success()) {
      ++successful_solve_count_;
    }
    return result;
  }

  void reset() override {
    backend_.reset();
    successful_solve_count_ = 0;
  }

  [[nodiscard]] bool isAnalyzed() const noexcept override {
    return backend_.isAnalyzed();
  }

  [[nodiscard]] bool isFactorized() const noexcept override {
    return backend_.isFactorized();
  }

 private:
  linalg::CholmodSolver backend_;
  std::size_t successful_solves_before_failure_{};
  std::size_t successful_solve_count_{};
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

class PlasmaSimulation64x64CheckpointTest : public ::testing::Test {
 protected:
  PlasmaSimulation64x64CheckpointTest() : mesh_(poisson64x64MeshPath()) {}

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

TEST(FixedStepClockTest, RestoresInitialStepAndDerivedTime) {
  const FixedStepClock clock(0.125, 8, 3);

  EXPECT_EQ(clock.step(), 3u);
  EXPECT_DOUBLE_EQ(clock.time(), 0.375);
  EXPECT_FALSE(clock.finished());
  EXPECT_THROW(auto _ = FixedStepClock(0.125, 8, 9), std::invalid_argument);
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
       AdvanceEvaluatesReactionAndUpdatesSpecies) {
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

  const auto result = simulation.advance();

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

  ASSERT_TRUE(simulation.advance().success());

  ASSERT_TRUE(simulation.advance().success());

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

  double evaluator_maximum{};
  ElectricFieldAwareEvaluator evaluator{.reaction = reaction,
                                        .observed_maximum = &evaluator_maximum};

  FixedStepPlasmaSimulation simulation(density, reactions, transport, evaluator,

                                       FixedStepClock(dt, 1));

  const auto result = simulation.advance();

  ASSERT_TRUE(result.success());
  EXPECT_NEAR(evaluator_maximum, 1.0, 1e-12);

  // --------------------------------------------------------
  // phi=x -> E=(-1,0)
  //
  // Therefore the max absolute face-normal field is 1.
  // --------------------------------------------------------

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

  EXPECT_EQ(simulation.state(), SimulationState::completed);
}

TEST_F(FixedStepPlasmaSimulationTest,
       StateMachineSupportsStartPauseResumeAndStop) {
  constexpr double dt = 0.01;
  physics::SpeciesSet species;
  (void)addElectronAndIon(species);
  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);
  physics::ReactionNetwork reactions(species);
  equation::FixedStepMultiSpeciesDriftDiffusionStepper transport(
      mesh_, species, 1.0, dt, makeZeroPotentialBoundaryConditions(),
      makeConstantSpeciesBoundaryConditions(species.size(), 1.0),
      std::make_unique<linalg::CholmodSolver>());
  std::size_t flush_count{};
  FixedStepPlasmaSimulation simulation(
      density, reactions, transport, NoReactionEvaluator{},
      FixedStepClock(dt, 4), FlushCountingTraceSink{&flush_count});

  EXPECT_EQ(simulation.state(), SimulationState::ready);
  EXPECT_THROW(simulation.pause(), std::logic_error);

  simulation.start();
  EXPECT_EQ(simulation.state(), SimulationState::running);
  ASSERT_TRUE(simulation.advance().success());
  EXPECT_EQ(simulation.step(), 1u);

  simulation.pause();
  EXPECT_EQ(simulation.state(), SimulationState::paused);
  EXPECT_EQ(flush_count, 1u);
  EXPECT_THROW((void)simulation.advance(), std::logic_error);

  simulation.start();
  ASSERT_TRUE(simulation.advance().success());
  EXPECT_EQ(simulation.step(), 2u);
  EXPECT_EQ(flush_count, 2u);
  simulation.stop();
  EXPECT_EQ(simulation.state(), SimulationState::stopped);
  EXPECT_EQ(flush_count, 3u);
  EXPECT_THROW(simulation.run(), std::logic_error);
  EXPECT_NO_THROW(simulation.stop());
}

TEST_F(FixedStepPlasmaSimulationTest, StopFailureTransitionsWorkflowToFailed) {
  constexpr double dt = 0.01;
  physics::SpeciesSet species;
  (void)addElectronAndIon(species);
  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);
  physics::ReactionNetwork reactions(species);
  equation::FixedStepMultiSpeciesDriftDiffusionStepper transport(
      mesh_, species, 1.0, dt, makeZeroPotentialBoundaryConditions(),
      makeConstantSpeciesBoundaryConditions(species.size(), 1.0),
      std::make_unique<linalg::CholmodSolver>());
  FailingFinishWriter writer;
  FixedStepPlasmaSimulation simulation(
      density, reactions, transport, NoReactionEvaluator{},
      FixedStepClock(dt, 2), trace::NullTraceSink{}, {},
      {.writer = &writer,
       .directory = "simulation-output",
       .file_stem = "failing-finish",
       .write_initial = true});

  ASSERT_TRUE(simulation.advance().success());
  ASSERT_EQ(simulation.state(), SimulationState::running);
  EXPECT_THROW(simulation.stop(), std::runtime_error);
  EXPECT_EQ(simulation.state(), SimulationState::failed);
}

TEST_F(FixedStepPlasmaSimulationTest,
       EmitsOrderedTraceEventsAtEachPipelineStage) {
  constexpr double dt = 0.01;
  constexpr std::size_t total_steps = 2;
  physics::SpeciesSet species;
  (void)addElectronAndIon(species);
  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);
  physics::ReactionNetwork reactions(species);

  equation::FixedStepMultiSpeciesDriftDiffusionStepper transport(
      mesh_, species, 1.0, dt, makeZeroPotentialBoundaryConditions(),
      makeConstantSpeciesBoundaryConditions(species.size(), 1.0),
      std::make_unique<linalg::CholmodSolver>());

  std::vector<std::string> names;
  std::vector<std::string> categories;
  std::vector<pemu::trace::DiagDomain> domains;
  std::vector<pemu::trace::EventKind> kinds;
  std::vector<double> selected_time_steps;
  std::vector<double> completed_times;
  names.reserve(12);
  categories.reserve(12);
  domains.reserve(12);
  kinds.reserve(12);
  completed_times.reserve(2);
  const RecordingTraceSink sink{.names = &names,
                                .categories = &categories,
                                .domains = &domains,
                                .kinds = &kinds,
                                .selected_time_steps = &selected_time_steps,
                                .completed_times = &completed_times};

  FixedStepPlasmaSimulation simulation(density, reactions, transport,
                                       NoReactionEvaluator{},
                                       FixedStepClock(dt, total_steps), sink);

  simulation.run();

  const std::vector<std::string> expected_names{
      "run.started",
      "step.started",
      "electrostatics.completed",
      "reaction_rates.completed",
      "sources.completed",
      "step.completed",
      "step.started",
      "electrostatics.completed",
      "reaction_rates.completed",
      "sources.completed",
      "step.completed",
      "run.completed",
  };
  EXPECT_EQ(names, expected_names);
  EXPECT_TRUE(std::ranges::all_of(categories, [](const std::string& category) {
    return category == "fixed_step";
  }));
  EXPECT_TRUE(std::ranges::all_of(domains, [](pemu::trace::DiagDomain domain) {
    return domain == pemu::trace::DiagDomain::simulation;
  }));
  EXPECT_TRUE(std::ranges::all_of(kinds, [](pemu::trace::EventKind kind) {
    return kind == pemu::trace::EventKind::Trace;
  }));
  EXPECT_TRUE(selected_time_steps.empty());
  ASSERT_EQ(completed_times.size(), 2u);
  EXPECT_NEAR(completed_times[0], 0.01, 1e-14);
  EXPECT_NEAR(completed_times[1], 0.02, 1e-14);
}

TEST_F(FixedStepPlasmaSimulationTest,
       WritesInitialPeriodicAndFinalFieldSnapshots) {
  constexpr double dt = 0.01;
  physics::SpeciesSet species;
  (void)addElectronAndIon(species);
  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);
  physics::ReactionNetwork reactions(species);
  equation::FixedStepMultiSpeciesDriftDiffusionStepper transport(
      mesh_, species, 1.0, dt, makeZeroPotentialBoundaryConditions(),
      makeConstantSpeciesBoundaryConditions(species.size(), 1.0),
      std::make_unique<linalg::CholmodSolver>());

  CapturingFieldOutputWriter writer;
  std::vector<std::string> names;
  std::vector<std::string> categories;
  std::vector<pemu::trace::DiagDomain> domains;
  std::vector<pemu::trace::EventKind> kinds;
  std::vector<double> selected_time_steps;
  std::vector<double> completed_times;
  const RecordingTraceSink sink{.names = &names,
                                .categories = &categories,
                                .domains = &domains,
                                .kinds = &kinds,
                                .selected_time_steps = &selected_time_steps,
                                .completed_times = &completed_times};
  const FieldOutputOptions output_options{
      .writer = &writer,
      .directory = "simulation-output",
      .file_stem = "snapshot",
      .every_steps = 1,
      .write_initial = true,
      .write_final = true,
  };
  FixedStepPlasmaSimulation simulation(
      density, reactions, transport, NoReactionEvaluator{},
      FixedStepClock(dt, 2), sink, {}, output_options);

  simulation.run();

  ASSERT_EQ(writer.requests.size(), 3u);
  for (std::size_t index = 0; index < writer.requests.size(); ++index) {
    const auto& request = writer.requests[index];
    EXPECT_EQ(request.stamp.step, index);
    EXPECT_NEAR(request.stamp.time, dt * static_cast<double>(index), 1e-14);
    EXPECT_EQ(request.path,
              std::filesystem::path{"simulation-output"} / "snapshot.vtkhdf");
    EXPECT_EQ(request.cell_names,
              (std::vector<std::string>{
                  "species_0_e_number_density", "species_1_Ar+_number_density",
                  "charge_density", "electric_potential"}));
    EXPECT_EQ(request.face_names,
              (std::vector<std::string>{"species_0_e_normal_drift_velocity",
                                        "species_1_Ar+_normal_drift_velocity",
                                        "normal_electric_field"}));
  }

  EXPECT_EQ(std::count(domains.begin(), domains.end(),
                       pemu::trace::DiagDomain::output),
            1);
  EXPECT_EQ(std::count(names.begin(), names.end(), "completed"), 1);
  EXPECT_EQ(writer.finish_count, 1u);
}

TEST_F(FixedStepPlasmaSimulationTest,
       EmitsReturnedDiagnosticBeforeStepFailureTrace) {
  constexpr double dt = 0.01;
  physics::SpeciesSet species;
  (void)addElectronAndIon(species);
  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);
  physics::ReactionNetwork reactions(species);
  equation::FixedStepMultiSpeciesDriftDiffusionStepper transport(
      mesh_, species, 1.0, dt, makeZeroPotentialBoundaryConditions(),
      makeConstantSpeciesBoundaryConditions(species.size(), 1.0),
      std::make_unique<DiagnosticFailingSolver>());
  std::ostringstream output;
  trace::OstreamTraceSink sink(output);
  FixedStepPlasmaSimulation simulation(density, reactions, transport,
                                       NoReactionEvaluator{},
                                       FixedStepClock(dt, 1), sink);

  const auto result = simulation.advance();

  EXPECT_EQ(result.status, linalg::SolverStatus::SolveFailed);
  EXPECT_NE(output.str().find("Error    | linalg.test_backend.solve.failed"),
            std::string::npos)
      << output.str();
  EXPECT_NE(output.str().find("injected linear solve failure"),
            std::string::npos)
      << output.str();
  EXPECT_LT(output.str().find("linalg.test_backend.solve.failed"),
            output.str().find("simulation.fixed_step.step.failed"));
  EXPECT_EQ(simulation.step(), 0u);
  EXPECT_EQ(simulation.state(), SimulationState::failed);
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

  ASSERT_TRUE(simulation.advance().success());

  ASSERT_TRUE(simulation.finished());

  EXPECT_THROW((void)simulation.advance(), std::logic_error);
}

TEST_F(FixedStepPlasmaSimulationTest,
       StatisticsRejectFieldsWithoutPhysicalMetadata) {
  using namespace mp_units;
  using namespace mp_units::si::unit_symbols;

  constexpr double dt = 0.01;
  physics::SpeciesSet species;
  (void)addElectronAndIon(species);
  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);
  physics::ReactionNetwork reactions(species);
  equation::FixedStepMultiSpeciesDriftDiffusionStepper transport(
      mesh_, species, 1.0, dt, makeZeroPotentialBoundaryConditions(),
      makeConstantSpeciesBoundaryConditions(species.size(), 1.0),
      std::make_unique<linalg::CholmodSolver>());

  EXPECT_THROW(
      FixedStepPlasmaSimulation(
          density, reactions, transport, NoReactionEvaluator{},
          FixedStepClock(dt, 1), trace::NullTraceSink{},
          trace::StatisticsOptions{true, 1, 1.0 * cm, isq::length[cm]}),
      std::invalid_argument);
}

TEST_F(FixedStepPlasmaSimulationTest,
       CheckpointRejectsReorderedSpeciesWithoutChangingDensity) {
  constexpr double dt = 0.01;
  physics::SpeciesSet original_species;
  const auto original_ids = addElectronAndIon(original_species);
  physics::SpeciesCellFields original_density(mesh_, original_species.size(),
                                              0.0);
  original_density[original_ids.electron].fill(1.0);
  original_density[original_ids.ion].fill(2.0);
  physics::ReactionNetwork original_reactions(original_species);
  equation::FixedStepMultiSpeciesDriftDiffusionStepper original_transport(
      mesh_, original_species, 1.0, dt, makeZeroPotentialBoundaryConditions(),
      makeConstantSpeciesBoundaryConditions(original_species.size(), 1.0),
      std::make_unique<linalg::CholmodSolver>());
  FixedStepPlasmaSimulation original(original_density, original_reactions,
                                     original_transport, NoReactionEvaluator{},
                                     FixedStepClock(dt, 2));
  const auto checkpoint_path =
      std::filesystem::path{PEMU_SIMULATION_TEST_OUTPUT_DIR} /
      "species-identity-checkpoint" / "state.vtkhdf";
  const output::checkpoint::VtkHdfWriter writer;
  (void)original.saveCheckpoint(writer, checkpoint_path, true);

  physics::SpeciesSet reordered_species;
  (void)reordered_species.add(
      {.name = "Ar+",
       .charge = +1.0,
       .mobility = 0.5,
       .diffusivity = 0.1,
       .transport_model = physics::SpeciesTransportModel::DriftDiffusion});
  (void)reordered_species.add(
      {.name = "e",
       .charge = -1.0,
       .mobility = 1.0,
       .diffusivity = 0.1,
       .transport_model = physics::SpeciesTransportModel::DriftDiffusion});
  physics::SpeciesCellFields reordered_density(mesh_, reordered_species.size(),
                                               7.0);
  physics::ReactionNetwork reordered_reactions(reordered_species);
  equation::FixedStepMultiSpeciesDriftDiffusionStepper reordered_transport(
      mesh_, reordered_species, 1.0, dt, makeZeroPotentialBoundaryConditions(),
      makeConstantSpeciesBoundaryConditions(reordered_species.size(), 1.0),
      std::make_unique<linalg::CholmodSolver>());
  FixedStepPlasmaSimulation reordered(
      reordered_density, reordered_reactions, reordered_transport,
      NoReactionEvaluator{}, FixedStepClock(dt, 2));
  const output::checkpoint::VtkHdfReader reader;

  EXPECT_THROW((void)reordered.restoreCheckpoint(reader, checkpoint_path),
               std::invalid_argument);
  for (const auto& species_density : reordered_density) {
    for (const double value : species_density) {
      EXPECT_DOUBLE_EQ(value, 7.0);
    }
  }
  EXPECT_EQ(reordered.state(), SimulationState::ready);
  EXPECT_EQ(reordered.step(), 0u);
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
       CheckpointWorkflowRestoresClockDensityAndTimeStepHistory) {
  physics::SpeciesSet species;
  const auto ids = addElectronAndIon(species);
  physics::ReactionNetwork reactions(species);
  const auto ionization = reactions.addReaction(
      {.name = "ionization",
       .stoichiometry = {{ids.electron, +1.0}, {ids.ion, +1.0}}});
  const ElectronImpactIonizationEvaluator evaluator{.electron = ids.electron,
                                                    .ionization = ionization,
                                                    .neutral_density = 4.0,
                                                    .rate_coefficient = 0.5};
  constexpr equation::time_integration::AdaptiveTimeStepConfig time_config{
      .safety = 0.9, .min_dt = 1e-8, .max_dt = 0.04, .max_growth = 2.0};

  physics::SpeciesCellFields uninterrupted_density(mesh_, species.size(), 1.0);
  equation::AdaptiveStepMultiSpeciesDriftDiffusionStepper
      uninterrupted_transport(
          mesh_, species, 1.0, makeZeroPotentialBoundaryConditions(),
          makeConstantSpeciesBoundaryConditions(species.size(), 1.0),
          std::make_unique<linalg::CholmodSolver>(), time_config);
  AdaptiveStepPlasmaSimulation uninterrupted(uninterrupted_density, reactions,
                                             uninterrupted_transport, evaluator,
                                             AdaptiveTimeClock(0.1));
  uninterrupted.run();

  const auto checkpoint_path =
      std::filesystem::path{PEMU_SIMULATION_TEST_OUTPUT_DIR} /
      "adaptive-checkpoint-workflow" / "state.vtkhdf";
  const output::checkpoint::VtkHdfWriter writer;
  physics::SpeciesCellFields interrupted_density(mesh_, species.size(), 1.0);
  equation::AdaptiveStepMultiSpeciesDriftDiffusionStepper interrupted_transport(
      mesh_, species, 1.0, makeZeroPotentialBoundaryConditions(),
      makeConstantSpeciesBoundaryConditions(species.size(), 1.0),
      std::make_unique<linalg::CholmodSolver>(), time_config);
  AdaptiveStepPlasmaSimulation interrupted(
      interrupted_density, reactions, interrupted_transport, evaluator,
      AdaptiveTimeClock(0.1), trace::NullTraceSink{}, {}, {},
      {.writer = &writer,
       .path = checkpoint_path,
       .every_steps = 1,
       .write_final = false,
       .overwrite = true});
  ASSERT_TRUE(interrupted.advance().success());
  interrupted.pause();
  ASSERT_TRUE(std::filesystem::is_regular_file(checkpoint_path));

  physics::SpeciesCellFields restored_density(mesh_, species.size(), 0.0);
  equation::AdaptiveStepMultiSpeciesDriftDiffusionStepper restored_transport(
      mesh_, species, 1.0, makeZeroPotentialBoundaryConditions(),
      makeConstantSpeciesBoundaryConditions(species.size(), 1.0),
      std::make_unique<linalg::CholmodSolver>(), time_config);
  AdaptiveStepPlasmaSimulation restored(restored_density, reactions,
                                        restored_transport, evaluator,
                                        AdaptiveTimeClock(0.1));
  const CountingCheckpointReader reader;
  const auto record = restored.restoreCheckpoint(reader, checkpoint_path);

  EXPECT_EQ(reader.inspect_count, 0u);
  EXPECT_EQ(reader.restore_count, 1u);
  EXPECT_EQ(restored.state(), SimulationState::ready);
  EXPECT_EQ(restored.step(), 1u);
  EXPECT_DOUBLE_EQ(restored.time(), record.stamp.time);
  EXPECT_DOUBLE_EQ(restored.lastTimeStep(), interrupted.lastTimeStep());
  EXPECT_FALSE(restored.hasLastTimeStepProposal());
  restored.run();

  EXPECT_EQ(restored.state(), SimulationState::completed);
  EXPECT_EQ(restored.step(), uninterrupted.step());
  EXPECT_DOUBLE_EQ(restored.time(), uninterrupted.time());
  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {
    EXPECT_DOUBLE_EQ(restored_density[ids.electron][cell],
                     uninterrupted_density[ids.electron][cell]);
    EXPECT_DOUBLE_EQ(restored_density[ids.ion][cell],
                     uninterrupted_density[ids.ion][cell]);
  }
}

TEST_F(AdaptiveStepPlasmaSimulationTest,
       WritesScheduledSnapshotsAtAdaptiveStateTimes) {
  physics::SpeciesSet species;
  (void)addElectronAndIon(species);
  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);
  physics::ReactionNetwork reactions(species);
  equation::AdaptiveStepMultiSpeciesDriftDiffusionStepper transport(
      mesh_, species, 1.0, makeZeroPotentialBoundaryConditions(),
      makeConstantSpeciesBoundaryConditions(species.size(), 1.0),
      std::make_unique<linalg::CholmodSolver>(),
      {.safety = 0.9, .min_dt = 1e-8, .max_dt = 0.04, .max_growth = 2.0});

  CapturingFieldOutputWriter writer;
  const FieldOutputOptions output_options{
      .writer = &writer,
      .directory = "adaptive-output",
      .file_stem = "snapshot",
      .every_steps = 2,
      .write_initial = true,
      .write_final = true,
  };
  AdaptiveStepPlasmaSimulation simulation(
      density, reactions, transport, NoReactionEvaluator{},
      AdaptiveTimeClock(0.1), {}, {}, output_options);

  simulation.run();

  ASSERT_EQ(writer.requests.size(), 3u);
  EXPECT_EQ(writer.requests[0].stamp.step, 0u);
  EXPECT_NEAR(writer.requests[0].stamp.time, 0.0, 1e-14);
  EXPECT_EQ(writer.requests[1].stamp.step, 2u);
  EXPECT_NEAR(writer.requests[1].stamp.time, 0.08, 1e-14);
  EXPECT_EQ(writer.requests[2].stamp.step, 3u);
  EXPECT_NEAR(writer.requests[2].stamp.time, 0.1, 1e-14);
  EXPECT_EQ(writer.finish_count, 1u);
}

TEST_F(AdaptiveStepPlasmaSimulationTest,
       EmitsOrderedTraceEventsWithTimeStepDiagnostics) {
  physics::SpeciesSet species;
  (void)addElectronAndIon(species);
  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);
  physics::ReactionNetwork reactions(species);

  equation::AdaptiveStepMultiSpeciesDriftDiffusionStepper transport(
      mesh_, species, 1.0, makeZeroPotentialBoundaryConditions(),
      makeConstantSpeciesBoundaryConditions(species.size(), 1.0),
      std::make_unique<linalg::CholmodSolver>(),
      {.safety = 0.9, .min_dt = 1e-8, .max_dt = 0.04, .max_growth = 2.0});

  std::vector<std::string> names;
  std::vector<std::string> categories;
  std::vector<pemu::trace::DiagDomain> domains;
  std::vector<pemu::trace::EventKind> kinds;
  std::vector<double> selected_time_steps;
  std::vector<double> completed_times;
  names.reserve(20);
  categories.reserve(20);
  domains.reserve(20);
  kinds.reserve(20);
  selected_time_steps.reserve(3);
  completed_times.reserve(3);
  const RecordingTraceSink sink{.names = &names,
                                .categories = &categories,
                                .domains = &domains,
                                .kinds = &kinds,
                                .selected_time_steps = &selected_time_steps,
                                .completed_times = &completed_times};

  AdaptiveStepPlasmaSimulation simulation(density, reactions, transport,
                                          NoReactionEvaluator{},
                                          AdaptiveTimeClock(0.1), sink);

  simulation.run();

  const std::vector<std::string> expected_names{
      "run.started",
      "step.started",
      "electrostatics.completed",
      "reaction_rates.completed",
      "sources.completed",
      "timestep.selected",
      "step.completed",
      "step.started",
      "electrostatics.completed",
      "reaction_rates.completed",
      "sources.completed",
      "timestep.selected",
      "step.completed",
      "step.started",
      "electrostatics.completed",
      "reaction_rates.completed",
      "sources.completed",
      "timestep.selected",
      "step.completed",
      "run.completed",
  };
  EXPECT_EQ(names, expected_names);
  EXPECT_TRUE(std::ranges::all_of(categories, [](const std::string& category) {
    return category == "adaptive_step";
  }));
  EXPECT_TRUE(std::ranges::all_of(domains, [](pemu::trace::DiagDomain domain) {
    return domain == pemu::trace::DiagDomain::simulation;
  }));
  EXPECT_TRUE(std::ranges::all_of(kinds, [](pemu::trace::EventKind kind) {
    return kind == pemu::trace::EventKind::Trace;
  }));
  ASSERT_EQ(selected_time_steps.size(), 3u);
  EXPECT_NEAR(selected_time_steps[0], 0.04, 1e-14);
  EXPECT_NEAR(selected_time_steps[1], 0.04, 1e-14);
  EXPECT_NEAR(selected_time_steps[2], 0.02, 1e-14);
  ASSERT_EQ(completed_times.size(), 3u);
  EXPECT_NEAR(completed_times[0], 0.04, 1e-14);
  EXPECT_NEAR(completed_times[1], 0.08, 1e-14);
  EXPECT_NEAR(completed_times[2], 0.1, 1e-14);
}

TEST_F(AdaptiveStepPlasmaSimulationTest,
       EmitsReturnedDiagnosticBeforeStepFailureTrace) {
  physics::SpeciesSet species;
  (void)addElectronAndIon(species);
  physics::SpeciesCellFields density(mesh_, species.size(), 1.0);
  physics::ReactionNetwork reactions(species);
  equation::AdaptiveStepMultiSpeciesDriftDiffusionStepper transport(
      mesh_, species, 1.0, makeZeroPotentialBoundaryConditions(),
      makeConstantSpeciesBoundaryConditions(species.size(), 1.0),
      std::make_unique<DiagnosticFailingSolver>(),
      {.safety = 0.9, .min_dt = 1e-8, .max_dt = 0.04, .max_growth = 2.0});
  std::ostringstream output;
  trace::OstreamTraceSink sink(output);
  AdaptiveStepPlasmaSimulation simulation(density, reactions, transport,
                                          NoReactionEvaluator{},
                                          AdaptiveTimeClock(0.1), sink);

  const auto result = simulation.advance();

  EXPECT_EQ(result.status, linalg::SolverStatus::SolveFailed);
  EXPECT_NE(output.str().find("Error    | linalg.test_backend.solve.failed"),
            std::string::npos)
      << output.str();
  EXPECT_NE(output.str().find("injected linear solve failure"),
            std::string::npos)
      << output.str();
  EXPECT_LT(output.str().find("linalg.test_backend.solve.failed"),
            output.str().find("simulation.adaptive_step.step.failed"));
  EXPECT_EQ(simulation.step(), 0u);
  EXPECT_EQ(simulation.state(), SimulationState::failed);
}

TEST_F(AdaptiveStepPlasmaSimulationTest,
       AdvanceEvaluatesReactionAndUpdatesSpecies) {
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

  ASSERT_TRUE(simulation.advance().success());
  ASSERT_TRUE(simulation.hasLastTimeStepProposal());

  EXPECT_NEAR(simulation.time(), 0.04, 1e-14);
  EXPECT_NEAR(simulation.lastTimeStep(), 0.04, 1e-14);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {
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

  double evaluator_maximum{};
  AdaptiveStepPlasmaSimulation simulation(
      density, reactions, transport,
      ElectricFieldAwareEvaluator{.reaction = reaction,
                                  .observed_maximum = &evaluator_maximum},
      AdaptiveTimeClock(1e-3));

  ASSERT_TRUE(simulation.advance().success());
  EXPECT_NEAR(evaluator_maximum, 1.0, 1e-12);
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

  ASSERT_TRUE(simulation.advance().success());
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
  ASSERT_TRUE(simulation.advance().success());
  ASSERT_TRUE(simulation.finished());
  EXPECT_THROW((void)simulation.advance(), std::logic_error);
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
    EXPECT_NEAR(density[ids.electron][cell], expected_density, 1e-12);
    EXPECT_NEAR(density[ids.ion][cell], expected_density, 1e-12);
    EXPECT_NEAR(transport.chargeDensity()[cell], 0.0, 1e-12);
    EXPECT_NEAR(transport.potential()[cell], 0.0, 1e-12);
  }

  for (const double electric_field : transport.electricFieldNormal()) {
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
  constexpr auto end_time = 1.0e-6 * s;
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
  EXPECT_EQ(*density[electron].metadata().physical_quantity,
            pemu::unit::bridgeReference(
                pemu::unit::plasma_quantity::particle_number_density
                    [number_density_unit]));
  EXPECT_EQ(*transport.potential().metadata().physical_quantity,
            pemu::unit::bridgeReference(isq::electric_potential[V]));
  EXPECT_EQ(
      *transport.electricFieldNormal().metadata().physical_quantity,
      pemu::unit::bridgeReference(
          pemu::unit::plasma_quantity::normal_electric_field_strength[V / cm]));

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
  // Keep the multi-step output out of the test runner's terminal output and
  // separate ordinary diagnostics from field statistics.
  std::ofstream diagnostic_output("test.diag.log");
  std::ofstream statistics_output("test.statistics.log");
  ASSERT_TRUE(diagnostic_output.is_open());
  ASSERT_TRUE(statistics_output.is_open());
  trace::SplitTraceSink trace_sink{trace::OstreamTraceSink{diagnostic_output},
                                   trace::OstreamTraceSink{statistics_output}};
  output::dump::VtkHdfWriter field_output_writer;
  const auto field_output_directory =
      std::filesystem::path{PEMU_SIMULATION_TEST_OUTPUT_DIR} /
      "parallel-plate-400v";
  const FieldOutputOptions field_output_options{
      .writer = &field_output_writer,
      .directory = field_output_directory,
      .file_stem = "parallel-plate-400v",
      .every_steps = 1,
      .write_initial = true,
      .write_final = true,
      .overwrite = true,
  };
  AdaptiveStepPlasmaSimulation simulation(
      density, reactions, transport, evaluator, AdaptiveTimeClock(end_time_s),
      trace_sink,
      trace::StatisticsOptions{true, 1, domain_length_cm * cm, isq::length[cm]},
      field_output_options);

  simulation.run();

  ASSERT_TRUE(simulation.finished());
  ASSERT_TRUE(simulation.hasLastTimeStepProposal());
  EXPECT_NEAR(simulation.time(), end_time_s, 1.0e-18);
  EXPECT_GT(simulation.step(), 1u);
  EXPECT_LT(simulation.lastTimeStep(), max_time_step_s);
  const auto field_output_path =
      field_output_directory / "parallel-plate-400v.vtkhdf";
  ASSERT_TRUE(std::filesystem::is_regular_file(field_output_path));
  vtkNew<vtkHDFReader> field_output_reader;
  ASSERT_TRUE(
      field_output_reader->CanReadFile(field_output_path.string().c_str()));
  field_output_reader->SetFileName(field_output_path.string().c_str());
  field_output_reader->Update();
  EXPECT_EQ(field_output_reader->GetNumberOfSteps(),
            static_cast<int>(simulation.step() + 1));
  field_output_reader->SetStep(static_cast<int>(simulation.step()));
  field_output_reader->Update();
  auto* final_output_grid = vtkUnstructuredGrid::SafeDownCast(
      field_output_reader->GetOutputDataObject(0));
  ASSERT_NE(final_output_grid, nullptr);
  auto* final_output_step =
      final_output_grid->GetFieldData()->GetArray("pemu_step");
  ASSERT_NE(final_output_step, nullptr);
  EXPECT_EQ(static_cast<std::uint64_t>(final_output_step->GetTuple1(0)),
            simulation.step());

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
  }

  diagnostic_output.close();
  statistics_output.close();
  std::ifstream diagnostic_input("test.diag.log");
  std::ifstream statistics_input("test.statistics.log");
  ASSERT_TRUE(diagnostic_input.is_open());
  ASSERT_TRUE(statistics_input.is_open());
  std::ostringstream diagnostic_contents;
  std::ostringstream statistics_contents;
  diagnostic_contents << diagnostic_input.rdbuf();
  statistics_contents << statistics_input.rdbuf();
  const std::string diagnostic_text = diagnostic_contents.str();
  const std::string statistics_text = statistics_contents.str();
  EXPECT_NE(diagnostic_text.find("simulation.adaptive_step.run.started"),
            std::string::npos);
  EXPECT_NE(diagnostic_text.find("Success"), std::string::npos)
      << diagnostic_text;
  EXPECT_EQ(diagnostic_text.find("STATISTICS"), std::string::npos);
  EXPECT_EQ(statistics_text.find("simulation.adaptive_step"),
            std::string::npos);
  EXPECT_NE(statistics_text.find("STATISTICS | STEP=0"), std::string::npos);
  EXPECT_NE(statistics_text.find("FIELD"), std::string::npos);
  EXPECT_NE(statistics_text.find("UNIT"), std::string::npos);
  EXPECT_NE(statistics_text.find("species[e]"), std::string::npos);
  EXPECT_NE(statistics_text.find("species[Ar+]"), std::string::npos);
  EXPECT_NE(statistics_text.find("charge_density"), std::string::npos);
  EXPECT_NE(statistics_text.find("potential"), std::string::npos);
  EXPECT_NE(statistics_text.find("electric_field_normal"), std::string::npos);
  EXPECT_NE(statistics_text.find("VOLUME_SEMANTICS=planar_extrusion"),
            std::string::npos);
  EXPECT_NE(statistics_text.find("PHYSICAL_VOLUME=1 [mL]"), std::string::npos);
  EXPECT_NE(statistics_text.find("PHYSICAL_FACE_AREA=130 [cm^2]"),
            std::string::npos);
  EXPECT_NE(statistics_text.find("1e+06 [1]"), std::string::npos);
  EXPECT_NE(statistics_text.find("0 [C]"), std::string::npos);
  EXPECT_EQ(statistics_text.find("DETAILS"), std::string::npos);
  EXPECT_EQ(statistics_text.find("non_finite=0"), std::string::npos);
  EXPECT_EQ(statistics_text.find("negative=0"), std::string::npos);
}

TEST_F(PlasmaSimulation64x64CheckpointTest,
       ParallelPlate400VCheckpointRestoresAndContinuesSimulation) {
  using namespace mp_units;
  using namespace mp_units::si::unit_symbols;

  constexpr auto number_density_unit = one / cubic(cm);
  constexpr auto mobility_unit = square(cm) / (V * s);
  constexpr auto diffusivity_unit = square(cm) / s;
  constexpr auto ionization_coefficient_unit = cubic(cm) / s;
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
  constexpr auto fixed_time_step = 1.0e-8 * s;
  constexpr std::size_t checkpoint_step = 3;
  constexpr std::size_t total_steps = 8;

  constexpr double left_voltage = 0.0;
  constexpr double right_voltage = 400.0;
  constexpr double elementary_charge_coulomb =
      elementary_charge.numerical_value_in(C);
  constexpr double initial_density_cm3 =
      initial_density.numerical_value_in(number_density_unit);
  constexpr double dt = fixed_time_step.numerical_value_in(s);
  const auto field_metadata = field::centimetrePlasmaFieldMetadata();

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

  physics::ReactionNetwork reactions(species);
  const auto ionization =
      reactions.addReaction({.name = "electron-impact ionization",
                             .stoichiometry = {{electron, +1.0}, {ion, +1.0}}});
  const ElectronImpactIonizationEvaluator evaluator{
      .electron = electron,
      .ionization = ionization,
      .neutral_density =
          neutral_density.numerical_value_in(number_density_unit),
      .rate_coefficient = ionization_rate_coefficient.numerical_value_in(
          ionization_coefficient_unit)};

  physics::SpeciesCellFields uninterrupted_density(
      mesh_, species.size(), initial_density_cm3,
      field_metadata.number_density);
  equation::FixedStepMultiSpeciesDriftDiffusionStepper uninterrupted_transport(
      mesh_, species, vacuum_permittivity.numerical_value_in(F / cm), dt,
      makeParallelPlatePotentialBoundaryConditions(left_voltage, right_voltage),
      makeConstantSpeciesBoundaryConditions(species.size(),
                                            initial_density_cm3),
      std::make_unique<linalg::CholmodSolver>(), field_metadata);
  FixedStepPlasmaSimulation uninterrupted(uninterrupted_density, reactions,
                                          uninterrupted_transport, evaluator,
                                          FixedStepClock(dt, total_steps));
  uninterrupted.run();

  physics::SpeciesCellFields interrupted_density(mesh_, species.size(),
                                                 initial_density_cm3,
                                                 field_metadata.number_density);
  equation::FixedStepMultiSpeciesDriftDiffusionStepper interrupted_transport(
      mesh_, species, vacuum_permittivity.numerical_value_in(F / cm), dt,
      makeParallelPlatePotentialBoundaryConditions(left_voltage, right_voltage),
      makeConstantSpeciesBoundaryConditions(species.size(),
                                            initial_density_cm3),
      std::make_unique<linalg::CholmodSolver>(), field_metadata);
  const auto checkpoint_path =
      std::filesystem::path{PEMU_SIMULATION_TEST_OUTPUT_DIR} /
      "parallel-plate-400v-checkpoint" / "state.vtkhdf";
  const output::checkpoint::VtkHdfWriter checkpoint_writer;
  FixedStepPlasmaSimulation interrupted(
      interrupted_density, reactions, interrupted_transport, evaluator,
      FixedStepClock(dt, total_steps), trace::NullTraceSink{}, {}, {},
      {.writer = &checkpoint_writer,
       .path = checkpoint_path,
       .every_steps = checkpoint_step,
       .write_final = false,
       .overwrite = true});
  for (std::size_t step = 0; step < checkpoint_step; ++step) {
    ASSERT_TRUE(interrupted.advance().success());
  }
  ASSERT_EQ(interrupted.step(), checkpoint_step);
  interrupted.pause();
  ASSERT_EQ(interrupted.state(), SimulationState::paused);
  ASSERT_TRUE(std::filesystem::is_regular_file(checkpoint_path));

  physics::SpeciesCellFields restored_density(mesh_, species.size(), 0.0,
                                              field_metadata.number_density);
  equation::FixedStepMultiSpeciesDriftDiffusionStepper restored_transport(
      mesh_, species, vacuum_permittivity.numerical_value_in(F / cm), dt,
      makeParallelPlatePotentialBoundaryConditions(left_voltage, right_voltage),
      makeConstantSpeciesBoundaryConditions(species.size(),
                                            initial_density_cm3),
      std::make_unique<linalg::CholmodSolver>(), field_metadata);
  FixedStepPlasmaSimulation restored(restored_density, reactions,
                                     restored_transport, evaluator,
                                     FixedStepClock(dt, total_steps));
  const CountingCheckpointReader checkpoint_reader;
  const auto restored_record =
      restored.restoreCheckpoint(checkpoint_reader, checkpoint_path);
  ASSERT_EQ(checkpoint_reader.inspect_count, 0u);
  ASSERT_EQ(checkpoint_reader.restore_count, 1u);
  ASSERT_EQ(restored_record.stamp.step, checkpoint_step);
  ASSERT_DOUBLE_EQ(restored_record.stamp.time,
                   static_cast<double>(checkpoint_step) * dt);
  ASSERT_EQ(restored.state(), SimulationState::ready);
  ASSERT_DOUBLE_EQ(restored.time(), restored_record.stamp.time);
  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {
    EXPECT_DOUBLE_EQ(restored_density[electron][cell],
                     interrupted_density[electron][cell]);
    EXPECT_DOUBLE_EQ(restored_density[ion][cell],
                     interrupted_density[ion][cell]);
  }
  restored.run();

  ASSERT_TRUE(uninterrupted.finished());
  ASSERT_TRUE(restored.finished());
  ASSERT_EQ(restored.state(), SimulationState::completed);
  ASSERT_EQ(restored.step(), uninterrupted.step());
  ASSERT_DOUBLE_EQ(restored.time(), uninterrupted.time());
  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {
    EXPECT_DOUBLE_EQ(restored_density[electron][cell],
                     uninterrupted_density[electron][cell]);
    EXPECT_DOUBLE_EQ(restored_density[ion][cell],
                     uninterrupted_density[ion][cell]);
    EXPECT_DOUBLE_EQ(restored_transport.chargeDensity()[cell],
                     uninterrupted_transport.chargeDensity()[cell]);
    EXPECT_DOUBLE_EQ(restored_transport.potential()[cell],
                     uninterrupted_transport.potential()[cell]);
  }
  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {
    EXPECT_DOUBLE_EQ(restored_transport.electricFieldNormal()[face],
                     uninterrupted_transport.electricFieldNormal()[face]);
    EXPECT_DOUBLE_EQ(
        restored_transport.driftVelocityNormal(electron)[face],
        uninterrupted_transport.driftVelocityNormal(electron)[face]);
    EXPECT_DOUBLE_EQ(restored_transport.driftVelocityNormal(ion)[face],
                     uninterrupted_transport.driftVelocityNormal(ion)[face]);
  }
}

TEST_F(AdaptiveStepPlasmaSimulation64x64Test,
       ParallelPlate400VPreservesLowerSolverDiagnosticAfterPhysicalStep) {
  using namespace mp_units;
  using namespace mp_units::si::unit_symbols;

  constexpr auto number_density_unit = one / cubic(cm);
  constexpr auto mobility_unit = square(cm) / (V * s);
  constexpr auto diffusivity_unit = square(cm) / s;
  constexpr auto ionization_coefficient_unit = cubic(cm) / s;

  // Keep this physical configuration equal to the successful parallel-plate
  // integration test above.  The only changed component is the injected
  // linear-solver backend.
  constexpr auto left_potential = 0.0 * V;
  constexpr auto right_potential = 400.0 * V;
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
      std::make_unique<FailAfterSuccessfulCholmodSolves>(1u),
      {.safety = 0.8,
       .min_dt = 1.0e-12,
       .max_dt = max_time_step_s,
       .max_growth = 1.5},
      field_metadata);

  ElectronImpactIonizationEvaluator evaluator{
      .electron = electron,
      .ionization = ionization,
      .neutral_density =
          neutral_density.numerical_value_in(number_density_unit),
      .rate_coefficient = ionization_rate_coefficient.numerical_value_in(
          ionization_coefficient_unit)};
  // Keep this injected failure separate from the successful 400 V run above.
  std::ofstream trace_output("test.failure.diag.log");
  ASSERT_TRUE(trace_output.is_open());
  trace::OstreamTraceSink trace_sink(trace_output);
  AdaptiveStepPlasmaSimulation simulation(
      density, reactions, transport, evaluator, AdaptiveTimeClock(end_time_s),
      trace_sink);

  EXPECT_THROW(simulation.run(), std::runtime_error);

  EXPECT_EQ(simulation.step(), 1u);
  EXPECT_GT(simulation.time(), 0.0);
  EXPECT_FALSE(simulation.finished());

  trace_output.close();
  std::ifstream trace_input("test.failure.diag.log");
  ASSERT_TRUE(trace_input.is_open());
  std::ostringstream trace_contents;
  trace_contents << trace_input.rdbuf();
  const std::string trace_text = trace_contents.str();
  const auto diagnostic_position =
      trace_text.find("Error    | linalg.test_backend.solve.failed");
  EXPECT_NE(diagnostic_position, std::string::npos) << trace_text;
  EXPECT_NE(trace_text.find("SolveFailed"), std::string::npos) << trace_text;
  EXPECT_NE(trace_text.find("injected failure after successful CHOLMOD solve"),
            std::string::npos)
      << trace_text;

  const auto failure_position =
      trace_text.find("simulation.adaptive_step.step.failed");
  EXPECT_NE(failure_position, std::string::npos) << trace_text;
  EXPECT_LT(diagnostic_position, failure_position) << trace_text;
}

}  // namespace pemu::simulation::test
