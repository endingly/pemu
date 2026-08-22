#include <pemu/output/dump/trace.hpp>
#include <pemu/simulation/adaptive_step_plasma_simulation.hpp>
#include <pemu/simulation/detail/plasma_checkpoint.hpp>
#include <pemu/simulation/detail/plasma_field_output.hpp>
#include <pemu/simulation/detail/plasma_statistics.hpp>
#include <pemu/simulation/fixed_step_plasma_simulation.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>

namespace pemu::simulation {
namespace {

class TraceAdapter {
 public:
  explicit TraceAdapter(const PlasmaTraceSink* sink) : sink_(sink) {}

  void operator()(const trace::TraceEvent& event) const noexcept {
    if (!sink_ || !*sink_) {
      return;
    }
    try {
      (*sink_)(event);
    } catch (...) {}
  }

 private:
  const PlasmaTraceSink* sink_{};
};

template <std::size_t N>
void emitTrace(
    const PlasmaTraceSink& sink, std::string_view category,
    std::string_view name, trace::Severity severity,
    const std::array<trace::TraceAttribute, N>& attributes) noexcept {
  TraceAdapter{&sink}({.domain = trace::DiagDomain::simulation,
                       .category = category,
                       .name = name,
                       .severity = severity,
                       .attributes = attributes});
}

void emitDiagnostic(const PlasmaTraceSink& sink,
                    const linalg::SolverResult& result, std::size_t step,
                    double time) noexcept {
  if (!result.diagnostic.has_value()) {
    return;
  }
  const std::array attributes{
      trace::TraceAttribute{"step", static_cast<std::uint64_t>(step)},
      trace::TraceAttribute{"solver_status", to_string(result.status)},
      trace::TraceAttribute{"time", time},
      trace::TraceAttribute{"residual_norm", result.residual_norm},
      trace::TraceAttribute{"relative_residual", result.relative_residual},
  };
  auto diagnostic = *result.diagnostic;
  diagnostic.attributes = attributes;
  TraceAdapter{&sink}(diagnostic);
}

void emitSolverResult(const PlasmaTraceSink& sink, std::string_view category,
                      std::string_view name, trace::Severity severity,
                      const linalg::SolverResult& result,
                      std::size_t step) noexcept {
  const std::array attributes{
      trace::TraceAttribute{"step", static_cast<std::uint64_t>(step)},
      trace::TraceAttribute{"solver_status", to_string(result.status)},
      trace::TraceAttribute{"residual_norm", result.residual_norm},
      trace::TraceAttribute{"relative_residual", result.relative_residual},
  };
  emitTrace(sink, category, name, severity, attributes);
}

}  // namespace

struct FixedStepPlasmaSimulation::Impl {
  physics::SpeciesCellFields* density;
  const physics::ReactionNetwork* reaction_network;
  Stepper* transport_stepper;
  PlasmaReactionRateEvaluator rate_evaluator;
  FixedStepClock clock;
  PlasmaTraceSink trace_sink;
  trace::StatisticsOptions statistics_options;
  FieldOutputOptions field_output_options;
  CheckpointOptions checkpoint_options;
  physics::ReactionRateFields reaction_rates;
  physics::SpeciesCellFields source;
  std::unique_ptr<output::dump::ISeries> field_output_series;
  bool field_output_has_snapshots{};
  SimulationState state{SimulationState::ready};

  Impl(physics::SpeciesCellFields& density_in,
       const physics::ReactionNetwork& reactions, Stepper& stepper,
       PlasmaReactionRateEvaluator evaluator, FixedStepClock clock_in,
       PlasmaTraceSink sink, trace::StatisticsOptions statistics,
       FieldOutputOptions field_output, CheckpointOptions checkpoint)
      : density(&density_in),
        reaction_network(&reactions),
        transport_stepper(&stepper),
        rate_evaluator(std::move(evaluator)),
        clock(std::move(clock_in)),
        trace_sink(std::move(sink)),
        statistics_options(statistics),
        field_output_options(std::move(field_output)),
        checkpoint_options(std::move(checkpoint)),
        reaction_rates(density_in.mesh(), reactions.size(), 0.0,
                       stepper.fieldMetadata().reaction_rate),
        source(density_in.mesh(), density_in.size(), 0.0,
               stepper.fieldMetadata().number_density_source) {
    if (!rate_evaluator) {
      throw std::invalid_argument(
          "simulation reaction evaluator must not be empty");
    }
    if (density->size() == 0) {
      throw std::invalid_argument("simulation density must contain species");
    }
    const double scale = std::max(
        {1.0, std::abs(clock.timeStep()), std::abs(stepper.timeStep())});
    if (std::abs(clock.timeStep() - stepper.timeStep()) > 1e-14 * scale) {
      throw std::invalid_argument(
          "FixedStepClock dt does not match FixedStep transport dt");
    }
    detail::validatePlasmaStatisticsConfiguration(
        statistics_options, *density, stepper.chargeDensity(),
        stepper.potential(), stepper.electricFieldNormal());
    validateFieldOutputOptions(field_output_options);
    validateCheckpointOptions(checkpoint_options);
  }

  void emitRun(std::string_view name, trace::Severity severity) const noexcept {
    const std::array attributes{
        trace::TraceAttribute{"step", static_cast<std::uint64_t>(clock.step())},
        trace::TraceAttribute{"time", clock.time()},
        trace::TraceAttribute{"end_time", clock.endTime()},
    };
    emitTrace(trace_sink, "fixed_step", name, severity, attributes);
  }

  void finishOutput() {
    if (!field_output_series) {
      return;
    }
    if (field_output_has_snapshots) {
      const auto record = field_output_series->finish();
      TraceAdapter output_sink{&trace_sink};
      output::dump::traceCompleted(output_sink, record);
    }
    field_output_series.reset();
  }

  void ensureOutput() {
    if (!field_output_options.enabled() || field_output_series) {
      return;
    }
    field_output_series = field_output_options.writer->openSeries(
        {.mesh = &transport_stepper->mesh(),
         .path = detail::plasmaFieldOutputPath(field_output_options),
         .overwrite = field_output_options.overwrite});
  }

  void writeOutput(bool terminal) {
    if (!field_output_options.enabled()) {
      return;
    }
    if (detail::shouldWritePlasmaFieldSnapshot(field_output_options,
                                               clock.step(), terminal)) {
      (void)detail::writePlasmaFieldSnapshot(
          *density, *transport_stepper, field_output_options,
          *field_output_series,
          {.step = static_cast<std::uint64_t>(clock.step()),
           .time = clock.time()});
      field_output_has_snapshots = true;
    }
    if (terminal && field_output_has_snapshots) {
      finishOutput();
    }
  }

  output::OutputRecord save(const output::checkpoint::IWriter& writer,
                            const std::filesystem::path& path,
                            bool overwrite) const {
    if (state == SimulationState::failed || state == SimulationState::stopped) {
      throw std::logic_error(
          "terminal simulation state cannot be checkpointed");
    }
    const auto sources = detail::plasmaCheckpointSources(*density);
    return writer.write(
        {.mesh = &transport_stepper->mesh(),
         .path = path,
         .stamp = {.step = static_cast<std::uint64_t>(clock.step()),
                   .time = clock.time()},
         .cell_fields = sources,
         .overwrite = overwrite});
  }

  void complete() {
    if (state != SimulationState::completed) {
      state = SimulationState::completed;
      emitRun("run.completed", trace::Severity::Info);
    }
  }

  void fail() noexcept {
    if (state != SimulationState::failed) {
      state = SimulationState::failed;
      emitRun("run.failed", trace::Severity::Error);
    }
  }

  linalg::SolverResult perform() {
    if (clock.finished()) {
      throw std::out_of_range("fixed-step simulation already finished");
    }
    const std::array started{
        trace::TraceAttribute{"step", static_cast<std::uint64_t>(clock.step())},
        trace::TraceAttribute{"time", clock.time()},
        trace::TraceAttribute{"dt", clock.timeStep()},
    };
    emitTrace(trace_sink, "fixed_step", "step.started", trace::Severity::Trace,
              started);
    const auto result = transport_stepper->updateElectrostatics(*density);
    if (!result.success()) {
      emitDiagnostic(trace_sink, result, clock.step(), clock.time());
      emitSolverResult(trace_sink, "fixed_step", "step.failed",
                       trace::Severity::Error, result, clock.step());
      return result;
    }
    emitSolverResult(trace_sink, "fixed_step", "electrostatics.completed",
                     trace::Severity::Trace, result, clock.step());
    writeOutput(false);
    reaction_rates.fill(0.0);
    rate_evaluator(*density, transport_stepper->potential(),
                   transport_stepper->electricFieldNormal(), reaction_rates);
    const std::array rate_attributes{
        trace::TraceAttribute{"step", static_cast<std::uint64_t>(clock.step())},
        trace::TraceAttribute{
            "field_count", static_cast<std::uint64_t>(reaction_rates.size())},
    };
    emitTrace(trace_sink, "fixed_step", "reaction_rates.completed",
              trace::Severity::Trace, rate_attributes);
    source.fill(0.0);
    reaction_network->accumulateSources(reaction_rates.span(), source);
    const std::array source_attributes{
        trace::TraceAttribute{"step", static_cast<std::uint64_t>(clock.step())},
        trace::TraceAttribute{"field_count",
                              static_cast<std::uint64_t>(source.size())},
    };
    emitTrace(trace_sink, "fixed_step", "sources.completed",
              trace::Severity::Trace, source_attributes);
    if (statistics_options.enabled()) {
      TraceAdapter statistics_sink{&trace_sink};
      detail::emitPlasmaStatistics(
          statistics_sink, transport_stepper->species(), *density,
          transport_stepper->chargeDensity(), transport_stepper->potential(),
          transport_stepper->electricFieldNormal(), statistics_options,
          clock.step(), clock.time());
    }
    transport_stepper->advanceTransport(*density, source);
    clock.advance();
    if (clock.finished() && field_output_options.enabled()) {
      const auto final_result =
          transport_stepper->updateElectrostatics(*density);
      if (!final_result.success()) {
        emitDiagnostic(trace_sink, final_result, clock.step(), clock.time());
        emitSolverResult(trace_sink, "fixed_step",
                         "final_output.electrostatics.failed",
                         trace::Severity::Error, final_result, clock.step());
        return final_result;
      }
      writeOutput(true);
    }
    const std::array complete_attributes{
        trace::TraceAttribute{"step", static_cast<std::uint64_t>(clock.step())},
        trace::TraceAttribute{"time", clock.time()},
        trace::TraceAttribute{"dt", clock.timeStep()},
        trace::TraceAttribute{"relative_residual", result.relative_residual},
    };
    emitTrace(trace_sink, "fixed_step", "step.completed", trace::Severity::Info,
              complete_attributes);
    return result;
  }
};

FixedStepPlasmaSimulation::FixedStepPlasmaSimulation(
    physics::SpeciesCellFields& density,
    const physics::ReactionNetwork& reactions, Stepper& transport,
    PlasmaReactionRateEvaluator evaluator, FixedStepClock clock,
    PlasmaTraceSink sink, trace::StatisticsOptions statistics,
    FieldOutputOptions field_output, CheckpointOptions checkpoint)
    : impl_(std::make_unique<Impl>(
          density, reactions, transport, std::move(evaluator), std::move(clock),
          std::move(sink), statistics, std::move(field_output),
          std::move(checkpoint))) {}
FixedStepPlasmaSimulation::~FixedStepPlasmaSimulation() = default;
FixedStepPlasmaSimulation::FixedStepPlasmaSimulation(
    FixedStepPlasmaSimulation&&) noexcept = default;
FixedStepPlasmaSimulation& FixedStepPlasmaSimulation::operator=(
    FixedStepPlasmaSimulation&&) noexcept = default;

void FixedStepPlasmaSimulation::start() {
  auto& self = *impl_;
  if (self.state == SimulationState::running)
    return;
  if (self.state != SimulationState::ready &&
      self.state != SimulationState::paused)
    throw std::logic_error("simulation cannot start from its current state");
  const bool resumed = self.state == SimulationState::paused;
  if (!resumed && !self.clock.finished()) {
    try {
      self.ensureOutput();
    } catch (...) {
      self.fail();
      throw;
    }
  }
  self.state = SimulationState::running;
  self.emitRun(resumed ? "run.resumed" : "run.started", trace::Severity::Info);
  if (self.clock.finished())
    self.complete();
}
void FixedStepPlasmaSimulation::pause() {
  auto& self = *impl_;
  if (self.state == SimulationState::paused)
    return;
  if (self.state != SimulationState::running)
    throw std::logic_error("only a running simulation can be paused");
  self.state = SimulationState::paused;
  self.emitRun("run.paused", trace::Severity::Info);
}
void FixedStepPlasmaSimulation::stop() {
  auto& self = *impl_;
  if (self.state == SimulationState::stopped)
    return;
  if (self.state == SimulationState::completed ||
      self.state == SimulationState::failed)
    throw std::logic_error("terminal simulation cannot be stopped");
  self.state = SimulationState::stopped;
  self.finishOutput();
  self.emitRun("run.stopped", trace::Severity::Info);
}
linalg::SolverResult FixedStepPlasmaSimulation::advance() {
  auto& self = *impl_;
  if (self.state == SimulationState::ready)
    start();
  if (self.state != SimulationState::running)
    throw std::logic_error("simulation step requires the running state");
  try {
    const auto result = self.perform();
    if (!result.success()) {
      self.fail();
      return result;
    }
    if (shouldWriteCheckpoint(self.checkpoint_options, self.clock.step(),
                              self.clock.finished()))
      (void)saveCheckpoint();
    if (self.clock.finished())
      self.complete();
    return result;
  } catch (...) {
    self.fail();
    throw;
  }
}
void FixedStepPlasmaSimulation::run() {
  if (impl_->state == SimulationState::ready ||
      impl_->state == SimulationState::paused)
    start();
  else if (impl_->state == SimulationState::completed)
    return;
  else if (impl_->state != SimulationState::running)
    throw std::logic_error("simulation cannot run from its current state");
  while (impl_->state == SimulationState::running) {
    if (!advance().success())
      throw std::runtime_error("fixed-step plasma simulation failed");
  }
}
SimulationState FixedStepPlasmaSimulation::state() const noexcept {
  return impl_->state;
}
double FixedStepPlasmaSimulation::time() const noexcept {
  return impl_->clock.time();
}
double FixedStepPlasmaSimulation::timeStep() const noexcept {
  return impl_->transport_stepper->timeStep();
}
std::size_t FixedStepPlasmaSimulation::step() const noexcept {
  return impl_->clock.step();
}
std::size_t FixedStepPlasmaSimulation::totalSteps() const noexcept {
  return impl_->clock.totalSteps();
}
bool FixedStepPlasmaSimulation::finished() const noexcept {
  return impl_->clock.finished();
}
output::OutputRecord FixedStepPlasmaSimulation::saveCheckpoint() const {
  if (!impl_->checkpoint_options.enabled())
    throw std::logic_error("simulation checkpoint output is disabled");
  return saveCheckpoint(*impl_->checkpoint_options.writer,
                        impl_->checkpoint_options.path,
                        impl_->checkpoint_options.overwrite);
}
output::OutputRecord FixedStepPlasmaSimulation::saveCheckpoint(
    const output::checkpoint::IWriter& writer,
    const std::filesystem::path& path, bool overwrite) const {
  return impl_->save(writer, path, overwrite);
}
output::OutputRecord FixedStepPlasmaSimulation::restoreCheckpoint(
    const output::checkpoint::IReader& reader,
    const std::filesystem::path& path) {
  auto& self = *impl_;
  if (self.state != SimulationState::ready)
    throw std::logic_error("checkpoint restore requires a ready simulation");
  const auto manifest = reader.inspect(path);
  if (manifest.stamp.step > self.clock.totalSteps())
    throw std::invalid_argument(
        "checkpoint step exceeds fixed simulation horizon");
  const auto restored_step = static_cast<std::size_t>(manifest.stamp.step);
  const double expected =
      static_cast<double>(restored_step) * self.clock.timeStep();
  if (std::abs(manifest.stamp.time - expected) >
      1e-14 *
          std::max({1.0, std::abs(expected), std::abs(manifest.stamp.time)}))
    throw std::invalid_argument(
        "checkpoint time does not match fixed simulation step");
  auto targets = detail::plasmaCheckpointTargets(*self.density);
  const auto record = reader.restore(
      path, {.mesh = &self.transport_stepper->mesh(), .cell_fields = targets});
  self.clock = FixedStepClock(self.clock.timeStep(), self.clock.totalSteps(),
                              restored_step);
  return record;
}

struct AdaptiveStepPlasmaSimulation::Impl {
  physics::SpeciesCellFields* density;
  const physics::ReactionNetwork* reaction_network;
  Stepper* transport_stepper;
  PlasmaReactionRateEvaluator rate_evaluator;
  AdaptiveTimeClock clock;
  PlasmaTraceSink trace_sink;
  trace::StatisticsOptions statistics_options;
  FieldOutputOptions field_output_options;
  CheckpointOptions checkpoint_options;
  physics::ReactionRateFields reaction_rates;
  physics::SpeciesCellFields source;
  std::unique_ptr<output::dump::ISeries> field_output_series;
  bool field_output_has_snapshots{};
  SimulationState state{SimulationState::ready};

  Impl(physics::SpeciesCellFields& density_in,
       const physics::ReactionNetwork& reactions, Stepper& stepper,
       PlasmaReactionRateEvaluator evaluator, AdaptiveTimeClock clock_in,
       PlasmaTraceSink sink, trace::StatisticsOptions statistics,
       FieldOutputOptions field_output, CheckpointOptions checkpoint)
      : density(&density_in),
        reaction_network(&reactions),
        transport_stepper(&stepper),
        rate_evaluator(std::move(evaluator)),
        clock(std::move(clock_in)),
        trace_sink(std::move(sink)),
        statistics_options(statistics),
        field_output_options(std::move(field_output)),
        checkpoint_options(std::move(checkpoint)),
        reaction_rates(density_in.mesh(), reactions.size(), 0.0,
                       stepper.fieldMetadata().reaction_rate),
        source(density_in.mesh(), density_in.size(), 0.0,
               stepper.fieldMetadata().number_density_source) {
    if (!rate_evaluator)
      throw std::invalid_argument(
          "simulation reaction evaluator must not be empty");
    if (density->size() == 0)
      throw std::invalid_argument("simulation density must contain species");
    detail::validatePlasmaStatisticsConfiguration(
        statistics_options, *density, stepper.chargeDensity(),
        stepper.potential(), stepper.electricFieldNormal());
    validateFieldOutputOptions(field_output_options);
    validateCheckpointOptions(checkpoint_options);
  }
  void emitRun(std::string_view name, trace::Severity severity) const noexcept {
    const std::array attributes{
        trace::TraceAttribute{"step", static_cast<std::uint64_t>(clock.step())},
        trace::TraceAttribute{"time", clock.time()},
        trace::TraceAttribute{"end_time", clock.endTime()},
    };
    emitTrace(trace_sink, "adaptive_step", name, severity, attributes);
  }
  void ensureOutput() {
    if (!field_output_options.enabled() || field_output_series)
      return;
    field_output_series = field_output_options.writer->openSeries(
        {.mesh = &transport_stepper->mesh(),
         .path = detail::plasmaFieldOutputPath(field_output_options),
         .overwrite = field_output_options.overwrite});
  }
  void finishOutput() {
    if (!field_output_series)
      return;
    if (field_output_has_snapshots) {
      const auto record = field_output_series->finish();
      TraceAdapter output_sink{&trace_sink};
      output::dump::traceCompleted(output_sink, record);
    }
    field_output_series.reset();
  }
  void writeOutput(bool terminal) {
    if (!field_output_options.enabled())
      return;
    if (detail::shouldWritePlasmaFieldSnapshot(field_output_options,
                                               clock.step(), terminal)) {
      (void)detail::writePlasmaFieldSnapshot(
          *density, *transport_stepper, field_output_options,
          *field_output_series,
          {.step = static_cast<std::uint64_t>(clock.step()),
           .time = clock.time()});
      field_output_has_snapshots = true;
    }
    if (terminal && field_output_has_snapshots)
      finishOutput();
  }
  output::OutputRecord save(const output::checkpoint::IWriter& writer,
                            const std::filesystem::path& path,
                            bool overwrite) const {
    if (state == SimulationState::failed || state == SimulationState::stopped)
      throw std::logic_error(
          "terminal simulation state cannot be checkpointed");
    const auto sources = detail::plasmaCheckpointSources(*density);
    const double previous_dt = transport_stepper->previousTimeStep();
    const std::array scalars{output::checkpoint::ScalarSource{
        .value = &previous_dt,
        .key = "adaptive_previous_time_step",
        .metadata = {.name = "adaptive_previous_time_step"}}};
    return writer.write(
        {.mesh = &transport_stepper->mesh(),
         .path = path,
         .stamp = {.step = static_cast<std::uint64_t>(clock.step()),
                   .time = clock.time()},
         .cell_fields = sources,
         .scalars = scalars,
         .overwrite = overwrite});
  }
  void complete() {
    if (state != SimulationState::completed) {
      state = SimulationState::completed;
      emitRun("run.completed", trace::Severity::Info);
    }
  }
  void fail() noexcept {
    if (state != SimulationState::failed) {
      state = SimulationState::failed;
      emitRun("run.failed", trace::Severity::Error);
    }
  }
  linalg::SolverResult perform() {
    if (clock.finished())
      throw std::out_of_range("adaptive simulation already finished");
    const std::array started{
        trace::TraceAttribute{"step", static_cast<std::uint64_t>(clock.step())},
        trace::TraceAttribute{"time", clock.time()},
        trace::TraceAttribute{"remaining_time", clock.remainingTime()},
    };
    emitTrace(trace_sink, "adaptive_step", "step.started",
              trace::Severity::Trace, started);
    const auto result = transport_stepper->prepareElectrostatics(*density);
    if (!result.success()) {
      emitDiagnostic(trace_sink, result, clock.step(), clock.time());
      emitSolverResult(trace_sink, "adaptive_step", "step.failed",
                       trace::Severity::Error, result, clock.step());
      return result;
    }
    emitSolverResult(trace_sink, "adaptive_step", "electrostatics.completed",
                     trace::Severity::Trace, result, clock.step());
    writeOutput(false);
    reaction_rates.fill(0.0);
    rate_evaluator(*density, transport_stepper->potential(),
                   transport_stepper->electricFieldNormal(), reaction_rates);
    const std::array rates{
        trace::TraceAttribute{"step", static_cast<std::uint64_t>(clock.step())},
        trace::TraceAttribute{
            "field_count", static_cast<std::uint64_t>(reaction_rates.size())},
    };
    emitTrace(trace_sink, "adaptive_step", "reaction_rates.completed",
              trace::Severity::Trace, rates);
    source.fill(0.0);
    reaction_network->accumulateSources(reaction_rates.span(), source);
    const std::array sources{
        trace::TraceAttribute{"step", static_cast<std::uint64_t>(clock.step())},
        trace::TraceAttribute{"field_count",
                              static_cast<std::uint64_t>(source.size())},
    };
    emitTrace(trace_sink, "adaptive_step", "sources.completed",
              trace::Severity::Trace, sources);
    if (statistics_options.enabled()) {
      TraceAdapter statistics_sink{&trace_sink};
      detail::emitPlasmaStatistics(
          statistics_sink, transport_stepper->species(), *density,
          transport_stepper->chargeDensity(), transport_stepper->potential(),
          transport_stepper->electricFieldNormal(), statistics_options,
          clock.step(), clock.time());
    }
    const auto proposal = transport_stepper->advancePrepared(
        *density, source, clock.remainingTime());
    const std::array selected{
        trace::TraceAttribute{"step", static_cast<std::uint64_t>(clock.step())},
        trace::TraceAttribute{"dt", proposal.dt},
        trace::TraceAttribute{"transport_limit", proposal.transport_limit},
        trace::TraceAttribute{"positivity_limit", proposal.positivity_limit},
        trace::TraceAttribute{"stability_limit", proposal.stability_limit},
    };
    emitTrace(trace_sink, "adaptive_step", "timestep.selected",
              trace::Severity::Debug, selected);
    clock.advance(proposal.dt);
    if (clock.finished() && field_output_options.enabled()) {
      const auto final_result =
          transport_stepper->prepareElectrostatics(*density);
      if (!final_result.success()) {
        emitDiagnostic(trace_sink, final_result, clock.step(), clock.time());
        emitSolverResult(trace_sink, "adaptive_step",
                         "final_output.electrostatics.failed",
                         trace::Severity::Error, final_result, clock.step());
        return final_result;
      }
      writeOutput(true);
    }
    const std::array completed{
        trace::TraceAttribute{"step", static_cast<std::uint64_t>(clock.step())},
        trace::TraceAttribute{"time", clock.time()},
        trace::TraceAttribute{"dt", proposal.dt},
        trace::TraceAttribute{"relative_residual", result.relative_residual},
    };
    emitTrace(trace_sink, "adaptive_step", "step.completed",
              trace::Severity::Info, completed);
    return result;
  }
};

AdaptiveStepPlasmaSimulation::AdaptiveStepPlasmaSimulation(
    physics::SpeciesCellFields& density,
    const physics::ReactionNetwork& reactions, Stepper& transport,
    PlasmaReactionRateEvaluator evaluator, AdaptiveTimeClock clock,
    PlasmaTraceSink sink, trace::StatisticsOptions statistics,
    FieldOutputOptions field_output, CheckpointOptions checkpoint)
    : impl_(std::make_unique<Impl>(
          density, reactions, transport, std::move(evaluator), std::move(clock),
          std::move(sink), statistics, std::move(field_output),
          std::move(checkpoint))) {}
AdaptiveStepPlasmaSimulation::~AdaptiveStepPlasmaSimulation() = default;
AdaptiveStepPlasmaSimulation::AdaptiveStepPlasmaSimulation(
    AdaptiveStepPlasmaSimulation&&) noexcept = default;
AdaptiveStepPlasmaSimulation& AdaptiveStepPlasmaSimulation::operator=(
    AdaptiveStepPlasmaSimulation&&) noexcept = default;
void AdaptiveStepPlasmaSimulation::start() {
  auto& self = *impl_;
  if (self.state == SimulationState::running)
    return;
  if (self.state != SimulationState::ready &&
      self.state != SimulationState::paused)
    throw std::logic_error("simulation cannot start from its current state");
  const bool resumed = self.state == SimulationState::paused;
  if (!resumed && !self.clock.finished())
    try {
      self.ensureOutput();
    } catch (...) {
      self.fail();
      throw;
    }
  self.state = SimulationState::running;
  self.emitRun(resumed ? "run.resumed" : "run.started", trace::Severity::Info);
  if (self.clock.finished())
    self.complete();
}
void AdaptiveStepPlasmaSimulation::pause() {
  auto& self = *impl_;
  if (self.state == SimulationState::paused)
    return;
  if (self.state != SimulationState::running)
    throw std::logic_error("only a running simulation can be paused");
  self.state = SimulationState::paused;
  self.emitRun("run.paused", trace::Severity::Info);
}
void AdaptiveStepPlasmaSimulation::stop() {
  auto& self = *impl_;
  if (self.state == SimulationState::stopped)
    return;
  if (self.state == SimulationState::completed ||
      self.state == SimulationState::failed)
    throw std::logic_error("terminal simulation cannot be stopped");
  self.state = SimulationState::stopped;
  self.finishOutput();
  self.emitRun("run.stopped", trace::Severity::Info);
}
linalg::SolverResult AdaptiveStepPlasmaSimulation::advance() {
  auto& self = *impl_;
  if (self.state == SimulationState::ready)
    start();
  if (self.state != SimulationState::running)
    throw std::logic_error("simulation step requires the running state");
  try {
    const auto result = self.perform();
    if (!result.success()) {
      self.fail();
      return result;
    }
    if (shouldWriteCheckpoint(self.checkpoint_options, self.clock.step(),
                              self.clock.finished()))
      (void)saveCheckpoint();
    if (self.clock.finished())
      self.complete();
    return result;
  } catch (...) {
    self.fail();
    throw;
  }
}
void AdaptiveStepPlasmaSimulation::run() {
  if (impl_->state == SimulationState::ready ||
      impl_->state == SimulationState::paused)
    start();
  else if (impl_->state == SimulationState::completed)
    return;
  else if (impl_->state != SimulationState::running)
    throw std::logic_error("simulation cannot run from its current state");
  while (impl_->state == SimulationState::running)
    if (!advance().success())
      throw std::runtime_error("adaptive plasma simulation failed");
}
SimulationState AdaptiveStepPlasmaSimulation::state() const noexcept {
  return impl_->state;
}
double AdaptiveStepPlasmaSimulation::time() const noexcept {
  return impl_->clock.time();
}
double AdaptiveStepPlasmaSimulation::endTime() const noexcept {
  return impl_->clock.endTime();
}
double AdaptiveStepPlasmaSimulation::remainingTime() const noexcept {
  return impl_->clock.remainingTime();
}
std::size_t AdaptiveStepPlasmaSimulation::step() const noexcept {
  return impl_->clock.step();
}
bool AdaptiveStepPlasmaSimulation::finished() const noexcept {
  return impl_->clock.finished();
}
double AdaptiveStepPlasmaSimulation::lastTimeStep() const noexcept {
  return impl_->transport_stepper->previousTimeStep();
}
bool AdaptiveStepPlasmaSimulation::hasLastTimeStepProposal() const noexcept {
  return impl_->transport_stepper->hasLastTimeStepProposal();
}
const AdaptiveStepPlasmaSimulation::TimeStepProposal&
AdaptiveStepPlasmaSimulation::lastTimeStepProposal() const {
  return impl_->transport_stepper->lastTimeStepProposal();
}
output::OutputRecord AdaptiveStepPlasmaSimulation::saveCheckpoint() const {
  if (!impl_->checkpoint_options.enabled())
    throw std::logic_error("simulation checkpoint output is disabled");
  return saveCheckpoint(*impl_->checkpoint_options.writer,
                        impl_->checkpoint_options.path,
                        impl_->checkpoint_options.overwrite);
}
output::OutputRecord AdaptiveStepPlasmaSimulation::saveCheckpoint(
    const output::checkpoint::IWriter& writer,
    const std::filesystem::path& path, bool overwrite) const {
  return impl_->save(writer, path, overwrite);
}
output::OutputRecord AdaptiveStepPlasmaSimulation::restoreCheckpoint(
    const output::checkpoint::IReader& reader,
    const std::filesystem::path& path) {
  auto& self = *impl_;
  if (self.state != SimulationState::ready)
    throw std::logic_error("checkpoint restore requires a ready simulation");
  const auto manifest = reader.inspect(path);
  if (manifest.stamp.step > std::numeric_limits<std::size_t>::max())
    throw std::overflow_error("checkpoint step is too large");
  const auto restored_step = static_cast<std::size_t>(manifest.stamp.step);
  AdaptiveTimeClock restored_clock(self.clock.endTime(), manifest.stamp.time,
                                   restored_step);
  physics::SpeciesCellFields restored_density(
      self.density->mesh(), self.density->size(), 0.0,
      (*self.density)[physics::SpeciesId{0}].metadata());
  auto targets = detail::plasmaCheckpointTargets(restored_density);
  double previous_dt{};
  const std::array scalars{output::checkpoint::ScalarTarget{
      .value = &previous_dt,
      .key = "adaptive_previous_time_step",
      .metadata = {.name = "adaptive_previous_time_step"}}};
  const auto record =
      reader.restore(path, {.mesh = &self.transport_stepper->mesh(),
                            .cell_fields = targets,
                            .scalars = scalars});
  self.transport_stepper->restoreTimeStepHistory(previous_dt);
  for (std::size_t i = 0; i < self.density->size(); ++i) {
    const physics::SpeciesId id{static_cast<std::uint32_t>(i)};
    (*self.density)[id] = std::move(restored_density[id]);
  }
  self.clock = restored_clock;
  return record;
}

}  // namespace pemu::simulation
