#pragma once

#include <pemu/output/dump/trace.hpp>
#include <pemu/simulation/checkpoint.hpp>
#include <pemu/simulation/detail/plasma_checkpoint.hpp>
#include <pemu/simulation/detail/plasma_field_output.hpp>
#include <pemu/simulation/detail/plasma_statistics.hpp>
#include <pemu/simulation/plasma_reaction_rate_evaluator.hpp>
#include <pemu/simulation/simulation_state.hpp>
#include <pemu/trace/any_trace_sink.hpp>

#include <array>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace pemu::simulation::detail {

inline constexpr std::size_t trace_flush_step_interval = 2;

template <std::size_t N>
void emitTrace(
    const trace::AnyTraceSink& sink, std::string_view category,
    std::string_view name, trace::Severity severity,
    const std::array<trace::TraceAttribute, N>& attributes) noexcept {
  sink({.domain = trace::DiagDomain::simulation,
        .category = category,
        .name = name,
        .severity = severity,
        .attributes = attributes});
}

inline void emitDiagnostic(const trace::AnyTraceSink& sink,
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
  sink(diagnostic);
}

inline void emitSolverResult(const trace::AnyTraceSink& sink,
                             std::string_view category, std::string_view name,
                             trace::Severity severity,
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

inline std::array<output::checkpoint::ScalarTarget, 2>
plasmaWorkflowScalarTargets(double& schema_version, double& workflow_kind) {
  return {
      output::checkpoint::ScalarTarget{
          .value = &schema_version,
          .key = std::string{plasma_checkpoint_schema_key},
          .metadata = {.name = std::string{plasma_checkpoint_schema_key}}},
      output::checkpoint::ScalarTarget{
          .value = &workflow_kind,
          .key = std::string{plasma_checkpoint_kind_key},
          .metadata = {.name = std::string{plasma_checkpoint_kind_key}}},
  };
}

inline void validatePlasmaWorkflowSchema(double schema_version,
                                         double workflow_kind,
                                         double expected_kind) {
  if (schema_version != plasma_checkpoint_schema_version ||
      workflow_kind != expected_kind) {
    throw std::invalid_argument("checkpoint workflow schema differs");
  }
}

template <typename Stepper, typename Clock>
struct PlasmaWorkflowCore {
  physics::SpeciesCellFields* density;
  const physics::ReactionNetwork* reaction_network;
  Stepper* transport_stepper;
  PlasmaReactionRateEvaluator rate_evaluator;
  Clock clock;
  trace::AnyTraceSink trace_sink;
  trace::StatisticsOptions statistics_options;
  FieldOutputOptions field_output_options;
  CheckpointOptions checkpoint_options;
  physics::ReactionRateFields reaction_rates;
  physics::SpeciesCellFields source;
  std::unique_ptr<output::dump::ISeries> field_output_series;
  bool field_output_has_snapshots{};
  SimulationState state{SimulationState::ready};
  std::string_view category;
  double workflow_kind;

  PlasmaWorkflowCore(physics::SpeciesCellFields& density_in,
                     const physics::ReactionNetwork& reactions,
                     Stepper& stepper, PlasmaReactionRateEvaluator evaluator,
                     Clock clock_in, trace::AnyTraceSink sink,
                     trace::StatisticsOptions statistics,
                     FieldOutputOptions field_output,
                     CheckpointOptions checkpoint,
                     std::string_view trace_category, double checkpoint_kind)
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
               stepper.fieldMetadata().number_density_source),
        category(trace_category),
        workflow_kind(checkpoint_kind) {
    if (!rate_evaluator) {
      throw std::invalid_argument(
          "simulation reaction evaluator must not be empty");
    }
    if (density->size() == 0) {
      throw std::invalid_argument("simulation density must contain species");
    }
    if (density->size() != stepper.species().size()) {
      throw std::invalid_argument(
          "simulation density and transport species counts differ");
    }
    validatePlasmaStatisticsConfiguration(
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
    emitTrace(trace_sink, category, name, severity, attributes);
  }

  void ensureOutput() {
    if (!field_output_options.enabled() || field_output_series) {
      return;
    }
    field_output_series = field_output_options.writer->openSeries(
        {.mesh = &transport_stepper->mesh(),
         .path = plasmaFieldOutputPath(field_output_options),
         .overwrite = field_output_options.overwrite});
  }

  void finishOutput() {
    if (!field_output_series) {
      return;
    }
    if (field_output_has_snapshots) {
      const auto record = field_output_series->finish();
      output::dump::traceCompleted(trace_sink, record);
    }
    field_output_series.reset();
  }

  void writeOutput(bool terminal) {
    if (!field_output_options.enabled()) {
      return;
    }
    if (shouldWritePlasmaFieldSnapshot(field_output_options, clock.step(),
                                       terminal)) {
      (void)writePlasmaFieldSnapshot(
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

  [[nodiscard]] output::OutputRecord save(
      const output::checkpoint::IWriter& writer,
      const std::filesystem::path& path, bool overwrite,
      std::span<const output::checkpoint::ScalarSource> extra_scalars = {})
      const {
    if (state == SimulationState::failed || state == SimulationState::stopped) {
      throw std::logic_error(
          "terminal simulation state cannot be checkpointed");
    }
    const auto sources =
        plasmaCheckpointSources(*density, transport_stepper->species());
    const double schema_version = plasma_checkpoint_schema_version;
    std::vector<output::checkpoint::ScalarSource> scalars{
        {.value = &schema_version,
         .key = std::string{plasma_checkpoint_schema_key},
         .metadata = {.name = std::string{plasma_checkpoint_schema_key}}},
        {.value = &workflow_kind,
         .key = std::string{plasma_checkpoint_kind_key},
         .metadata = {.name = std::string{plasma_checkpoint_kind_key}}},
    };
    scalars.insert(scalars.end(), extra_scalars.begin(), extra_scalars.end());
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
      trace_sink.flush();
    }
  }

  void fail() noexcept {
    if (state != SimulationState::failed) {
      state = SimulationState::failed;
      emitRun("run.failed", trace::Severity::Error);
      trace_sink.flush();
    }
  }
};

template <typename Core>
void startWorkflow(Core& core) {
  if (core.state == SimulationState::running) {
    return;
  }
  if (core.state != SimulationState::ready &&
      core.state != SimulationState::paused) {
    throw std::logic_error("simulation cannot start from its current state");
  }
  const bool resumed = core.state == SimulationState::paused;
  if (!resumed && !core.clock.finished()) {
    try {
      core.ensureOutput();
    } catch (...) {
      core.fail();
      throw;
    }
  }
  core.state = SimulationState::running;
  core.emitRun(resumed ? "run.resumed" : "run.started", trace::Severity::Info);
  if (core.clock.finished()) {
    core.complete();
  }
}

template <typename Core>
void pauseWorkflow(Core& core) {
  if (core.state == SimulationState::paused) {
    return;
  }
  if (core.state != SimulationState::running) {
    throw std::logic_error("only a running simulation can be paused");
  }
  core.state = SimulationState::paused;
  core.emitRun("run.paused", trace::Severity::Info);
  core.trace_sink.flush();
}

template <typename Core>
void stopWorkflow(Core& core) {
  if (core.state == SimulationState::stopped) {
    return;
  }
  if (isTerminal(core.state)) {
    throw std::logic_error("terminal simulation cannot be stopped");
  }
  try {
    core.finishOutput();
  } catch (...) {
    core.fail();
    throw;
  }
  core.state = SimulationState::stopped;
  core.emitRun("run.stopped", trace::Severity::Info);
  core.trace_sink.flush();
}

template <typename Core, typename Perform, typename SaveCheckpoint>
[[nodiscard]] linalg::SolverResult advanceWorkflow(
    Core& core, Perform&& perform, SaveCheckpoint&& save_checkpoint) {
  if (core.state == SimulationState::ready) {
    startWorkflow(core);
  }
  if (core.state != SimulationState::running) {
    throw std::logic_error("simulation step requires the running state");
  }
  try {
    const auto result = std::invoke(std::forward<Perform>(perform));
    if (!result.success()) {
      core.fail();
      return result;
    }
    // Flush cadence is a workflow policy. Keeping it here prevents generic
    // sinks and channel routers from parsing simulation-specific event names.
    if (!core.clock.finished() &&
        core.clock.step() % trace_flush_step_interval == 0) {
      core.trace_sink.flush();
    }
    if (shouldWriteCheckpoint(core.checkpoint_options, core.clock.step(),
                              core.clock.finished())) {
      std::invoke(std::forward<SaveCheckpoint>(save_checkpoint));
    }
    if (core.clock.finished()) {
      core.complete();
    }
    return result;
  } catch (...) {
    core.fail();
    throw;
  }
}

template <typename Core, typename Advance>
void runWorkflow(Core& core, Advance&& advance,
                 std::string_view failure_message) {
  if (core.state == SimulationState::ready ||
      core.state == SimulationState::paused) {
    startWorkflow(core);
  } else if (core.state == SimulationState::completed) {
    return;
  } else if (core.state != SimulationState::running) {
    throw std::logic_error("simulation cannot run from its current state");
  }
  while (core.state == SimulationState::running) {
    if (!std::invoke(advance).success()) {
      throw std::runtime_error(std::string{failure_message});
    }
  }
}

}  // namespace pemu::simulation::detail
