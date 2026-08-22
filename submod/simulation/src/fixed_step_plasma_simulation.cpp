#include <pemu/simulation/fixed_step_plasma_simulation.hpp>

#include "plasma_workflow_common.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <utility>

namespace pemu::simulation {

struct FixedStepPlasmaSimulation::Impl {
  using Core = detail::PlasmaWorkflowCore<Stepper, FixedStepClock>;

  Core core;

  Impl(physics::SpeciesCellFields& density,
       const physics::ReactionNetwork& reactions, Stepper& stepper,
       PlasmaReactionRateEvaluator evaluator, FixedStepClock clock,
       trace::AnyTraceSink sink, trace::StatisticsOptions statistics,
       FieldOutputOptions field_output, CheckpointOptions checkpoint)
      : core(density, reactions, stepper, std::move(evaluator),
             std::move(clock), std::move(sink), statistics,
             std::move(field_output), std::move(checkpoint), "fixed_step",
             detail::fixed_step_checkpoint_kind) {
    const double scale = std::max(
        {1.0, std::abs(core.clock.timeStep()), std::abs(stepper.timeStep())});
    if (std::abs(core.clock.timeStep() - stepper.timeStep()) > 1e-14 * scale) {
      throw std::invalid_argument(
          "FixedStepClock dt does not match FixedStep transport dt");
    }
  }

  [[nodiscard]] linalg::SolverResult perform() {
    auto& self = core;
    if (self.clock.finished()) {
      throw std::out_of_range("fixed-step simulation already finished");
    }
    const std::array started{
        trace::TraceAttribute{"step",
                              static_cast<std::uint64_t>(self.clock.step())},
        trace::TraceAttribute{"time", self.clock.time()},
        trace::TraceAttribute{"dt", self.clock.timeStep()},
    };
    detail::emitTrace(self.trace_sink, self.category, "step.started",
                      trace::Severity::Trace, started);
    const auto result =
        self.transport_stepper->updateElectrostatics(*self.density);
    if (!result.success()) {
      detail::emitDiagnostic(self.trace_sink, result, self.clock.step(),
                             self.clock.time());
      detail::emitSolverResult(self.trace_sink, self.category, "step.failed",
                               trace::Severity::Error, result,
                               self.clock.step());
      return result;
    }
    detail::emitSolverResult(self.trace_sink, self.category,
                             "electrostatics.completed", trace::Severity::Trace,
                             result, self.clock.step());
    self.writeOutput(false);
    self.reaction_rates.fill(0.0);
    self.rate_evaluator(*self.density, self.transport_stepper->potential(),
                        self.transport_stepper->electricFieldNormal(),
                        self.reaction_rates);
    const std::array rate_attributes{
        trace::TraceAttribute{"step",
                              static_cast<std::uint64_t>(self.clock.step())},
        trace::TraceAttribute{"field_count", static_cast<std::uint64_t>(
                                                 self.reaction_rates.size())},
    };
    detail::emitTrace(self.trace_sink, self.category,
                      "reaction_rates.completed", trace::Severity::Trace,
                      rate_attributes);
    self.source.fill(0.0);
    self.reaction_network->accumulateSources(self.reaction_rates.span(),
                                             self.source);
    const std::array source_attributes{
        trace::TraceAttribute{"step",
                              static_cast<std::uint64_t>(self.clock.step())},
        trace::TraceAttribute{"field_count",
                              static_cast<std::uint64_t>(self.source.size())},
    };
    detail::emitTrace(self.trace_sink, self.category, "sources.completed",
                      trace::Severity::Trace, source_attributes);
    if (self.statistics_options.enabled()) {
      detail::emitPlasmaStatistics(
          self.trace_sink, self.transport_stepper->species(), *self.density,
          self.transport_stepper->chargeDensity(),
          self.transport_stepper->potential(),
          self.transport_stepper->electricFieldNormal(),
          self.statistics_options, self.clock.step(), self.clock.time());
    }
    self.transport_stepper->advanceTransport(*self.density, self.source);
    self.clock.advance();
    if (self.clock.finished() && self.field_output_options.enabled()) {
      const auto final_result =
          self.transport_stepper->updateElectrostatics(*self.density);
      if (!final_result.success()) {
        detail::emitDiagnostic(self.trace_sink, final_result, self.clock.step(),
                               self.clock.time());
        detail::emitSolverResult(self.trace_sink, self.category,
                                 "final_output.electrostatics.failed",
                                 trace::Severity::Error, final_result,
                                 self.clock.step());
        return final_result;
      }
      self.writeOutput(true);
    }
    const std::array complete_attributes{
        trace::TraceAttribute{"step",
                              static_cast<std::uint64_t>(self.clock.step())},
        trace::TraceAttribute{"time", self.clock.time()},
        trace::TraceAttribute{"dt", self.clock.timeStep()},
        trace::TraceAttribute{"relative_residual", result.relative_residual},
    };
    detail::emitTrace(self.trace_sink, self.category, "step.completed",
                      trace::Severity::Info, complete_attributes);
    return result;
  }
};

FixedStepPlasmaSimulation::FixedStepPlasmaSimulation(
    physics::SpeciesCellFields& density,
    const physics::ReactionNetwork& reactions, Stepper& transport,
    PlasmaReactionRateEvaluator evaluator, FixedStepClock clock,
    trace::AnyTraceSink sink, trace::StatisticsOptions statistics,
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
  detail::startWorkflow(impl_->core);
}

void FixedStepPlasmaSimulation::pause() {
  detail::pauseWorkflow(impl_->core);
}

void FixedStepPlasmaSimulation::stop() {
  detail::stopWorkflow(impl_->core);
}

linalg::SolverResult FixedStepPlasmaSimulation::advance() {
  return detail::advanceWorkflow(
      impl_->core, [this] { return impl_->perform(); },
      [this] { (void)saveCheckpoint(); });
}

void FixedStepPlasmaSimulation::run() {
  detail::runWorkflow(
      impl_->core, [this] { return advance(); },
      "fixed-step plasma simulation failed");
}

SimulationState FixedStepPlasmaSimulation::state() const noexcept {
  return impl_->core.state;
}

double FixedStepPlasmaSimulation::time() const noexcept {
  return impl_->core.clock.time();
}

double FixedStepPlasmaSimulation::timeStep() const noexcept {
  return impl_->core.transport_stepper->timeStep();
}

std::size_t FixedStepPlasmaSimulation::step() const noexcept {
  return impl_->core.clock.step();
}

std::size_t FixedStepPlasmaSimulation::totalSteps() const noexcept {
  return impl_->core.clock.totalSteps();
}

bool FixedStepPlasmaSimulation::finished() const noexcept {
  return impl_->core.clock.finished();
}

output::OutputRecord FixedStepPlasmaSimulation::saveCheckpoint() const {
  if (!impl_->core.checkpoint_options.enabled()) {
    throw std::logic_error("simulation checkpoint output is disabled");
  }
  return saveCheckpoint(*impl_->core.checkpoint_options.writer,
                        impl_->core.checkpoint_options.path,
                        impl_->core.checkpoint_options.overwrite);
}

output::OutputRecord FixedStepPlasmaSimulation::saveCheckpoint(
    const output::checkpoint::IWriter& writer,
    const std::filesystem::path& path, bool overwrite) const {
  return impl_->core.save(writer, path, overwrite);
}

output::OutputRecord FixedStepPlasmaSimulation::restoreCheckpoint(
    const output::checkpoint::IReader& reader,
    const std::filesystem::path& path) {
  auto& self = impl_->core;
  if (self.state != SimulationState::ready) {
    throw std::logic_error("checkpoint restore requires a ready simulation");
  }
  physics::SpeciesCellFields restored_density(
      self.density->mesh(), self.density->size(), 0.0,
      (*self.density)[physics::SpeciesId{0}].metadata());
  auto targets = detail::plasmaCheckpointTargets(
      restored_density, self.transport_stepper->species());
  double schema_version{};
  double workflow_kind{};
  const auto scalars =
      detail::plasmaWorkflowScalarTargets(schema_version, workflow_kind);
  const auto record =
      reader.restore(path, {.mesh = &self.transport_stepper->mesh(),
                            .cell_fields = targets,
                            .scalars = scalars});
  detail::validatePlasmaWorkflowSchema(schema_version, workflow_kind,
                                       detail::fixed_step_checkpoint_kind);
  if (record.stamp.step > self.clock.totalSteps()) {
    throw std::invalid_argument(
        "checkpoint step exceeds fixed simulation horizon");
  }
  const auto restored_step = static_cast<std::size_t>(record.stamp.step);
  const double expected =
      static_cast<double>(restored_step) * self.clock.timeStep();
  if (!std::isfinite(record.stamp.time) ||
      std::abs(record.stamp.time - expected) >
          1e-14 * std::max(
                      {1.0, std::abs(expected), std::abs(record.stamp.time)})) {
    throw std::invalid_argument(
        "checkpoint time does not match fixed simulation step");
  }
  for (std::size_t i = 0; i < self.density->size(); ++i) {
    const physics::SpeciesId id{static_cast<std::uint32_t>(i)};
    (*self.density)[id] = std::move(restored_density[id]);
  }
  self.clock = FixedStepClock(self.clock.timeStep(), self.clock.totalSteps(),
                              restored_step);
  return record;
}

}  // namespace pemu::simulation
