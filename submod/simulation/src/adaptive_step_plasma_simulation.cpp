#include <pemu/simulation/adaptive_step_plasma_simulation.hpp>

#include "plasma_workflow_common.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>

namespace pemu::simulation {

struct AdaptiveStepPlasmaSimulation::Impl {
  using Core = detail::PlasmaWorkflowCore<Stepper, AdaptiveTimeClock>;

  Core core;

  Impl(physics::SpeciesCellFields& density,
       const physics::reaction::ReactionNetwork& reactions, Stepper& stepper,
       PlasmaReactionRateEvaluator evaluator,
       ElectronEnergyConfiguration electron_energy, AdaptiveTimeClock clock,
       trace::AnyTraceSink sink, trace::StatisticsOptions statistics,
       FieldOutputOptions field_output, CheckpointOptions checkpoint)
      : core(density, reactions, stepper, std::move(evaluator),
             std::move(electron_energy),
             std::move(clock), std::move(sink), statistics,
             std::move(field_output), std::move(checkpoint), "adaptive_step",
             detail::adaptive_step_checkpoint_kind) {}

  [[nodiscard]] linalg::SolverResult perform() {
    auto& self = core;
    if (self.clock.finished()) {
      throw std::out_of_range("adaptive simulation already finished");
    }
    const std::array started{
        trace::TraceAttribute{"step",
                              static_cast<std::uint64_t>(self.clock.step())},
        trace::TraceAttribute{"time", self.clock.time()},
        trace::TraceAttribute{"remaining_time", self.clock.remainingTime()},
    };
    detail::emitTrace(self.trace_sink, self.category, "step.started",
                      trace::Severity::Trace, started);
    const auto result =
        self.transport_stepper->prepareElectrostatics(*self.density);
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
    self.electron_energy.refreshMeanEnergy(*self.density);
    self.writeOutput(false);
    self.reaction_rates.fill(0.0);
    self.rate_evaluator(
        {.density = *self.density,
         .potential = self.transport_stepper->potential(),
         .electric_field_normal =
             self.transport_stepper->electricFieldNormal(),
         .electron_mean_energy = self.electron_energy.mean_energy},
        self.reaction_rates);
    const std::array rates{
        trace::TraceAttribute{"step",
                              static_cast<std::uint64_t>(self.clock.step())},
        trace::TraceAttribute{"field_count", static_cast<std::uint64_t>(
                                                 self.reaction_rates.size())},
    };
    detail::emitTrace(self.trace_sink, self.category,
                      "reaction_rates.completed", trace::Severity::Trace,
                      rates);
    self.source.fill(0.0);
    self.reaction_network->accumulateSources(self.reaction_rates.span(),
                                             self.source);
    const std::array sources{
        trace::TraceAttribute{"step",
                              static_cast<std::uint64_t>(self.clock.step())},
        trace::TraceAttribute{"field_count",
                              static_cast<std::uint64_t>(self.source.size())},
    };
    detail::emitTrace(self.trace_sink, self.category, "sources.completed",
                      trace::Severity::Trace, sources);
    self.evaluateElectronEnergySource();
    if (self.statistics_options.enabled()) {
      detail::emitPlasmaStatistics(
          self.trace_sink, self.transport_stepper->species(), *self.density,
          self.transport_stepper->chargeDensity(),
          self.transport_stepper->potential(),
          self.transport_stepper->electricFieldNormal(),
          *self.electron_energy.energy_density,
          self.electron_energy.mean_energy, self.electron_energy.source,
          self.statistics_options, self.clock.step(), self.clock.time());
    }
    const auto proposal = self.transport_stepper->proposeTimeStep(
        *self.density, self.source, self.clock.remainingTime(),
        self.electron_energy.maxStableTransportTimeStep(),
        self.electron_energy.maxPositiveTimeStep());
    const std::array selected{
        trace::TraceAttribute{"step",
                              static_cast<std::uint64_t>(self.clock.step())},
        trace::TraceAttribute{"dt", proposal.dt},
        trace::TraceAttribute{"transport_limit", proposal.transport_limit},
        trace::TraceAttribute{"positivity_limit", proposal.positivity_limit},
        trace::TraceAttribute{"stability_limit", proposal.stability_limit},
    };
    detail::emitTrace(self.trace_sink, self.category, "timestep.selected",
                      trace::Severity::Debug, selected);
    self.electron_energy.prepareIncrement(proposal.dt);
    self.transport_stepper->advancePrepared(*self.density, self.source,
                                            proposal);
    self.commitElectronEnergy(proposal.dt);
    self.clock.advance(proposal.dt);
    if (self.clock.finished() && self.field_output_options.enabled()) {
      const auto final_result =
          self.transport_stepper->prepareElectrostatics(*self.density);
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
    const std::array completed{
        trace::TraceAttribute{"step",
                              static_cast<std::uint64_t>(self.clock.step())},
        trace::TraceAttribute{"time", self.clock.time()},
        trace::TraceAttribute{"dt", proposal.dt},
        trace::TraceAttribute{"relative_residual", result.relative_residual},
    };
    detail::emitTrace(self.trace_sink, self.category, "step.completed",
                      trace::Severity::Info, completed);
    return result;
  }
};

AdaptiveStepPlasmaSimulation::AdaptiveStepPlasmaSimulation(
    physics::SpeciesCellFields& density,
    const physics::reaction::ReactionNetwork& reactions, Stepper& transport,
    PlasmaReactionRateEvaluator evaluator,
    ElectronEnergyConfiguration electron_energy, AdaptiveTimeClock clock,
    trace::AnyTraceSink sink, trace::StatisticsOptions statistics,
    FieldOutputOptions field_output, CheckpointOptions checkpoint)
    : impl_(std::make_unique<Impl>(
          density, reactions, transport, std::move(evaluator),
          std::move(electron_energy), std::move(clock),
          std::move(sink), statistics, std::move(field_output),
          std::move(checkpoint))) {}

AdaptiveStepPlasmaSimulation::~AdaptiveStepPlasmaSimulation() = default;
AdaptiveStepPlasmaSimulation::AdaptiveStepPlasmaSimulation(
    AdaptiveStepPlasmaSimulation&&) noexcept = default;
AdaptiveStepPlasmaSimulation& AdaptiveStepPlasmaSimulation::operator=(
    AdaptiveStepPlasmaSimulation&&) noexcept = default;

void AdaptiveStepPlasmaSimulation::start() {
  detail::startWorkflow(impl_->core);
}

void AdaptiveStepPlasmaSimulation::pause() {
  detail::pauseWorkflow(impl_->core);
}

void AdaptiveStepPlasmaSimulation::stop() {
  detail::stopWorkflow(impl_->core);
}

linalg::SolverResult AdaptiveStepPlasmaSimulation::advance() {
  return detail::advanceWorkflow(
      impl_->core, [this] { return impl_->perform(); },
      [this] { (void)saveCheckpoint(); });
}

void AdaptiveStepPlasmaSimulation::run() {
  detail::runWorkflow(
      impl_->core, [this] { return advance(); },
      "adaptive plasma simulation failed");
}

SimulationState AdaptiveStepPlasmaSimulation::state() const noexcept {
  return impl_->core.state;
}

double AdaptiveStepPlasmaSimulation::time() const noexcept {
  return impl_->core.clock.time();
}

double AdaptiveStepPlasmaSimulation::endTime() const noexcept {
  return impl_->core.clock.endTime();
}

double AdaptiveStepPlasmaSimulation::remainingTime() const noexcept {
  return impl_->core.clock.remainingTime();
}

std::size_t AdaptiveStepPlasmaSimulation::step() const noexcept {
  return impl_->core.clock.step();
}

bool AdaptiveStepPlasmaSimulation::finished() const noexcept {
  return impl_->core.clock.finished();
}

double AdaptiveStepPlasmaSimulation::lastTimeStep() const noexcept {
  return impl_->core.transport_stepper->previousTimeStep();
}

bool AdaptiveStepPlasmaSimulation::hasLastTimeStepProposal() const noexcept {
  return impl_->core.transport_stepper->hasLastTimeStepProposal();
}

const AdaptiveStepPlasmaSimulation::TimeStepProposal&
AdaptiveStepPlasmaSimulation::lastTimeStepProposal() const {
  return impl_->core.transport_stepper->lastTimeStepProposal();
}

const field::CellField<double>&
AdaptiveStepPlasmaSimulation::electronEnergyDensity() const noexcept {
  return *impl_->core.electron_energy.energy_density;
}

const field::CellField<double>&
AdaptiveStepPlasmaSimulation::electronMeanEnergy() const noexcept {
  return impl_->core.electron_energy.mean_energy;
}

const field::CellField<double>&
AdaptiveStepPlasmaSimulation::electronEnergySource() const noexcept {
  return impl_->core.electron_energy.source;
}

output::OutputRecord AdaptiveStepPlasmaSimulation::saveCheckpoint() const {
  if (!impl_->core.checkpoint_options.enabled()) {
    throw std::logic_error("simulation checkpoint output is disabled");
  }
  return saveCheckpoint(*impl_->core.checkpoint_options.writer,
                        impl_->core.checkpoint_options.path,
                        impl_->core.checkpoint_options.overwrite);
}

output::OutputRecord AdaptiveStepPlasmaSimulation::saveCheckpoint(
    const output::checkpoint::IWriter& writer,
    const std::filesystem::path& path, bool overwrite) const {
  const double previous_dt = impl_->core.transport_stepper->previousTimeStep();
  const std::array extra_scalars{output::checkpoint::ScalarSource{
      .value = &previous_dt,
      .key = "adaptive_previous_time_step",
      .metadata = {.name = "adaptive_previous_time_step"}}};
  return impl_->core.save(writer, path, overwrite, extra_scalars);
}

output::OutputRecord AdaptiveStepPlasmaSimulation::restoreCheckpoint(
    const output::checkpoint::IReader& reader,
    const std::filesystem::path& path) {
  auto& self = impl_->core;
  if (self.state != SimulationState::ready) {
    throw std::logic_error("checkpoint restore requires a ready simulation");
  }
  physics::SpeciesCellFields restored_density(
      self.density->mesh(), self.density->size(), 0.0,
      (*self.density)[physics::SpeciesId{0}].metadata());
  field::CellField<double> restored_energy(
      self.density->mesh(), 0.0,
      self.electron_energy.energy_density->metadata());
  auto targets = detail::plasmaCheckpointTargets(
      restored_density, self.transport_stepper->species());
  targets.push_back(
      {.field = &restored_energy,
       .key = detail::electronEnergyDensityCheckpointKey(
           self.electron_energy.electron,
           self.transport_stepper->species()
               .at(self.electron_energy.electron)
               .name)});
  double schema_version{};
  double workflow_kind{};
  double previous_dt{};
  const auto workflow_scalars =
      detail::plasmaWorkflowScalarTargets(schema_version, workflow_kind);
  const std::array scalars{
      workflow_scalars[0],
      workflow_scalars[1],
      output::checkpoint::ScalarTarget{
          .value = &previous_dt,
          .key = "adaptive_previous_time_step",
          .metadata = {.name = "adaptive_previous_time_step"}},
  };
  const auto record =
      reader.restore(path, {.mesh = &self.transport_stepper->mesh(),
                            .cell_fields = targets,
                            .scalars = scalars});
  detail::validatePlasmaWorkflowSchema(schema_version, workflow_kind,
                                       detail::adaptive_step_checkpoint_kind);
  if (record.stamp.step > std::numeric_limits<std::size_t>::max()) {
    throw std::overflow_error("checkpoint step is too large");
  }
  const auto restored_step = static_cast<std::size_t>(record.stamp.step);
  AdaptiveTimeClock restored_clock(self.clock.endTime(), record.stamp.time,
                                   restored_step);
  field::CellField<double> restored_mean_energy(
      self.density->mesh(), 0.0,
      self.transport_stepper->fieldMetadata().electron_mean_energy);
  physics::computeElectronMeanEnergy(
      restored_energy, restored_density[self.electron_energy.electron],
      self.electron_energy.density_floor, restored_mean_energy);
  self.transport_stepper->restoreTimeStepHistory(previous_dt);
  for (std::size_t i = 0; i < self.density->size(); ++i) {
    const physics::SpeciesId id{static_cast<std::uint32_t>(i)};
    (*self.density)[id] = std::move(restored_density[id]);
  }
  *self.electron_energy.energy_density = std::move(restored_energy);
  self.electron_energy.mean_energy = std::move(restored_mean_energy);
  self.clock = restored_clock;
  return record;
}

}  // namespace pemu::simulation
