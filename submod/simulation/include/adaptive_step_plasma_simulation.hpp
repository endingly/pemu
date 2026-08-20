#pragma once

#include <pemu/equation/adaptive_step_multi_species_drift_diffusion_stepper.hpp>
#include <pemu/equation/time_integration/adaptive_time_step_controller.hpp>
#include <pemu/trace/trace.hpp>
#include <pemu/physics/reaction.hpp>
#include <pemu/physics/species.hpp>
#include <pemu/simulation/adaptive_time_clock.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>

namespace pemu::simulation {

template <typename ReactionRateEvaluator,
          pemu::trace::TraceSink TraceSink = pemu::trace::NullTraceSink>
class AdaptiveStepPlasmaSimulation {
 public:
  using Stepper = equation::AdaptiveStepMultiSpeciesDriftDiffusionStepper;

  using TimeStepProposal = equation::time_integration::TimeStepProposal;

  AdaptiveStepPlasmaSimulation(physics::SpeciesCellFields& density,

                               const physics::ReactionNetwork& reaction_network,

                               Stepper& transport_stepper,

                               ReactionRateEvaluator rate_evaluator,

                               AdaptiveTimeClock clock,

                               TraceSink trace_sink = {})

      : density_(&density),

        reaction_network_(&reaction_network),

        transport_stepper_(&transport_stepper),

        rate_evaluator_(std::move(rate_evaluator)),

        clock_(std::move(clock)),

        trace_sink_(std::move(trace_sink)),

        reaction_rates_(density.mesh(), reaction_network.size(), 0.0,
                        transport_stepper.fieldMetadata().reaction_rate),

        source_(density.mesh(), density.size(), 0.0,
                transport_stepper.fieldMetadata().number_density_source) {
    if (density_->size() == 0) {

      throw std::invalid_argument(
          "simulation density must "
          "contain species");
    }
  }

  AdaptiveStepPlasmaSimulation(const AdaptiveStepPlasmaSimulation&) = delete;

  AdaptiveStepPlasmaSimulation& operator=(const AdaptiveStepPlasmaSimulation&) =
      delete;

  // ========================================================
  // Execute one complete adaptive timestep:
  //
  // n^k
  //   │
  //   ▼
  // electrostatics
  //   │
  //   ▼
  // R^k
  //   │
  //   ▼
  // S^k
  //   │
  //   ▼
  // determine dt_k
  //   │
  //   ▼
  // n^(k+1)
  //   │
  //   ▼
  // t_(k+1) = t_k + dt_k
  // ========================================================

  [[nodiscard]]
  linalg::SolverResult advanceOneStep() {
    if (clock_.finished()) {

      throw std::out_of_range(
          "adaptive simulation "
          "already finished");
    }

    if constexpr (tracing_enabled_) {
      emitStepStarted();
    }

    // ----------------------------------------------------
    // 1.
    //
    // n^k -> rho^k -> phi^k -> E^k
    // ----------------------------------------------------

    const auto result = transport_stepper_->prepareElectrostatics(*density_);

    if (!result.success()) {
      if constexpr (tracing_enabled_) {
        emitSolverResult("step.failed", pemu::trace::Severity::Error, result);
      }
      return result;
    }

    if constexpr (tracing_enabled_) {
      emitSolverResult("electrostatics.completed", pemu::trace::Severity::Trace,
                       result);
    }

    // ----------------------------------------------------
    // 2.
    //
    // Evaluate reaction rate using:
    //
    //     n^k
    //     phi^k
    //     E^k
    // ----------------------------------------------------

    reaction_rates_.fill(0.0);

    rate_evaluator_(*density_,

                    transport_stepper_->potential(),

                    transport_stepper_->electricFieldNormal(),

                    reaction_rates_);

    if constexpr (tracing_enabled_) {
      emitWorkspaceCompleted("reaction_rates.completed",
                             reaction_rates_.size());
    }

    // ----------------------------------------------------
    // 3.
    //
    // R -> S
    // ----------------------------------------------------

    source_.fill(0.0);

    reaction_network_->accumulateSources(reaction_rates_.span(), source_);

    if constexpr (tracing_enabled_) {
      emitWorkspaceCompleted("sources.completed", source_.size());
    }

    // ----------------------------------------------------
    // 4.
    //
    // Adaptive stepper:
    //
    // transport limit
    // reaction/positivity limit
    // safety
    // growth limiter
    // remaining time
    //
    // -> dt_k
    //
    // and then:
    //
    // n^k -> n^(k+1)
    // ----------------------------------------------------

    const auto proposal = transport_stepper_->advancePrepared(
        *density_, source_, clock_.remainingTime());

    if constexpr (tracing_enabled_) {
      emitTimeStepSelected(proposal);
    }

    // ----------------------------------------------------
    // 5.
    //
    // Commit time AFTER physical state update succeeds.
    // ----------------------------------------------------

    clock_.advance(proposal.dt);

    if constexpr (tracing_enabled_) {
      emitStepCompleted(proposal.dt, result);
    }

    return result;
  }

  // ========================================================
  // Run until t_end.
  // ========================================================

  void run() {
    if constexpr (tracing_enabled_) {
      emitRunEvent("run.started", pemu::trace::Severity::Info);
    }

    while (!clock_.finished()) {

      const auto result = advanceOneStep();

      if (!result.success()) {

        if constexpr (tracing_enabled_) {
          emitRunEvent("run.failed", pemu::trace::Severity::Error);
        }

        throw std::runtime_error(
            "adaptive plasma "
            "simulation failed");
      }
    }

    if constexpr (tracing_enabled_) {
      emitRunEvent("run.completed", pemu::trace::Severity::Info);
    }
  }

  // ========================================================
  // Time
  // ========================================================

  [[nodiscard]]
  double time() const noexcept {
    return clock_.time();
  }

  [[nodiscard]]
  double endTime() const noexcept {
    return clock_.endTime();
  }

  [[nodiscard]]
  double remainingTime() const noexcept {
    return clock_.remainingTime();
  }

  [[nodiscard]]
  std::size_t step() const noexcept {
    return clock_.step();
  }

  [[nodiscard]]
  bool finished() const noexcept {
    return clock_.finished();
  }

  // ========================================================
  // Adaptive timestep diagnostics
  // ========================================================

  [[nodiscard]]
  double lastTimeStep() const noexcept {
    return transport_stepper_->previousTimeStep();
  }

  [[nodiscard]]
  bool hasLastTimeStepProposal() const noexcept {
    return transport_stepper_->hasLastTimeStepProposal();
  }

  [[nodiscard]]
  const TimeStepProposal& lastTimeStepProposal() const {
    return transport_stepper_->lastTimeStepProposal();
  }

  // ========================================================
  // State
  // ========================================================

  [[nodiscard]]
  physics::SpeciesCellFields& density() noexcept {
    return *density_;
  }

  [[nodiscard]]
  const physics::SpeciesCellFields& density() const noexcept {
    return *density_;
  }

  // ========================================================
  // Reaction workspace
  // ========================================================

  [[nodiscard]]
  const physics::ReactionRateFields& reactionRates() const noexcept {
    return reaction_rates_;
  }

  [[nodiscard]]
  const physics::SpeciesCellFields& source() const noexcept {
    return source_;
  }

  // ========================================================
  // Electrostatic diagnostics
  // ========================================================

  [[nodiscard]]
  const field::CellField<double>& chargeDensity() const noexcept {
    return transport_stepper_->chargeDensity();
  }

  [[nodiscard]]
  const field::CellField<double>& potential() const noexcept {
    return transport_stepper_->potential();
  }

  [[nodiscard]]
  const field::FaceField<double>& electricFieldNormal() const noexcept {
    return transport_stepper_->electricFieldNormal();
  }

  [[nodiscard]]
  const field::FaceField<double>& driftVelocityNormal(
      physics::SpeciesId id) const {
    return transport_stepper_->driftVelocityNormal(id);
  }

 private:
  static constexpr bool tracing_enabled_ =
      !std::same_as<std::remove_cvref_t<TraceSink>, pemu::trace::NullTraceSink>;

  void emitRunEvent(std::string_view name,
                    pemu::trace::Severity severity) noexcept {
    const std::array attributes{
        pemu::trace::TraceAttribute{"step",
                                  static_cast<std::uint64_t>(clock_.step())},
        pemu::trace::TraceAttribute{"time", clock_.time()},
        pemu::trace::TraceAttribute{"end_time", clock_.endTime()},
    };
    emitTrace(name, severity, attributes);
  }

  void emitStepStarted() noexcept {
    const std::array attributes{
        pemu::trace::TraceAttribute{"step",
                                  static_cast<std::uint64_t>(clock_.step())},
        pemu::trace::TraceAttribute{"time", clock_.time()},
        pemu::trace::TraceAttribute{"remaining_time", clock_.remainingTime()},
    };
    emitTrace("step.started", pemu::trace::Severity::Trace, attributes);
  }

  void emitSolverResult(std::string_view name, pemu::trace::Severity severity,
                        const linalg::SolverResult& result) noexcept {
    const std::array attributes{
        pemu::trace::TraceAttribute{"step",
                                  static_cast<std::uint64_t>(clock_.step())},
        pemu::trace::TraceAttribute{"solver_status",
                                  static_cast<std::int64_t>(result.status)},
        pemu::trace::TraceAttribute{"residual_norm", result.residual_norm},
        pemu::trace::TraceAttribute{"relative_residual",
                                  result.relative_residual},
    };
    emitTrace(name, severity, attributes);
  }

  void emitWorkspaceCompleted(std::string_view name,
                              std::size_t field_count) noexcept {
    const std::array attributes{
        pemu::trace::TraceAttribute{"step",
                                  static_cast<std::uint64_t>(clock_.step())},
        pemu::trace::TraceAttribute{"field_count",
                                  static_cast<std::uint64_t>(field_count)},
    };
    emitTrace(name, pemu::trace::Severity::Trace, attributes);
  }

  void emitTimeStepSelected(const TimeStepProposal& proposal) noexcept {
    const std::array attributes{
        pemu::trace::TraceAttribute{"step",
                                  static_cast<std::uint64_t>(clock_.step())},
        pemu::trace::TraceAttribute{"dt", proposal.dt},
        pemu::trace::TraceAttribute{"transport_limit", proposal.transport_limit},
        pemu::trace::TraceAttribute{"positivity_limit",
                                  proposal.positivity_limit},
        pemu::trace::TraceAttribute{"stability_limit", proposal.stability_limit},
    };
    emitTrace("timestep.selected", pemu::trace::Severity::Debug, attributes);
  }

  void emitStepCompleted(double dt,
                         const linalg::SolverResult& result) noexcept {
    const std::array attributes{
        pemu::trace::TraceAttribute{"step",
                                  static_cast<std::uint64_t>(clock_.step())},
        pemu::trace::TraceAttribute{"time", clock_.time()},
        pemu::trace::TraceAttribute{"dt", dt},
        pemu::trace::TraceAttribute{"relative_residual",
                                  result.relative_residual},
    };
    emitTrace("step.completed", pemu::trace::Severity::Info, attributes);
  }

  template <std::size_t N>
  void emitTrace(
      std::string_view name, pemu::trace::Severity severity,
      const std::array<pemu::trace::TraceAttribute, N>& attributes) noexcept {
    trace_sink_({.category = "simulation.adaptive_step",
                 .name = name,
                 .severity = severity,
                 .attributes = attributes});
  }

  // --------------------------------------------------------
  // External state
  // --------------------------------------------------------

  physics::SpeciesCellFields* density_;

  const physics::ReactionNetwork* reaction_network_;

  Stepper* transport_stepper_;

  // --------------------------------------------------------
  // Physics evaluator
  // --------------------------------------------------------

  ReactionRateEvaluator rate_evaluator_;

  // --------------------------------------------------------
  // Adaptive simulation clock
  // --------------------------------------------------------

  AdaptiveTimeClock clock_;

  [[no_unique_address]] TraceSink trace_sink_;

  // --------------------------------------------------------
  // Workspaces
  // --------------------------------------------------------

  physics::ReactionRateFields reaction_rates_;

  physics::SpeciesCellFields source_;
};

}  // namespace pemu::simulation
