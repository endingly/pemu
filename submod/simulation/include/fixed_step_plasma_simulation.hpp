#pragma once

#include <pemu/equation/fixed_step_multi_species_drift_diffusion_stepper.hpp>
#include <pemu/trace/trace.hpp>
#include <pemu/physics/reaction.hpp>
#include <pemu/physics/species.hpp>
#include <pemu/simulation/fixed_step_clock.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>

namespace pemu::simulation {

template <typename ReactionRateEvaluator,
          pemu::trace::TraceSink TraceSink = pemu::trace::NullTraceSink>
class FixedStepPlasmaSimulation {
 public:
  using Stepper = equation::FixedStepMultiSpeciesDriftDiffusionStepper;

  FixedStepPlasmaSimulation(physics::SpeciesCellFields& density,

                            const physics::ReactionNetwork& reaction_network,

                            Stepper& transport_stepper,

                            ReactionRateEvaluator rate_evaluator,

                            FixedStepClock clock,

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
    validateConfiguration();
  }

  FixedStepPlasmaSimulation(const FixedStepPlasmaSimulation&) = delete;

  FixedStepPlasmaSimulation& operator=(const FixedStepPlasmaSimulation&) =
      delete;

  // ========================================================
  // Execute one complete fixed-dt plasma timestep:
  //
  // n^k
  //   │
  //   ▼
  // rho^k
  //   │
  //   ▼
  // Poisson
  //   │
  //   ▼
  // phi^k, E^k
  //   │
  //   ▼
  // reaction rate R^k(n^k, phi^k, E^k)
  //   │
  //   ▼
  // S^k = nu R^k
  //   │
  //   ▼
  // fixed-dt SG transport
  //   │
  //   ▼
  // n^(k+1)
  // ========================================================

  [[nodiscard]]
  linalg::SolverResult advanceOneStep() {
    if (clock_.finished()) {

      throw std::out_of_range(
          "fixed-step simulation "
          "already finished");
    }

    if constexpr (tracing_enabled_) {
      emitStepStarted();
    }

    // ----------------------------------------------------
    // 1.
    //
    // n^k
    //
    // -> rho^k
    // -> phi^k
    // -> E^k
    // -> v_s^k
    // ----------------------------------------------------

    const auto result = transport_stepper_->updateElectrostatics(*density_);

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
    // Evaluate reaction rates from exactly the same
    // time-level:
    //
    //     n^k
    //     phi^k
    //     E^k
    //
    // The evaluator overwrites/fills reaction_rates_.
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
    // Reaction rates -> species source:
    //
    //     S_s = sum_r nu_sr R_r
    // ----------------------------------------------------

    source_.fill(0.0);

    reaction_network_->accumulateSources(reaction_rates_.span(), source_);

    if constexpr (tracing_enabled_) {
      emitWorkspaceCompleted("sources.completed", source_.size());
    }

    // ----------------------------------------------------
    // 4.
    //
    // Advance:
    //
    //     n^k -> n^(k+1)
    //
    // using the FIXED timestep stored in the stepper.
    //
    // Electrostatic fields remain frozen at time level k.
    // ----------------------------------------------------

    transport_stepper_->advanceTransport(*density_, source_);

    // ----------------------------------------------------
    // 5.
    //
    // Commit simulation time only after the physical state
    // has been updated successfully.
    // ----------------------------------------------------

    clock_.advance();

    if constexpr (tracing_enabled_) {
      emitStepCompleted(result);
    }

    return result;
  }

  // ========================================================
  // Run until FixedStepClock is finished.
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
            "fixed-step plasma "
            "simulation failed");
      }
    }

    if constexpr (tracing_enabled_) {
      emitRunEvent("run.completed", pemu::trace::Severity::Info);
    }
  }

  // ========================================================
  // Time state
  // ========================================================

  [[nodiscard]]
  double time() const noexcept {
    return clock_.time();
  }

  [[nodiscard]]
  double timeStep() const noexcept {
    return transport_stepper_->timeStep();
  }

  [[nodiscard]]
  std::size_t step() const noexcept {
    return clock_.step();
  }

  [[nodiscard]]
  std::size_t totalSteps() const noexcept {
    return clock_.totalSteps();
  }

  [[nodiscard]]
  bool finished() const noexcept {
    return clock_.finished();
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
  // Current reaction workspace.
  //
  // After advanceOneStep():
  //
  // these correspond to the reaction state evaluated from
  // n^k, before n was advanced to n^(k+1).
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
  // Electrostatic state.
  //
  // After a completed timestep these correspond to n^k,
  // i.e. the fields actually used to advance to n^(k+1).
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
        pemu::trace::TraceAttribute{"dt", clock_.timeStep()},
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

  void emitStepCompleted(const linalg::SolverResult& result) noexcept {
    const std::array attributes{
        pemu::trace::TraceAttribute{"step",
                                  static_cast<std::uint64_t>(clock_.step())},
        pemu::trace::TraceAttribute{"time", clock_.time()},
        pemu::trace::TraceAttribute{"dt", clock_.timeStep()},
        pemu::trace::TraceAttribute{"relative_residual",
                                  result.relative_residual},
    };
    emitTrace("step.completed", pemu::trace::Severity::Info, attributes);
  }

  template <std::size_t N>
  void emitTrace(
      std::string_view name, pemu::trace::Severity severity,
      const std::array<pemu::trace::TraceAttribute, N>& attributes) noexcept {
    trace_sink_({.category = "simulation.fixed_step",
                 .name = name,
                 .severity = severity,
                 .attributes = attributes});
  }

  void validateConfiguration() const {
    // ----------------------------------------------------
    // FixedStepClock and FixedStepStepper MUST use exactly
    // the same timestep.
    // ----------------------------------------------------

    const double clock_dt = clock_.timeStep();

    const double solver_dt = transport_stepper_->timeStep();

    const double scale =
        std::max({1.0, std::abs(clock_dt), std::abs(solver_dt)});

    if (std::abs(clock_dt - solver_dt)

        > 1e-14 * scale) {

      throw std::invalid_argument(
          "FixedStepClock dt does not "
          "match FixedStep transport dt");
    }

    if (density_->size() == 0) {

      throw std::invalid_argument(
          "simulation density must "
          "contain species");
    }
  }

 private:
  // --------------------------------------------------------
  // External/non-owning simulation state.
  // --------------------------------------------------------

  physics::SpeciesCellFields* density_;

  const physics::ReactionNetwork* reaction_network_;

  Stepper* transport_stepper_;

  // --------------------------------------------------------
  // Rate law.
  //
  // Kept as a concrete template type instead of std::function
  // so later GPU/backend implementations are not forced
  // through type erasure.
  // --------------------------------------------------------

  ReactionRateEvaluator rate_evaluator_;

  // --------------------------------------------------------
  // Fixed time integration state.
  // --------------------------------------------------------

  FixedStepClock clock_;

  [[no_unique_address]] TraceSink trace_sink_;

  // --------------------------------------------------------
  // Workspaces.
  // --------------------------------------------------------

  physics::ReactionRateFields reaction_rates_;

  physics::SpeciesCellFields source_;
};

}  // namespace pemu::simulation
