#pragma once

#include <pemu/equation/adaptive_step_multi_species_drift_diffusion_stepper.hpp>
#include <pemu/equation/time_integration/adaptive_time_step_controller.hpp>
#include <pemu/physics/reaction.hpp>
#include <pemu/physics/species.hpp>
#include <pemu/simulation/adaptive_time_clock.hpp>

#include <cstddef>
#include <stdexcept>
#include <utility>

namespace pemu::simulation {

template <typename ReactionRateEvaluator>
class AdaptiveStepPlasmaSimulation {
 public:
  using Stepper = equation::AdaptiveStepMultiSpeciesDriftDiffusionStepper;

  using TimeStepProposal = equation::time_integration::TimeStepProposal;

  AdaptiveStepPlasmaSimulation(physics::SpeciesCellFields& density,

                               const physics::ReactionNetwork& reaction_network,

                               Stepper& transport_stepper,

                               ReactionRateEvaluator rate_evaluator,

                               AdaptiveTimeClock clock)

      : density_(&density),

        reaction_network_(&reaction_network),

        transport_stepper_(&transport_stepper),

        rate_evaluator_(std::move(rate_evaluator)),

        clock_(std::move(clock)),

        reaction_rates_(density.mesh(), reaction_network.size(), 0.0),

        source_(density.mesh(), density.size(), 0.0) {
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

    // ----------------------------------------------------
    // 1.
    //
    // n^k -> rho^k -> phi^k -> E^k
    // ----------------------------------------------------

    const auto result = transport_stepper_->prepareElectrostatics(*density_);

    if (!result.success()) {
      return result;
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

    // ----------------------------------------------------
    // 3.
    //
    // R -> S
    // ----------------------------------------------------

    source_.fill(0.0);

    reaction_network_->accumulateSources(reaction_rates_.span(), source_);

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

    // ----------------------------------------------------
    // 5.
    //
    // Commit time AFTER physical state update succeeds.
    // ----------------------------------------------------

    clock_.advance(proposal.dt);

    return result;
  }

  // ========================================================
  // Run until t_end.
  // ========================================================

  void run() {
    while (!clock_.finished()) {

      const auto result = advanceOneStep();

      if (!result.success()) {

        throw std::runtime_error(
            "adaptive plasma "
            "simulation failed");
      }
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

  // --------------------------------------------------------
  // Workspaces
  // --------------------------------------------------------

  physics::ReactionRateFields reaction_rates_;

  physics::SpeciesCellFields source_;
};

}  // namespace pemu::simulation