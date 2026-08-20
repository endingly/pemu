#pragma once

#include <pemu/equation/multi_species_drift_diffusion_stepper.hpp>
#include <pemu/physics/reaction.hpp>
#include <pemu/physics/species.hpp>
#include <pemu/simulation/fixed_step_clock.hpp>

#include <cmath>
#include <stdexcept>
#include <utility>

namespace pemu::simulation {

template <typename ReactionRateEvaluator>
class PlasmaSimulation {
 public:
  PlasmaSimulation(physics::SpeciesCellFields& density,
                   const physics::ReactionNetwork& reaction_network,
                   equation::MultiSpeciesDriftDiffusionStepper& transport,
                   ReactionRateEvaluator rate_evaluator, FixedStepClock clock)
      : density_(&density),
        reaction_network_(&reaction_network),
        transport_(&transport),
        rate_evaluator_(std::move(rate_evaluator)),
        clock_(std::move(clock)),
        reaction_rates_(density.mesh(), reaction_network.size(), 0.0),
        source_(density.mesh(), density.size(), 0.0) {
    const double scale = std::max(1.0, std::abs(clock_.timeStep()));

    if (std::abs(clock_.timeStep() - transport_->timeStep()) > 1e-14 * scale) {

      throw std::invalid_argument(
          "simulation clock dt does not "
          "match transport solver dt");
    }
  }

  // ========================================================
  // One complete plasma timestep:
  //
  // n^k
  //   ↓
  // electrostatics
  //   ↓
  // phi^k, E^k
  //   ↓
  // reaction rates R^k
  //   ↓
  // source S^k = nu R
  //   ↓
  // transport
  //   ↓
  // n^(k+1)
  // ========================================================

  linalg::SolverResult advanceOneStep() {
    if (clock_.finished()) {
      throw std::out_of_range("simulation already finished");
    }

    // ----------------------------------------------------
    // 1. n^k -> rho^k -> phi^k -> E^k
    // ----------------------------------------------------

    const auto result = transport_->updateElectrostatics(*density_);

    if (!result.success()) {
      return result;
    }

    // ----------------------------------------------------
    // 2. Evaluate reaction rates using the SAME state:
    //
    //     n^k, phi^k, E^k
    // ----------------------------------------------------

    reaction_rates_.fill(0.0);

    rate_evaluator_(*density_,

                    transport_->potential(),

                    transport_->electricFieldNormal(),

                    reaction_rates_);

    // ----------------------------------------------------
    // 3. R -> S
    // ----------------------------------------------------

    source_.fill(0.0);

    reaction_network_->accumulateSources(reaction_rates_.span(), source_);

    // ----------------------------------------------------
    // 4. n^k -> n^(k+1)
    //
    // using the already prepared E^k.
    // ----------------------------------------------------

    transport_->advanceTransport(*density_, source_);

    // ----------------------------------------------------
    // 5. Commit simulation time.
    // ----------------------------------------------------

    clock_.advance();

    return result;
  }

  void run() {
    while (!clock_.finished()) {

      const auto result = advanceOneStep();

      if (!result.success()) {
        throw std::runtime_error(
            "plasma simulation "
            "timestep failed");
      }
    }
  }

  [[nodiscard]]
  double time() const noexcept {
    return clock_.time();
  }

  [[nodiscard]]
  std::size_t step() const noexcept {
    return clock_.step();
  }

  [[nodiscard]]
  bool finished() const noexcept {
    return clock_.finished();
  }

  [[nodiscard]]
  const physics::SpeciesCellFields& density() const noexcept {
    return *density_;
  }

  [[nodiscard]]
  const physics::SpeciesCellFields& source() const noexcept {
    return source_;
  }

  [[nodiscard]]
  const physics::ReactionRateFields& reactionRates() const noexcept {
    return reaction_rates_;
  }

 private:
  physics::SpeciesCellFields* density_;
  const physics::ReactionNetwork* reaction_network_;
  equation::MultiSpeciesDriftDiffusionStepper* transport_;
  ReactionRateEvaluator rate_evaluator_;
  FixedStepClock clock_;
  physics::ReactionRateFields reaction_rates_;
  physics::SpeciesCellFields source_;
};

}  // namespace pemu::simulation