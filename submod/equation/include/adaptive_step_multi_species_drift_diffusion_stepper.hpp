#pragma once

#include <cmath>
#include <memory>
#include <optional>
#include <pemu/equation/detail/explicit_multi_species_drift_diffusion_operator.hpp>
#include <pemu/equation/time_integration/adaptive_time_step_controller.hpp>
#include <pemu/field/plasma_field_metadata.hpp>
#include <pemu/linalg/i_solver.hpp>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pemu::equation {

class AdaptiveStepMultiSpeciesDriftDiffusionStepper {
 public:
  /**
   * @brief Creates an adaptive multi-species drift-diffusion stepper.
   * @param pure_neumann_options Gauge and compatibility tolerances required
   * when every potential boundary is Neumann.
   */
  AdaptiveStepMultiSpeciesDriftDiffusionStepper(
      const mesh::IMesh& mesh, const physics::SpeciesSet& species,

      double permittivity,

      boundary::BoundaryConditionSet potential_bc,

      std::vector<boundary::BoundaryConditionSet> species_bc,

      std::unique_ptr<linalg::ISolver> poisson_backend,

      time_integration::AdaptiveTimeStepConfig time_step_config,

      field::PlasmaFieldMetadata field_metadata = {},

      std::optional<PureNeumannOptions> pure_neumann_options = std::nullopt)

      : transport_operator_(
            mesh, species, permittivity, std::move(potential_bc),
            std::move(species_bc), std::move(poisson_backend),
            std::move(field_metadata), std::move(pure_neumann_options)),

        controller_(time_step_config) {}

  AdaptiveStepMultiSpeciesDriftDiffusionStepper(
      const AdaptiveStepMultiSpeciesDriftDiffusionStepper&) = delete;

  AdaptiveStepMultiSpeciesDriftDiffusionStepper& operator=(
      const AdaptiveStepMultiSpeciesDriftDiffusionStepper&) = delete;

  AdaptiveStepMultiSpeciesDriftDiffusionStepper(
      AdaptiveStepMultiSpeciesDriftDiffusionStepper&&) = delete;

  AdaptiveStepMultiSpeciesDriftDiffusionStepper& operator=(
      AdaptiveStepMultiSpeciesDriftDiffusionStepper&&) = delete;

  // ========================================================
  // Phase 1:
  //
  // n^k -> rho^k -> phi^k -> E^k
  //
  // Reaction evaluator should run AFTER this call and BEFORE
  // advancePrepared().
  // ========================================================

  [[nodiscard]]
  linalg::SolverResult prepareElectrostatics(
      const physics::SpeciesCellFields& density) {
    return transport_operator_.updateElectrostatics(density);
  }

  // ========================================================
  // Determine dt WITHOUT advancing.
  //
  // Requires prepareElectrostatics() first.
  // ========================================================

  [[nodiscard]]
  time_integration::TimeStepProposal proposeTimeStep(
      const physics::SpeciesCellFields& density,
      const physics::SpeciesCellFields& source, double remaining_time) {
    const double transport_limit =
        transport_operator_.maxStableTransportTimeStep();

    const double positivity_limit =
        transport_operator_.positivityTimeStepLimit(density, source);

    return controller_.propose(transport_limit, positivity_limit, previous_dt_,
                               remaining_time);
  }

  // ========================================================
  // Phase 2:
  //
  // Select adaptive dt and advance using prepared E^k.
  //
  // Returns the dt actually used.
  // ========================================================

  [[nodiscard]]
  time_integration::TimeStepProposal advancePrepared(
      physics::SpeciesCellFields& density,
      const physics::SpeciesCellFields& source, double remaining_time) {
    const auto proposal = proposeTimeStep(density, source, remaining_time);

    transport_operator_.advanceTransport(density, source, proposal.dt);

    previous_dt_ = proposal.dt;

    last_proposal_ = proposal;

    has_last_proposal_ = true;

    return proposal;
  }

  // ========================================================
  // Convenience API.
  //
  // Use only when source is already known independently of
  // freshly computed E^k.
  //
  // For E-dependent chemistry use:
  //
  //     prepareElectrostatics()
  //     evaluate reactions
  //     advancePrepared()
  // ========================================================

  [[nodiscard]]
  linalg::SolverResult step(physics::SpeciesCellFields& density,
                            const physics::SpeciesCellFields& source,
                            double remaining_time) {
    const auto result = prepareElectrostatics(density);

    if (!result.success()) {
      return result;
    }

    (void)advancePrepared(density, source, remaining_time);

    return result;
  }

  // ========================================================
  // Adaptive state
  // ========================================================

  [[nodiscard]]
  double previousTimeStep() const noexcept {
    return previous_dt_;
  }

  [[nodiscard]]
  bool hasLastTimeStepProposal() const noexcept {
    return has_last_proposal_;
  }

  [[nodiscard]]
  const time_integration::TimeStepProposal& lastTimeStepProposal() const {
    if (!has_last_proposal_) {

      throw std::logic_error(
          "adaptive stepper has "
          "not advanced yet");
    }

    return last_proposal_;
  }

  void resetTimeStepHistory() noexcept {
    previous_dt_ = 0.0;

    has_last_proposal_ = false;
  }

  /** @brief Restores the growth-limiter history required by a restart. */
  void restoreTimeStepHistory(double previous_dt) {
    if (!std::isfinite(previous_dt) || previous_dt < 0.0) {
      throw std::invalid_argument(
          "previous timestep must be finite and non-negative");
    }
    previous_dt_ = previous_dt;
    has_last_proposal_ = false;
  }

  // ========================================================
  // Shared physics accessors
  // ========================================================

  [[nodiscard]]
  const field::CellField<double>& chargeDensity() const noexcept {
    return transport_operator_.chargeDensity();
  }

  [[nodiscard]]
  const field::CellField<double>& potential() const noexcept {
    return transport_operator_.potential();
  }

  [[nodiscard]]
  const field::FaceField<double>& electricFieldNormal() const noexcept {
    return transport_operator_.electricFieldNormal();
  }

  [[nodiscard]]
  const field::FaceField<double>& driftVelocityNormal(
      physics::SpeciesId id) const {
    return transport_operator_.driftVelocityNormal(id);
  }

  [[nodiscard]]
  bool electrostaticsReady() const noexcept {
    return transport_operator_.electrostaticsReady();
  }

  [[nodiscard]]
  const field::PlasmaFieldMetadata& fieldMetadata() const noexcept {
    return transport_operator_.fieldMetadata();
  }

  [[nodiscard]]
  const physics::SpeciesSet& species() const noexcept {
    return transport_operator_.species();
  }

  [[nodiscard]]
  const mesh::IMesh& mesh() const noexcept {
    return transport_operator_.mesh();
  }

 private:
  detail::ExplicitMultiSpeciesDriftDiffusionOperator transport_operator_;

  time_integration::AdaptiveTimeStepController controller_;

  double previous_dt_{0.0};

  time_integration::TimeStepProposal last_proposal_{};

  bool has_last_proposal_{false};
};

}  // namespace pemu::equation
