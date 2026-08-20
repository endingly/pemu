#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace pemu::equation::time_integration {

struct AdaptiveTimeStepConfig {
  //
  // Applied to stability/positivity limits.
  //
  // Typical:
  //
  //     0.8 ~ 0.95
  //
  double safety{0.9};

  //
  // Hard lower bound except the final remainder step.
  //
  double min_dt{1e-12};

  //
  // User/global upper bound.
  //
  double max_dt{std::numeric_limits<double>::infinity()};

  //
  // Prevent:
  //
  //     dt_new >> dt_old
  //
  // after a temporarily restrictive timestep.
  //
  double max_growth{1.5};
};

struct TimeStepProposal {
  double dt{};

  double transport_limit{std::numeric_limits<double>::infinity()};

  double positivity_limit{std::numeric_limits<double>::infinity()};

  double stability_limit{std::numeric_limits<double>::infinity()};
};

class AdaptiveTimeStepController {
 public:
  explicit AdaptiveTimeStepController(AdaptiveTimeStepConfig config)
      : config_(config) {
    validateConfig();
  }

  [[nodiscard]]
  const AdaptiveTimeStepConfig& config() const noexcept {
    return config_;
  }

  [[nodiscard]]
  TimeStepProposal propose(double transport_limit, double positivity_limit,
                           double previous_dt, double remaining_time) const {
    validateLimit(transport_limit, "transport limit");

    validateLimit(positivity_limit, "positivity limit");

    if (remaining_time <= 0.0) {
      throw std::invalid_argument(
          "remaining time "
          "must be positive");
    }

    if (previous_dt < 0.0) {
      throw std::invalid_argument(
          "previous dt must "
          "not be negative");
    }

    if (positivity_limit == 0.0 || transport_limit == 0.0) {

      throw std::runtime_error(
          "no positive timestep "
          "is admissible");
    }

    // ----------------------------------------------------
    // Both restrictions matter.
    //
    // positivity_limit includes reaction depletion,
    // while transport_limit remains necessary even in
    // zero-density cells.
    // ----------------------------------------------------

    const double stability_limit = std::min(transport_limit, positivity_limit);

    double dt = config_.max_dt;

    if (std::isfinite(stability_limit)) {

      dt = std::min(dt,

                    config_.safety * stability_limit);
    }

    // ----------------------------------------------------
    // Growth limiter.
    // ----------------------------------------------------

    if (previous_dt > 0.0) {

      dt = std::min(dt,

                    previous_dt * config_.max_growth);
    }

    // ----------------------------------------------------
    // Final timestep lands exactly at t_end.
    // ----------------------------------------------------

    dt = std::min(dt, remaining_time);

    if (!std::isfinite(dt) || dt <= 0.0) {

      throw std::runtime_error(
          "failed to determine "
          "adaptive timestep");
    }

    // ----------------------------------------------------
    // A final remainder smaller than min_dt is allowed.
    //
    // Otherwise dt < min_dt means the simulation has
    // become too stiff for the configured explicit method.
    // ----------------------------------------------------

    if (dt < config_.min_dt && remaining_time > config_.min_dt) {

      throw std::runtime_error(
          "adaptive timestep "
          "fell below min_dt");
    }

    return {.dt = dt,

            .transport_limit = transport_limit,

            .positivity_limit = positivity_limit,

            .stability_limit = stability_limit};
  }

 private:
  void validateConfig() const {
    if (config_.safety <= 0.0 || config_.safety > 1.0) {

      throw std::invalid_argument(
          "safety must lie "
          "in (0, 1]");
    }

    if (config_.min_dt <= 0.0) {

      throw std::invalid_argument("min_dt must be positive");
    }

    if (config_.max_dt < config_.min_dt) {

      throw std::invalid_argument(
          "max_dt must be "
          ">= min_dt");
    }

    if (config_.max_growth < 1.0) {

      throw std::invalid_argument(
          "max_growth must "
          "be >= 1");
    }
  }

  static void validateLimit(double value, const char* name) {
    if (std::isnan(value) || value < 0.0) {

      throw std::invalid_argument(name);
    }
  }

 private:
  AdaptiveTimeStepConfig config_;
};

}  // namespace pemu::equation::time_integration