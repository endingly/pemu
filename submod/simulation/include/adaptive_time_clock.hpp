#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace pemu::simulation {

class AdaptiveTimeClock {
 public:
  explicit AdaptiveTimeClock(double end_time)

      : end_time_(end_time) {
    if (end_time_ <= 0.0) {

      throw std::invalid_argument("end time must be positive");
    }
  }

  [[nodiscard]]
  double time() const noexcept {
    return static_cast<double>(time_);
  }

  [[nodiscard]]
  double endTime() const noexcept {
    return end_time_;
  }

  [[nodiscard]]
  double remainingTime() const noexcept {
    const double remaining = end_time_ - time();

    return std::max(0.0, remaining);
  }

  [[nodiscard]]
  std::size_t step() const noexcept {
    return step_;
  }

  [[nodiscard]]
  bool finished() const noexcept {
    const double tolerance = 1e-14 * std::max(1.0, end_time_);

    return remainingTime() <= tolerance;
  }

  void advance(double dt) {
    if (dt <= 0.0) {

      throw std::invalid_argument("time step must be positive");
    }

    if (finished()) {

      throw std::out_of_range(
          "adaptive clock "
          "already finished");
    }

    const double remaining = remainingTime();

    if (dt > remaining * (1.0 + 1e-12)) {

      throw std::out_of_range(
          "adaptive timestep "
          "exceeds remaining time");
    }

    time_ += static_cast<long double>(dt);

    ++step_;

    // ----------------------------------------------------
    // Snap to t_end to eliminate tiny floating remainder.
    // ----------------------------------------------------

    if (finished()) {

      time_ = static_cast<long double>(end_time_);
    }
  }

 private:
  long double time_{0.0L};

  double end_time_{};

  std::size_t step_{0};
};

}  // namespace pemu::simulation