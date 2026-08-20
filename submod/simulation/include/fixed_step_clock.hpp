#pragma once

#include <cstddef>
#include <stdexcept>

namespace pemu::simulation {

class FixedStepClock {
 public:
  FixedStepClock(double dt, std::size_t total_steps)
      : dt_(dt), total_steps_(total_steps) {
    if (dt <= 0.0) {
      throw std::invalid_argument("time step must be positive");
    }
  }

  [[nodiscard]]
  double timeStep() const noexcept {
    return dt_;
  }

  [[nodiscard]]
  std::size_t step() const noexcept {
    return step_;
  }

  [[nodiscard]]
  std::size_t totalSteps() const noexcept {
    return total_steps_;
  }

  [[nodiscard]]
  double time() const noexcept {
    return static_cast<double>(step_) * dt_;
  }

  [[nodiscard]]
  double endTime() const noexcept {
    return static_cast<double>(total_steps_) * dt_;
  }

  [[nodiscard]]
  bool finished() const noexcept {
    return step_ >= total_steps_;
  }

  void advance() {
    if (finished()) {
      throw std::out_of_range("simulation clock already finished");
    }

    ++step_;
  }

 private:
  double dt_;

  std::size_t step_{0};
  std::size_t total_steps_;
};

}  // namespace pemu::simulation