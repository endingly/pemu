#pragma once

namespace pemu::simulation {

/** @brief Lifecycle state shared by synchronous simulation workflows. */
enum class SimulationState {
  ready,
  running,
  paused,
  stopped,
  completed,
  failed,
};

/** @brief Returns whether no further workflow transition can resume a run. */
[[nodiscard]] constexpr bool isTerminal(SimulationState state) noexcept {
  return state == SimulationState::stopped ||
         state == SimulationState::completed ||
         state == SimulationState::failed;
}

}  // namespace pemu::simulation
