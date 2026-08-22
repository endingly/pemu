#pragma once

#include <pemu/output/checkpoint/i_writer.hpp>

#include <cstddef>
#include <filesystem>
#include <stdexcept>

namespace pemu::simulation {

/** @brief Configures checkpoints committed by a simulation workflow. */
struct CheckpointOptions {
  const output::checkpoint::IWriter* writer{};
  std::filesystem::path path{"checkpoint/plasma.vtkhdf"};
  std::size_t every_steps{};
  bool write_final{true};
  bool overwrite{true};

  [[nodiscard]] bool enabled() const noexcept { return writer != nullptr; }
};

inline void validateCheckpointOptions(const CheckpointOptions& options) {
  if (!options.enabled()) {
    return;
  }
  if (options.path.empty()) {
    throw std::invalid_argument("checkpoint path must not be empty");
  }
  if (options.every_steps == 0 && !options.write_final) {
    throw std::invalid_argument(
        "checkpoint must select periodic or final states");
  }
  if (options.every_steps != 0 && !options.overwrite) {
    throw std::invalid_argument(
        "periodic checkpoint output requires atomic overwrite");
  }
}

[[nodiscard]] inline bool shouldWriteCheckpoint(
    const CheckpointOptions& options, std::size_t step,
    bool terminal) noexcept {
  return options.enabled() &&
         ((options.every_steps != 0 && step % options.every_steps == 0) ||
          (terminal && options.write_final));
}

}  // namespace pemu::simulation
