#pragma once

#include <pemu/output/i_field_output_writer.hpp>

#include <cstddef>
#include <filesystem>
#include <stdexcept>
#include <string>

namespace pemu::simulation {

/**
 * @brief Configures scheduled field snapshots emitted by a plasma simulation.
 *
 * The writer is non-owning and must outlive the simulation. A null writer
 * disables all field output without changing the numerical stepping path.
 */
struct FieldOutputOptions {
  const output::IFieldOutputWriter* writer{};
  std::filesystem::path directory{"output"};
  std::string file_stem{"plasma"};
  std::size_t every_steps{};
  bool write_initial{true};
  bool write_final{true};
  bool overwrite{};
  bool include_face_cell_visualization{true};

  /** @brief Returns whether a field-output writer has been configured. */
  [[nodiscard]]
  bool enabled() const noexcept {
    return writer != nullptr;
  }
};

/**
 * @brief Validates a field-output configuration before a simulation starts.
 * @param options Field-output configuration to validate.
 * @throws std::invalid_argument if an enabled configuration has no file stem
 * or selects no snapshot times.
 */
inline void validateFieldOutputOptions(const FieldOutputOptions& options) {
  if (!options.enabled()) {
    return;
  }
  if (options.file_stem.empty()) {
    throw std::invalid_argument("field output file stem must not be empty");
  }
  if (options.every_steps == 0 && !options.write_initial &&
      !options.write_final) {
    throw std::invalid_argument(
        "field output must select an initial, periodic, or final snapshot");
  }
}

}  // namespace pemu::simulation
