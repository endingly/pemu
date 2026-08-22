#pragma once

#include <cstdint>
#include <filesystem>

namespace pemu::output {

/** @brief Identifies the simulation state persisted by an output operation. */
struct OutputStamp {
  std::uint64_t step{};
  double time{};
};

/** @brief Describes one completed output artifact. */
struct OutputRecord {
  std::filesystem::path path;
  OutputStamp stamp;
};

}  // namespace pemu::output
