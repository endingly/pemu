#pragma once

#include <pemu/output/output_types.hpp>

#include <filesystem>

namespace pemu::output {

// Checkpoint output intentionally has a separate interface from field dumps.
// The concrete simulation State type and persistence format are deferred until
// checkpoint/restart is implemented.
template <typename State>
class ICheckpointWriter {
 public:
  virtual ~ICheckpointWriter() = default;

  [[nodiscard]]
  virtual OutputRecord write(const State& state,
                             const std::filesystem::path& path,
                             OutputStamp stamp) const = 0;
};

}  // namespace pemu::output
