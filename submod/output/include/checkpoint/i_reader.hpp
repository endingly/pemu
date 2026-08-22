#pragma once

#include <pemu/output/checkpoint/types.hpp>

namespace pemu::output::checkpoint {

class IReader {
 public:
  virtual ~IReader() = default;

  [[nodiscard]] virtual Manifest inspect(
      const std::filesystem::path& path) const = 0;

  [[nodiscard]] virtual OutputRecord restore(
      const std::filesystem::path& path,
      const RestoreRequest& request) const = 0;
};

}  // namespace pemu::output::checkpoint
