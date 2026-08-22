#pragma once

#include <pemu/output/checkpoint/types.hpp>

namespace pemu::output::checkpoint {

class IWriter {
 public:
  virtual ~IWriter() = default;
  [[nodiscard]] virtual OutputRecord write(
      const WriteRequest& request) const = 0;
};

}  // namespace pemu::output::checkpoint
