#pragma once

#include <pemu/output/output_types.hpp>

namespace pemu::output {

class IFieldOutputWriter {
 public:
  virtual ~IFieldOutputWriter() = default;

  [[nodiscard]]
  virtual OutputRecord write(const FieldDumpRequest& request) const = 0;
};

}  // namespace pemu::output
