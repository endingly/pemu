#pragma once

#include <pemu/output/i_field_output_writer.hpp>
#include <pemu/output/output_trace.hpp>

namespace pemu::output {

class VtkHdfWriter final : public IFieldOutputWriter {
 public:
  [[nodiscard]]
  OutputRecord write(const FieldDumpRequest& request) const override;

  template <trace::TraceSink Sink>
  [[nodiscard]] OutputRecord writeAndTrace(const FieldDumpRequest& request,
                                           Sink& sink) const {
    auto record = write(request);
    traceOutputCompleted(sink, record);
    return record;
  }
};

}  // namespace pemu::output
