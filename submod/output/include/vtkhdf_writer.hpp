#pragma once

#include <pemu/output/i_field_output_writer.hpp>
#include <pemu/output/output_trace.hpp>

namespace pemu::output {

/** @brief Creates single-file VTKHDF field time series. */
class VtkHdfWriter final : public IFieldOutputWriter {
 public:
  /** @copydoc IFieldOutputWriter::openSeries */
  [[nodiscard]]
  std::unique_ptr<IFieldOutputSeries> openSeries(
      const FieldSeriesRequest& request) const override;

  /**
   * @brief Writes one snapshot and emits its completed-file trace event.
   * @tparam Sink Trace-sink type.
   * @param request Complete standalone field-dump request.
   * @param sink Trace sink receiving the completion event.
   * @return Record identifying the completed file.
   */
  template <trace::TraceSink Sink>
  [[nodiscard]] OutputRecord writeAndTrace(const FieldDumpRequest& request,
                                           Sink& sink) const {
    auto record = write(request);
    traceOutputCompleted(sink, record);
    return record;
  }
};

}  // namespace pemu::output
