#pragma once

#include <pemu/output/dump/i_writer.hpp>
#include <pemu/output/dump/trace.hpp>

namespace pemu::output::dump {

/** @brief Creates single-file VTKHDF visualization time series. */
class VtkHdfWriter final : public IWriter {
 public:
  [[nodiscard]] std::unique_ptr<ISeries> openSeries(
      const SeriesRequest& request) const override;

  template <trace::TraceSink Sink>
  [[nodiscard]] OutputRecord writeAndTrace(const Request& request,
                                           Sink& sink) const {
    auto record = write(request);
    traceCompleted(sink, record);
    return record;
  }
};

}  // namespace pemu::output::dump
