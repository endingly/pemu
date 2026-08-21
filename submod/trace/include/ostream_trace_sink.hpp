#pragma once

#include <pemu/trace/trace.hpp>

#include <cstddef>
#include <iosfwd>

namespace pemu::trace {

class OstreamTraceSink {
 public:
  explicit OstreamTraceSink(std::ostream& stream) noexcept;

  void operator()(const TraceEvent& event) noexcept;

 [[nodiscard]] bool failed() const noexcept;

 private:
  [[nodiscard]] bool shouldFlushAfter(const TraceEvent& event) noexcept;

  std::ostream* stream_;
  bool header_written_{false};
  std::size_t completed_steps_since_flush_{0};
  bool failed_{false};
};

static_assert(TraceSink<OstreamTraceSink>);

}  // namespace pemu::trace
