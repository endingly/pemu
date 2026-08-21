#pragma once

#include <pemu/trace/trace.hpp>

#include <iosfwd>

namespace pemu::trace {

class OstreamTraceSink {
 public:
  explicit OstreamTraceSink(std::ostream& stream) noexcept;

  void operator()(const TraceEvent& event) noexcept;

  [[nodiscard]] bool failed() const noexcept;

 private:
  std::ostream* stream_;
  bool failed_{false};
};

static_assert(TraceSink<OstreamTraceSink>);

}  // namespace pemu::trace
