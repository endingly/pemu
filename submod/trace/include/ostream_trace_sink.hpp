#pragma once

#include <pemu/trace/renderer/diagnostic_event_renderer.hpp>
#include <pemu/trace/renderer/ordinary_event_renderer.hpp>
#include <pemu/trace/renderer/statistics_event_renderer.hpp>

#include <iosfwd>

namespace pemu::trace {

// Concrete synchronous rendering backend. This class owns formatting and
// stream-error state; it is not the runtime-polymorphism abstraction (see
// AnyTraceSink) and not a channel router (see SplitTraceSink).
class OstreamTraceSink {
 public:
  explicit OstreamTraceSink(std::ostream& stream) noexcept;

  void operator()(const TraceEvent& event) noexcept;

  void flush() noexcept;

  [[nodiscard]] bool failed() const noexcept;

 private:
  std::ostream* stream_;
  renderer::OrdinaryEventRenderer ordinary_renderer_;
  renderer::DiagnosticEventRenderer diagnostic_renderer_;
  renderer::StatisticsEventRenderer statistics_renderer_;
  bool event_header_written_{false};
  bool failed_{false};
};

static_assert(TraceSink<OstreamTraceSink>);

}  // namespace pemu::trace
