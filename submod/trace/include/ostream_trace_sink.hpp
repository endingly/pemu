#pragma once

#include <pemu/trace/renderer/diagnostic_event_renderer.hpp>
#include <pemu/trace/renderer/ordinary_event_renderer.hpp>
#include <pemu/trace/renderer/statistics_event_renderer.hpp>

#include <cstddef>
#include <iosfwd>

namespace pemu::trace {

class OstreamTraceSink {
 public:
  explicit OstreamTraceSink(std::ostream& stream) noexcept;

  void operator()(const TraceEvent& event) noexcept;

  void flush() noexcept;

  [[nodiscard]] bool failed() const noexcept;

 private:
  [[nodiscard]] bool shouldFlushAfter(const TraceEvent& event) noexcept;

  std::ostream* stream_;
  renderer::OrdinaryEventRenderer ordinary_renderer_;
  renderer::DiagnosticEventRenderer diagnostic_renderer_;
  renderer::StatisticsEventRenderer statistics_renderer_;
  bool event_header_written_{false};
  std::size_t completed_steps_since_flush_{0};
  bool failed_{false};
};

static_assert(TraceSink<OstreamTraceSink>);

}  // namespace pemu::trace
