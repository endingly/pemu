#pragma once

#include <pemu/trace/trace.hpp>

#include <type_traits>
#include <utility>

namespace pemu::trace {

// Compile-time composition policy: routes statistics to their own sink while
// ordinary trace and diagnostic events share the diagnostic sink. Both child
// types stay visible to the compiler, so composition adds no virtual dispatch.
// EventKind remains free to describe an abnormal statistics sample; routing
// never depends on names or severity.
template <TraceSink DiagnosticSink, TraceSink StatisticsSink>
class SplitTraceSink {
  static constexpr bool nothrow_move_constructible =
      std::is_nothrow_move_constructible_v<DiagnosticSink> &&
      std::is_nothrow_move_constructible_v<StatisticsSink>;

 public:
  SplitTraceSink(
      DiagnosticSink diagnostic_sink,
      StatisticsSink statistics_sink) noexcept(nothrow_move_constructible)
      : diagnostic_sink_(std::move(diagnostic_sink)),
        statistics_sink_(std::move(statistics_sink)) {}

  void operator()(const TraceEvent& event) noexcept {
    if (event.output_channel == OutputChannel::Statistics) {
      statistics_sink_(event);
      return;
    }

    diagnostic_sink_(event);
  }

  void flush() noexcept {
    diagnostic_sink_.flush();
    statistics_sink_.flush();
  }

  [[nodiscard]] DiagnosticSink& diagnosticSink() noexcept {
    return diagnostic_sink_;
  }

  [[nodiscard]] const DiagnosticSink& diagnosticSink() const noexcept {
    return diagnostic_sink_;
  }

  [[nodiscard]] StatisticsSink& statisticsSink() noexcept {
    return statistics_sink_;
  }

  [[nodiscard]] const StatisticsSink& statisticsSink() const noexcept {
    return statistics_sink_;
  }

 private:
  [[no_unique_address]] DiagnosticSink diagnostic_sink_;
  [[no_unique_address]] StatisticsSink statistics_sink_;
};

template <TraceSink DiagnosticSink, TraceSink StatisticsSink>
SplitTraceSink(DiagnosticSink, StatisticsSink)
    -> SplitTraceSink<DiagnosticSink, StatisticsSink>;

static_assert(TraceSink<SplitTraceSink<NullTraceSink, NullTraceSink>>);

}  // namespace pemu::trace
