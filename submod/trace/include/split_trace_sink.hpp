#pragma once

#include <pemu/trace/trace.hpp>

#include <cstddef>
#include <type_traits>
#include <utility>

namespace pemu::trace {

// Routes statistics to their own sink while all ordinary trace and diagnostic
// events share the diagnostic sink.  EventKind remains free to describe an
// abnormal statistics sample; routing never depends on names or severity.
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
      statistics_dirty_ = true;
      return;
    }

    diagnostic_sink_(event);
    flushStatisticsAfter(event);
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
  template <typename Sink>
  static void flushIfSupported(Sink& sink) noexcept {
    if constexpr (requires {
                    { sink.flush() } noexcept -> std::same_as<void>;
                  }) {
      sink.flush();
    }
  }

  void flushStatisticsAfter(const TraceEvent& event) noexcept {
    if (event.domain != DiagDomain::simulation) {
      return;
    }

    if (event.name == "step.completed") {
      ++completed_steps_since_flush_;
      if (completed_steps_since_flush_ == 2) {
        completed_steps_since_flush_ = 0;
        if (statistics_dirty_) {
          flushIfSupported(statistics_sink_);
          statistics_dirty_ = false;
        }
      }
      return;
    }

    if (event.name == "run.completed" || event.name == "run.failed") {
      completed_steps_since_flush_ = 0;
      if (statistics_dirty_) {
        flushIfSupported(statistics_sink_);
        statistics_dirty_ = false;
      }
    }
  }

  [[no_unique_address]] DiagnosticSink diagnostic_sink_;
  [[no_unique_address]] StatisticsSink statistics_sink_;
  std::size_t completed_steps_since_flush_{0};
  bool statistics_dirty_{false};
};

template <TraceSink DiagnosticSink, TraceSink StatisticsSink>
SplitTraceSink(DiagnosticSink, StatisticsSink)
    -> SplitTraceSink<DiagnosticSink, StatisticsSink>;

static_assert(TraceSink<SplitTraceSink<NullTraceSink, NullTraceSink>>);

}  // namespace pemu::trace
