#pragma once

#include <pemu/trace/trace.hpp>

#include <ostream>
#include <string_view>
#include <type_traits>
#include <variant>

namespace pemu::trace {

class OstreamTraceSink {
 public:
  explicit OstreamTraceSink(std::ostream& stream) noexcept : stream_(&stream) {}

  void operator()(const TraceEvent& event) noexcept {
    try {
      *stream_ << '[' << severityName(event.severity) << "] " << event.category
               << '.' << event.name;

      for (const auto& attribute : event.attributes) {
        *stream_ << ' ' << attribute.name << '=';
        std::visit([this](const auto& value) { *stream_ << value; },
                   attribute.value);
      }

      *stream_ << '\n';
      failed_ = failed_ || !stream_->good();
    } catch (...) {
      failed_ = true;
    }
  }

  [[nodiscard]] bool failed() const noexcept { return failed_; }

 private:
  [[nodiscard]] static constexpr std::string_view severityName(
      Severity severity) noexcept {
    switch (severity) {
      case Severity::Trace:
        return "trace";
      case Severity::Debug:
        return "debug";
      case Severity::Info:
        return "info";
      case Severity::Warning:
        return "warning";
      case Severity::Error:
        return "error";
      case Severity::Critical:
        return "critical";
    }

    return "unknown";
  }

  std::ostream* stream_;
  bool failed_{false};
};

static_assert(TraceSink<OstreamTraceSink>);

}  // namespace pemu::trace
