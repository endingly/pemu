#pragma once

#include <pemu/trace/diag.hpp>

#include <llnl-units/units.hpp>

#include <concepts>
#include <cstdint>
#include <optional>
#include <span>
#include <string_view>
#include <variant>

namespace pemu::trace {

enum class Severity : std::uint8_t {
  Trace,
  Debug,
  Info,
  Warning,
  Error,
  Critical,
};

}  // namespace pemu::trace

namespace pemu {

[[nodiscard]] constexpr std::string_view to_string(
    trace::Severity severity) noexcept {
  switch (severity) {
    case trace::Severity::Trace:
      return "Trace";
    case trace::Severity::Debug:
      return "Debug";
    case trace::Severity::Info:
      return "Info";
    case trace::Severity::Warning:
      return "Warning";
    case trace::Severity::Error:
      return "Error";
    case trace::Severity::Critical:
      return "Critical";
  }
  return "Unknown";
}

}  // namespace pemu

namespace pemu::trace {

enum class EventKind : std::uint8_t {
  Trace,
  Diagnostic,
};

enum class OutputChannel : std::uint8_t {
  Diagnostic,
  Statistics,
};

}  // namespace pemu::trace

namespace pemu {

[[nodiscard]] constexpr std::string_view to_string(
    trace::EventKind kind) noexcept {
  switch (kind) {
    case trace::EventKind::Trace:
      return "Trace";
    case trace::EventKind::Diagnostic:
      return "Diagnostic";
  }
  return "Unknown";
}

[[nodiscard]] constexpr std::string_view to_string(
    trace::OutputChannel channel) noexcept {
  switch (channel) {
    case trace::OutputChannel::Diagnostic:
      return "Diagnostic";
    case trace::OutputChannel::Statistics:
      return "Statistics";
  }
  return "Unknown";
}

}  // namespace pemu

namespace pemu::trace {

using TraceValue =
    std::variant<bool, std::int64_t, std::uint64_t, double, std::string_view>;

struct TraceAttribute {
  std::string_view name;
  TraceValue value;
  std::optional<units::precise_unit> unit{};
};

// TraceEvent and its attributes are non-owning views. A sink must consume them
// synchronously. A sink that retains events must copy names and string values.
struct TraceEvent {
  EventKind kind{EventKind::Trace};
  OutputChannel output_channel{OutputChannel::Diagnostic};
  DiagDomain domain{DiagDomain::trace};
  std::string_view category;
  std::string_view name;
  std::string_view message{};
  Severity severity{Severity::Trace};
  std::span<const TraceAttribute> attributes{};
};

// Creates a diagnostic event without attributes. This form is suitable for
// attaching to a returned result when its string views refer to static storage.
// Runtime numeric context remains in the result and can be added by the
// consuming layer when the event is emitted.
[[nodiscard]] constexpr TraceEvent makeDiagnosticEvent(
    DiagDomain domain, std::string_view category, std::string_view name,
    std::string_view message, Severity severity = Severity::Error) noexcept {
  return {
      .kind = EventKind::Diagnostic,
      .domain = domain,
      .category = category,
      .name = name,
      .message = message,
      .severity = severity,
  };
}

// Complete synchronous sink contract. operator() consumes the non-owning event
// before returning; flush publishes any buffered events. Both operations are
// noexcept so observability failures never alter numerical control flow.
template <typename Sink>
concept TraceSink = requires(Sink& sink, const TraceEvent& event) {
  { sink(event) } noexcept -> std::same_as<void>;
  { sink.flush() } noexcept -> std::same_as<void>;
};

// Static no-op policy for template/generic code. Unlike AnyTraceSink's empty
// runtime state, this type lets the compiler remove both calls and storage
// (when held with [[no_unique_address]]) without virtual dispatch.
struct NullTraceSink {
  constexpr void operator()(const TraceEvent&) const noexcept {}
  constexpr void flush() const noexcept {}
};

static_assert(TraceSink<NullTraceSink>);

}  // namespace pemu::trace
