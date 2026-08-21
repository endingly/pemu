#pragma once

#include <pemu/trace/diag.hpp>

#include <concepts>
#include <cstdint>
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

using TraceValue =
    std::variant<bool, std::int64_t, std::uint64_t, double, std::string_view>;

struct TraceAttribute {
  std::string_view name;
  TraceValue value;
};

// TraceEvent and its attributes are non-owning views. A sink must consume them
// synchronously. A sink that retains events must copy names and string values.
struct TraceEvent {
  EventKind kind{EventKind::Trace};
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

template <typename Sink>
concept TraceSink = requires(Sink& sink, const TraceEvent& event) {
  { sink(event) } noexcept -> std::same_as<void>;
};

// The default sink is intentionally empty. When it is held with
// [[no_unique_address]], tracing adds neither storage nor virtual dispatch.
struct NullTraceSink {
  constexpr void operator()(const TraceEvent&) const noexcept {}
};

static_assert(TraceSink<NullTraceSink>);

}  // namespace pemu::trace
