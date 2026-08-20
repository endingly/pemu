#pragma once

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

using TraceValue =
    std::variant<bool, std::int64_t, std::uint64_t, double, std::string_view>;

struct TraceAttribute {
  std::string_view name;
  TraceValue value;
};

// TraceEvent and its attributes are non-owning views. A sink must consume them
// synchronously. A sink that retains events must copy names and string values.
struct TraceEvent {
  std::string_view category;
  std::string_view name;
  Severity severity{Severity::Trace};
  std::span<const TraceAttribute> attributes{};
};

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
