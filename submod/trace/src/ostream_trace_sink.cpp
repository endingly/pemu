#include <magic_enum/magic_enum.hpp>
#include <pemu/trace/ostream_trace_sink.hpp>

#include <fmt/format.h>

#include <iterator>
#include <ostream>
#include <string_view>
#include <type_traits>
#include <variant>

namespace pemu::trace {
namespace {

[[nodiscard]] constexpr std::string_view severityName(
    Severity severity) noexcept {
  return magic_enum::enum_name(severity);
}

}  // namespace

OstreamTraceSink::OstreamTraceSink(std::ostream& stream) noexcept
    : stream_(&stream) {}

void OstreamTraceSink::operator()(const TraceEvent& event) noexcept {
  try {
    fmt::memory_buffer output;
    fmt::format_to(std::back_inserter(output), "[{}] {}",
                   severityName(event.severity), diagDomainName(event.domain));

    if (!event.category.empty()) {
      fmt::format_to(std::back_inserter(output), ".{}", event.category);
    }
    if (!event.name.empty()) {
      fmt::format_to(std::back_inserter(output), ".{}", event.name);
    }
    if (!event.message.empty()) {
      fmt::format_to(std::back_inserter(output), ": {}", event.message);
    }

    for (const auto& attribute : event.attributes) {
      std::visit(
          [&output, name = attribute.name](const auto& value) {
            using Value = std::remove_cvref_t<decltype(value)>;
            if constexpr (std::is_same_v<Value, bool>) {
              fmt::format_to(std::back_inserter(output), " {}={:d}", name,
                             value);
            } else {
              fmt::format_to(std::back_inserter(output), " {}={}", name, value);
            }
          },
          attribute.value);
    }

    fmt::format_to(std::back_inserter(output), "\n");
    stream_->write(output.data(), static_cast<std::streamsize>(output.size()));
    failed_ = failed_ || !stream_->good();
  } catch (...) {
    failed_ = true;
  }
}

bool OstreamTraceSink::failed() const noexcept {
  return failed_;
}

}  // namespace pemu::trace
