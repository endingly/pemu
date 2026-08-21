#pragma once

#include <pemu/trace/trace.hpp>

#include <fmt/format.h>

#include <iterator>
#include <ostream>
#include <string_view>
#include <type_traits>
#include <variant>

namespace pemu::trace::renderer::detail {

[[nodiscard]] inline const TraceAttribute* findAttribute(
    const TraceEvent& event, std::string_view name) noexcept {
  for (const auto& attribute : event.attributes) {
    if (attribute.name == name) {
      return &attribute;
    }
  }
  return nullptr;
}

[[nodiscard]] constexpr bool isTableColumn(std::string_view name) noexcept {
  return name == "step" || name == "time" || name == "dt" ||
         name == "end_time" || name == "solver_status" ||
         name == "residual_norm" || name == "relative_residual";
}

template <typename Buffer>
void appendValue(Buffer& output, const TraceValue& value) {
  std::visit(
      [&output](const auto& item) {
        using Value = std::remove_cvref_t<decltype(item)>;
        if constexpr (std::is_same_v<Value, bool>) {
          fmt::format_to(std::back_inserter(output), "{}",
                         item ? "true" : "false");
        } else if constexpr (std::is_same_v<Value, double>) {
          fmt::format_to(std::back_inserter(output), "{:.6g}", item);
        } else {
          fmt::format_to(std::back_inserter(output), "{}", item);
        }
      },
      value);
}

template <typename Buffer>
void appendUnit(Buffer& output, units::precise_unit unit) {
  if (unit == units::precise::one) {
    fmt::format_to(std::back_inserter(output), " [1]");
  } else {
    fmt::format_to(std::back_inserter(output), " [{}]", units::to_string(unit));
  }
}

template <typename Buffer>
void appendAttributeValue(Buffer& output, const TraceAttribute* attribute) {
  if (attribute == nullptr) {
    return;
  }
  appendValue(output, attribute->value);
  if (attribute->unit.has_value()) {
    appendUnit(output, *attribute->unit);
  }
}

template <typename Buffer>
void appendEventName(Buffer& output, const TraceEvent& event) {
  fmt::format_to(std::back_inserter(output), "{}",
                 pemu::to_string(event.domain));
  if (!event.category.empty()) {
    fmt::format_to(std::back_inserter(output), ".{}", event.category);
  }
  if (!event.name.empty()) {
    fmt::format_to(std::back_inserter(output), ".{}", event.name);
  }
}

template <typename Buffer>
void appendDetails(Buffer& output, const TraceEvent& event) {
  bool has_detail = false;
  if (!event.message.empty()) {
    fmt::format_to(std::back_inserter(output), "{}", event.message);
    has_detail = true;
  }
  for (const auto& attribute : event.attributes) {
    if (isTableColumn(attribute.name)) {
      continue;
    }
    fmt::format_to(std::back_inserter(output), "{}{}=", has_detail ? "; " : "",
                   attribute.name);
    appendValue(output, attribute.value);
    if (attribute.unit.has_value()) {
      appendUnit(output, *attribute.unit);
    }
    has_detail = true;
  }
}

[[nodiscard]] inline fmt::string_view asStringView(
    const fmt::memory_buffer& buffer) {
  return {buffer.data(), buffer.size()};
}

inline void renderEventTable(std::ostream& stream, const TraceEvent& event,
                             bool write_header) {
  fmt::memory_buffer output;
  if (write_header) {
    fmt::format_to(
        std::back_inserter(output),
        "{:<8} | {:<58} | {:>4} | {:>15} | {:>15} | {:>15} | {:>24} | "
        "{:>12} | {:>12} | {}\n",
        "LEVEL", "EVENT", "STEP", "TIME", "DT", "END TIME", "SOLVER",
        "RESIDUAL", "REL RESIDUAL", "DETAILS");
    fmt::format_to(std::back_inserter(output),
                   "{:-<8}-+-{:-<58}-+-{:-<4}-+-{:-<15}-+-{:-<15}-+-{:-<15}-+-"
                   "{:-<24}-+-{:-<12}-+-{:-<12}-+-{}\n",
                   "", "", "", "", "", "", "", "", "", "");
  }

  fmt::memory_buffer event_name;
  fmt::memory_buffer step;
  fmt::memory_buffer time;
  fmt::memory_buffer dt;
  fmt::memory_buffer end_time;
  fmt::memory_buffer solver_status;
  fmt::memory_buffer residual_norm;
  fmt::memory_buffer relative_residual;
  fmt::memory_buffer details;
  appendEventName(event_name, event);
  appendAttributeValue(step, findAttribute(event, "step"));
  appendAttributeValue(time, findAttribute(event, "time"));
  appendAttributeValue(dt, findAttribute(event, "dt"));
  appendAttributeValue(end_time, findAttribute(event, "end_time"));
  appendAttributeValue(solver_status, findAttribute(event, "solver_status"));
  appendAttributeValue(residual_norm, findAttribute(event, "residual_norm"));
  appendAttributeValue(relative_residual,
                       findAttribute(event, "relative_residual"));
  appendDetails(details, event);

  fmt::format_to(std::back_inserter(output),
                 "{:<8} | {:<58} | {:>4} | {:>15} | {:>15} | {:>15} | {:>24} | "
                 "{:>12} | {:>12} | {}\n",
                 pemu::to_string(event.severity), asStringView(event_name),
                 asStringView(step), asStringView(time), asStringView(dt),
                 asStringView(end_time), asStringView(solver_status),
                 asStringView(residual_norm), asStringView(relative_residual),
                 asStringView(details));
  stream.write(output.data(), static_cast<std::streamsize>(output.size()));
}

}  // namespace pemu::trace::renderer::detail
