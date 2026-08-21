#include <pemu/trace/renderer/statistics_event_renderer.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <iterator>
#include <ostream>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>

namespace pemu::trace::renderer {
namespace {

[[nodiscard]] const TraceAttribute* findAttribute(
    const TraceEvent& event, std::string_view name) noexcept {
  for (const auto& attribute : event.attributes) {
    if (attribute.name == name) {
      return &attribute;
    }
  }
  return nullptr;
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
void appendUnit(Buffer& output, units::precise_unit unit, bool brackets) {
  const auto unit_name =
      unit == units::precise::one ? std::string{"1"} : units::to_string(unit);
  if (brackets) {
    fmt::format_to(std::back_inserter(output), " [{}]", unit_name);
  } else {
    fmt::format_to(std::back_inserter(output), "{}", unit_name);
  }
}

template <typename Buffer>
void appendAttributeValue(Buffer& output, const TraceAttribute* attribute,
                          bool include_unit = true) {
  if (attribute == nullptr) {
    return;
  }
  appendValue(output, attribute->value);
  if (include_unit && attribute->unit.has_value()) {
    appendUnit(output, *attribute->unit, true);
  }
}

[[nodiscard]] std::string attributeValue(const TraceAttribute* attribute,
                                         bool include_unit = true) {
  fmt::memory_buffer output;
  appendAttributeValue(output, attribute, include_unit);
  return {output.data(), output.size()};
}

[[nodiscard]] std::string attributeUnit(const TraceAttribute* attribute) {
  if (attribute == nullptr || !attribute->unit.has_value()) {
    return {};
  }
  fmt::memory_buffer output;
  appendUnit(output, *attribute->unit, false);
  return {output.data(), output.size()};
}

[[nodiscard]] std::optional<std::uint64_t> unsignedAttribute(
    const TraceAttribute* attribute) noexcept {
  if (attribute == nullptr) {
    return std::nullopt;
  }
  if (const auto* value = std::get_if<std::uint64_t>(&attribute->value)) {
    return *value;
  }
  return std::nullopt;
}

[[nodiscard]] std::string statisticsField(const TraceEvent& event) {
  if (event.domain == DiagDomain::physics && event.category == "species") {
    const auto species = attributeValue(findAttribute(event, "species"));
    return species.empty() ? "species" : fmt::format("species[{}]", species);
  }
  if (event.domain == DiagDomain::physics && event.category == "charge") {
    return "charge_density";
  }
  return std::string{event.category};
}

[[nodiscard]] const TraceAttribute* statisticsMean(
    const TraceEvent& event) noexcept {
  if (const auto* mean = findAttribute(event, "volume_mean")) {
    return mean;
  }
  return findAttribute(event, "area_mean");
}

[[nodiscard]] const TraceAttribute* statisticsIntegral(
    const TraceEvent& event) noexcept {
  constexpr std::string_view names[]{"total_number", "net_charge",
                                     "volume_integral", "surface_integral"};
  for (const auto name : names) {
    if (const auto* integral = findAttribute(event, name)) {
      return integral;
    }
  }
  return nullptr;
}

}  // namespace

void StatisticsEventRenderer::render(std::ostream& stream,
                                     const TraceEvent& event) {
  const auto step = unsignedAttribute(findAttribute(event, "step"));
  if (group_.active && group_.step != step) {
    renderGroup(stream);
  }

  if (!group_.active) {
    group_.active = true;
    group_.step = step;
    group_.time = attributeValue(findAttribute(event, "time"));
  }

  const auto samples = unsignedAttribute(findAttribute(event, "samples"));
  if (event.category == "electric_field_normal") {
    group_.face_samples = samples;
  } else if (!group_.cell_samples.has_value()) {
    group_.cell_samples = samples;
  }

  if (group_.volume_semantics.empty()) {
    group_.volume_semantics =
        attributeValue(findAttribute(event, "volume_semantics"));
  }
  if (group_.physical_volume.empty()) {
    group_.physical_volume =
        attributeValue(findAttribute(event, "physical_volume"));
  }
  if (group_.physical_face_area.empty()) {
    group_.physical_face_area =
        attributeValue(findAttribute(event, "physical_face_area"));
  }

  const auto* minimum = findAttribute(event, "minimum");
  Row row{
      .field = statisticsField(event),
      .unit = attributeUnit(minimum),
      .minimum = attributeValue(minimum, false),
      .maximum = attributeValue(findAttribute(event, "maximum"), false),
      .mean = attributeValue(statisticsMean(event), false),
      .rms = attributeValue(findAttribute(event, "rms"), false),
      .integral = attributeValue(statisticsIntegral(event)),
      .severity = event.severity,
      .kind = event.kind,
      .non_finite =
          unsignedAttribute(findAttribute(event, "non_finite")).value_or(0),
      .negative =
          unsignedAttribute(findAttribute(event, "negative")).value_or(0),
      .invalid_weight =
          unsignedAttribute(findAttribute(event, "invalid_weight")).value_or(0),
  };
  group_.rows.push_back(std::move(row));
}

void StatisticsEventRenderer::flush(std::ostream& stream) {
  renderGroup(stream);
}

void StatisticsEventRenderer::renderGroup(std::ostream& stream) {
  if (!group_.active) {
    return;
  }

  fmt::memory_buffer output;
  if (has_output_) {
    fmt::format_to(std::back_inserter(output), "\n");
  }
  fmt::format_to(std::back_inserter(output), "STATISTICS");
  if (group_.step.has_value()) {
    fmt::format_to(std::back_inserter(output), " | STEP={}", *group_.step);
  }
  if (!group_.time.empty()) {
    fmt::format_to(std::back_inserter(output), " | TIME={}", group_.time);
  }
  fmt::format_to(std::back_inserter(output), "\n");

  if (!metadata_written_) {
    bool has_metadata = false;
    const auto appendMetadata = [&](std::string_view name, const auto& value) {
      if (!value.has_value()) {
        return;
      }
      fmt::format_to(std::back_inserter(output), "{}{}={}",
                     has_metadata ? " | " : "", name, *value);
      has_metadata = true;
    };
    appendMetadata("CELL_SAMPLES", group_.cell_samples);
    appendMetadata("FACE_SAMPLES", group_.face_samples);
    const auto appendStringMetadata = [&](std::string_view name,
                                          const std::string& value) {
      if (value.empty()) {
        return;
      }
      fmt::format_to(std::back_inserter(output), "{}{}={}",
                     has_metadata ? " | " : "", name, value);
      has_metadata = true;
    };
    appendStringMetadata("VOLUME_SEMANTICS", group_.volume_semantics);
    appendStringMetadata("PHYSICAL_VOLUME", group_.physical_volume);
    appendStringMetadata("PHYSICAL_FACE_AREA", group_.physical_face_area);
    if (has_metadata) {
      fmt::format_to(std::back_inserter(output), "\n");
      metadata_written_ = true;
    }
  }

  fmt::format_to(std::back_inserter(output),
                 "{:<30} | {:<16} | {:>14} | {:>14} | {:>14} | {:>14} | "
                 "{:>22}\n",
                 "FIELD", "UNIT", "MIN", "MAX", "MEAN", "RMS", "INTEGRAL");
  fmt::format_to(std::back_inserter(output),
                 "{:-<30}-+-{:-<16}-+-{:-<14}-+-{:-<14}-+-{:-<14}-+-"
                 "{:-<14}-+-{:-<22}\n",
                 "", "", "", "", "", "", "");
  for (const auto& row : group_.rows) {
    fmt::format_to(std::back_inserter(output),
                   "{:<30} | {:<16} | {:>14} | {:>14} | {:>14} | {:>14} | "
                   "{:>22}\n",
                   row.field, row.unit, row.minimum, row.maximum, row.mean,
                   row.rms, row.integral);
  }

  const bool has_diagnostics =
      std::ranges::any_of(group_.rows, [](const auto& row) {
        return row.non_finite != 0 || row.negative != 0 ||
               row.invalid_weight != 0;
      });
  if (has_diagnostics) {
    fmt::format_to(std::back_inserter(output), "\nSTATISTICS DIAGNOSTICS\n");
    fmt::format_to(std::back_inserter(output), "{:<8} | {:<12} | {:<30} | {}\n",
                   "LEVEL", "KIND", "FIELD", "CONTEXT");
    for (const auto& row : group_.rows) {
      if (row.non_finite == 0 && row.negative == 0 && row.invalid_weight == 0) {
        continue;
      }
      fmt::memory_buffer context;
      bool has_context = false;
      const auto appendCount = [&](std::string_view name, std::uint64_t value) {
        if (value == 0) {
          return;
        }
        fmt::format_to(std::back_inserter(context), "{}{}={}",
                       has_context ? "; " : "", name, value);
        has_context = true;
      };
      appendCount("non_finite", row.non_finite);
      appendCount("negative", row.negative);
      appendCount("invalid_weight", row.invalid_weight);
      fmt::format_to(std::back_inserter(output), "{:<8} | {:<12} | {:<30} | ",
                     pemu::to_string(row.severity), pemu::to_string(row.kind),
                     row.field);
      output.append(context.data(), context.data() + context.size());
      fmt::format_to(std::back_inserter(output), "\n");
    }
  }

  stream.write(output.data(), static_cast<std::streamsize>(output.size()));
  group_ = {};
  has_output_ = true;
}

}  // namespace pemu::trace::renderer
