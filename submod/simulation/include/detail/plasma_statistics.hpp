#pragma once

#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>
#include <pemu/physics/species.hpp>
#include <pemu/trace/statistics.hpp>
#include <pemu/trace/trace.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string_view>

namespace pemu::simulation::detail {

[[nodiscard]] inline trace::ScalarFieldStatistics cellStatistics(
    const field::CellField<double>& values,
    const trace::PhysicalVolumeSemantics& physical_volume) noexcept {
  trace::ScalarFieldStatisticsAccumulator accumulator;
  const auto& mesh = values.mesh();
  for (mesh::CellId cell = 0; cell < mesh.numCells(); ++cell) {
    accumulator.add(values[cell],
                    physical_volume.cellVolume(mesh.cellVolume(cell)));
  }
  return accumulator.finish(values.metadata().physical_quantity->unit(),
                            physical_volume.cellWeightUnit());
}

[[nodiscard]] inline trace::ScalarFieldStatistics faceStatistics(
    const field::FaceField<double>& values,
    const trace::PhysicalVolumeSemantics& physical_volume) noexcept {
  trace::ScalarFieldStatisticsAccumulator accumulator;
  const auto& mesh = values.mesh();
  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {
    accumulator.add(values[face],
                    physical_volume.faceArea(mesh.faceArea(face)));
  }
  return accumulator.finish(values.metadata().physical_quantity->unit(),
                            physical_volume.faceWeightUnit());
}

template <typename Field>
[[nodiscard]] inline bool hasQuantityKind(
    const Field& values, unit::QuantityKind expected_kind) noexcept {
  return values.metadata().physical_quantity.has_value() &&
         values.metadata().physical_quantity->kind() == expected_kind;
}

inline void validatePlasmaStatisticsConfiguration(
    const trace::StatisticsOptions& options,
    const physics::SpeciesCellFields& density,
    const field::CellField<double>& charge_density,
    const field::CellField<double>& potential,
    const field::FaceField<double>& electric_field_normal) {
  if (!options.enabled()) {
    return;
  }
  if (!options.valid() ||
      (density.mesh().dimension() != 2 && density.mesh().dimension() != 3)) {
    throw std::invalid_argument("invalid plasma statistics configuration");
  }
  for (const auto& species_density : density) {
    if (!hasQuantityKind(species_density,
                         unit::QuantityKind::particle_number_density)) {
      throw std::invalid_argument(
          "species statistics require particle-number-density metadata");
    }
  }
  if (!hasQuantityKind(charge_density,
                       unit::QuantityKind::electric_charge_density) ||
      !hasQuantityKind(potential, unit::QuantityKind::electric_potential) ||
      !hasQuantityKind(electric_field_normal,
                       unit::QuantityKind::normal_electric_field_strength)) {
    throw std::invalid_argument(
        "plasma statistics require physical field metadata");
  }
}

[[nodiscard]] inline trace::Severity statisticsSeverity(
    const trace::ScalarFieldStatistics& statistics,
    bool negative_values_are_invalid) noexcept {
  if (statistics.non_finite_count != 0 ||
      statistics.invalid_weight_count != 0 || !statistics.valid()) {
    return trace::Severity::Error;
  }
  if (negative_values_are_invalid && statistics.negative_count != 0) {
    return trace::Severity::Warning;
  }
  return trace::Severity::Debug;
}

[[nodiscard]] inline trace::EventKind statisticsKind(
    trace::Severity severity) noexcept {
  return severity == trace::Severity::Warning ||
                 severity == trace::Severity::Error ||
                 severity == trace::Severity::Critical
             ? trace::EventKind::Diagnostic
             : trace::EventKind::Trace;
}

template <trace::TraceSink Sink>
void emitSpeciesStatistics(
    Sink& sink, const physics::SpeciesSet& species,
    const physics::SpeciesCellFields& density,
    const trace::PhysicalVolumeSemantics& physical_volume, std::size_t step,
    double time) noexcept {
  for (std::size_t i = 0; i < species.size(); ++i) {
    const physics::SpeciesId id{static_cast<std::uint32_t>(i)};
    const auto statistics = cellStatistics(density[id], physical_volume);
    const auto severity = statisticsSeverity(statistics, true);
    const auto value_unit = statistics.value_unit;
    const auto weight_unit = statistics.weight_unit;
    const auto integral_unit = statistics.integralUnit();
    const std::array attributes{
        trace::TraceAttribute{"step", static_cast<std::uint64_t>(step)},
        trace::TraceAttribute{"time", time, units::precise::s},
        trace::TraceAttribute{"species_id", static_cast<std::uint64_t>(i)},
        trace::TraceAttribute{"species", std::string_view{species.at(id).name}},
        trace::TraceAttribute{"samples", statistics.sample_count},
        trace::TraceAttribute{"non_finite", statistics.non_finite_count},
        trace::TraceAttribute{"negative", statistics.negative_count},
        trace::TraceAttribute{"minimum", statistics.minimum, value_unit},
        trace::TraceAttribute{"maximum", statistics.maximum, value_unit},
        trace::TraceAttribute{"max_abs", statistics.max_abs, value_unit},
        trace::TraceAttribute{"physical_volume", statistics.weight_sum,
                              weight_unit},
        trace::TraceAttribute{"volume_mean", statistics.weighted_mean,
                              value_unit},
        trace::TraceAttribute{"total_number", statistics.integral,
                              integral_unit},
        trace::TraceAttribute{"l1_total_number", statistics.l1_integral,
                              integral_unit},
        trace::TraceAttribute{"rms", statistics.weighted_rms, value_unit},
        trace::TraceAttribute{"volume_semantics",
                              std::string_view{physical_volume.name()}},
    };
    sink({.kind = statisticsKind(severity),
          .output_channel = trace::OutputChannel::Statistics,
          .domain = trace::DiagDomain::physics,
          .category = "species",
          .name = "statistics",
          .severity = severity,
          .attributes = attributes});
  }
}

template <trace::TraceSink Sink>
void emitChargeStatistics(Sink& sink,
                          const field::CellField<double>& charge_density,
                          const trace::PhysicalVolumeSemantics& physical_volume,
                          std::size_t step, double time) noexcept {
  const auto statistics = cellStatistics(charge_density, physical_volume);
  const auto severity = statisticsSeverity(statistics, false);
  const double relative_imbalance =
      statistics.l1_integral > 0.0
          ? std::abs(statistics.integral) / statistics.l1_integral
          : 0.0;
  const auto value_unit = statistics.value_unit;
  const auto weight_unit = statistics.weight_unit;
  const auto integral_unit = statistics.integralUnit();
  const std::array attributes{
      trace::TraceAttribute{"step", static_cast<std::uint64_t>(step)},
      trace::TraceAttribute{"time", time, units::precise::s},
      trace::TraceAttribute{"samples", statistics.sample_count},
      trace::TraceAttribute{"non_finite", statistics.non_finite_count},
      trace::TraceAttribute{"minimum", statistics.minimum, value_unit},
      trace::TraceAttribute{"maximum", statistics.maximum, value_unit},
      trace::TraceAttribute{"max_abs", statistics.max_abs, value_unit},
      trace::TraceAttribute{"physical_volume", statistics.weight_sum,
                            weight_unit},
      trace::TraceAttribute{"volume_mean", statistics.weighted_mean,
                            value_unit},
      trace::TraceAttribute{"net_charge", statistics.integral, integral_unit},
      trace::TraceAttribute{"absolute_charge", statistics.l1_integral,
                            integral_unit},
      trace::TraceAttribute{"rms", statistics.weighted_rms, value_unit},
      trace::TraceAttribute{"relative_imbalance", relative_imbalance,
                            units::precise::one},
      trace::TraceAttribute{"volume_semantics",
                            std::string_view{physical_volume.name()}},
  };
  sink({.kind = statisticsKind(severity),
        .output_channel = trace::OutputChannel::Statistics,
        .domain = trace::DiagDomain::physics,
        .category = "charge",
        .name = "statistics",
        .severity = severity,
        .attributes = attributes});
}

template <trace::TraceSink Sink>
void emitPotentialStatistics(
    Sink& sink, const field::CellField<double>& potential,
    const trace::PhysicalVolumeSemantics& physical_volume, std::size_t step,
    double time) noexcept {
  const auto statistics = cellStatistics(potential, physical_volume);
  const auto severity = statisticsSeverity(statistics, false);
  const auto value_unit = statistics.value_unit;
  const auto weight_unit = statistics.weight_unit;
  const auto integral_unit = statistics.integralUnit();
  const std::array attributes{
      trace::TraceAttribute{"step", static_cast<std::uint64_t>(step)},
      trace::TraceAttribute{"time", time, units::precise::s},
      trace::TraceAttribute{"samples", statistics.sample_count},
      trace::TraceAttribute{"non_finite", statistics.non_finite_count},
      trace::TraceAttribute{"minimum", statistics.minimum, value_unit},
      trace::TraceAttribute{"maximum", statistics.maximum, value_unit},
      trace::TraceAttribute{"max_abs", statistics.max_abs, value_unit},
      trace::TraceAttribute{"physical_volume", statistics.weight_sum,
                            weight_unit},
      trace::TraceAttribute{"volume_mean", statistics.weighted_mean,
                            value_unit},
      trace::TraceAttribute{"volume_integral", statistics.integral,
                            integral_unit},
      trace::TraceAttribute{"l1_volume_integral", statistics.l1_integral,
                            integral_unit},
      trace::TraceAttribute{"rms", statistics.weighted_rms, value_unit},
      trace::TraceAttribute{"volume_semantics",
                            std::string_view{physical_volume.name()}},
  };
  sink({.kind = statisticsKind(severity),
        .output_channel = trace::OutputChannel::Statistics,
        .domain = trace::DiagDomain::field,
        .category = "potential",
        .name = "statistics",
        .severity = severity,
        .attributes = attributes});
}

template <trace::TraceSink Sink>
void emitElectricFieldStatistics(
    Sink& sink, const field::FaceField<double>& electric_field_normal,
    const trace::PhysicalVolumeSemantics& physical_volume, std::size_t step,
    double time) noexcept {
  const auto statistics =
      faceStatistics(electric_field_normal, physical_volume);
  const auto severity = statisticsSeverity(statistics, false);
  const auto value_unit = statistics.value_unit;
  const auto weight_unit = statistics.weight_unit;
  const auto integral_unit = statistics.integralUnit();
  const std::array attributes{
      trace::TraceAttribute{"step", static_cast<std::uint64_t>(step)},
      trace::TraceAttribute{"time", time, units::precise::s},
      trace::TraceAttribute{"samples", statistics.sample_count},
      trace::TraceAttribute{"non_finite", statistics.non_finite_count},
      trace::TraceAttribute{"minimum", statistics.minimum, value_unit},
      trace::TraceAttribute{"maximum", statistics.maximum, value_unit},
      trace::TraceAttribute{"max_abs", statistics.max_abs, value_unit},
      trace::TraceAttribute{"physical_face_area", statistics.weight_sum,
                            weight_unit},
      trace::TraceAttribute{"area_mean", statistics.weighted_mean, value_unit},
      trace::TraceAttribute{"surface_integral", statistics.integral,
                            integral_unit},
      trace::TraceAttribute{"l1_surface_integral", statistics.l1_integral,
                            integral_unit},
      trace::TraceAttribute{"rms", statistics.weighted_rms, value_unit},
      trace::TraceAttribute{"volume_semantics",
                            std::string_view{physical_volume.name()}},
  };
  sink({.kind = statisticsKind(severity),
        .output_channel = trace::OutputChannel::Statistics,
        .domain = trace::DiagDomain::field,
        .category = "electric_field_normal",
        .name = "statistics",
        .severity = severity,
        .attributes = attributes});
}

template <trace::TraceSink Sink>
void emitPlasmaStatistics(Sink& sink, const physics::SpeciesSet& species,
                          const physics::SpeciesCellFields& density,
                          const field::CellField<double>& charge_density,
                          const field::CellField<double>& potential,
                          const field::FaceField<double>& electric_field_normal,
                          const trace::StatisticsOptions& options,
                          std::size_t step, double time) noexcept {
  if (!options.shouldSample(step)) {
    return;
  }

  const trace::PhysicalVolumeSemantics physical_volume{
      .mesh_dimension = density.mesh().dimension(),
      .planar_depth = options.planarDepth(),
      .mesh_length_unit = options.meshLengthUnit(),
  };
  emitSpeciesStatistics(sink, species, density, physical_volume, step, time);
  emitChargeStatistics(sink, charge_density, physical_volume, step, time);
  emitPotentialStatistics(sink, potential, physical_volume, step, time);
  emitElectricFieldStatistics(sink, electric_field_normal, physical_volume,
                              step, time);
}

}  // namespace pemu::simulation::detail
