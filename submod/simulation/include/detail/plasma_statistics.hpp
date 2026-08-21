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
  return accumulator.finish();
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
  return accumulator.finish();
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
    const std::array attributes{
        trace::TraceAttribute{"step", static_cast<std::uint64_t>(step)},
        trace::TraceAttribute{"time", time},
        trace::TraceAttribute{"species_id", static_cast<std::uint64_t>(i)},
        trace::TraceAttribute{"species", std::string_view{species.at(id).name}},
        trace::TraceAttribute{"samples", statistics.sample_count},
        trace::TraceAttribute{"non_finite", statistics.non_finite_count},
        trace::TraceAttribute{"negative", statistics.negative_count},
        trace::TraceAttribute{"minimum", statistics.minimum},
        trace::TraceAttribute{"maximum", statistics.maximum},
        trace::TraceAttribute{"max_abs", statistics.max_abs},
        trace::TraceAttribute{"physical_volume", statistics.weight_sum},
        trace::TraceAttribute{"volume_mean", statistics.weighted_mean},
        trace::TraceAttribute{"total_number", statistics.integral},
        trace::TraceAttribute{"rms", statistics.weighted_rms},
        trace::TraceAttribute{"volume_semantics",
                              std::string_view{physical_volume.name()}},
    };
    sink({.kind = statisticsKind(severity),
          .domain = trace::DiagDomain::physics,
          .category = "species",
          .name = "statistics",
          .severity = severity,
          .attributes = attributes});
  }
}

template <trace::TraceSink Sink>
void emitChargeStatistics(
    Sink& sink, const field::CellField<double>& charge_density,
    const trace::PhysicalVolumeSemantics& physical_volume, std::size_t step,
    double time) noexcept {
  const auto statistics = cellStatistics(charge_density, physical_volume);
  const auto severity = statisticsSeverity(statistics, false);
  const double relative_imbalance =
      statistics.l1_integral > 0.0
          ? std::abs(statistics.integral) / statistics.l1_integral
          : 0.0;
  const std::array attributes{
      trace::TraceAttribute{"step", static_cast<std::uint64_t>(step)},
      trace::TraceAttribute{"time", time},
      trace::TraceAttribute{"samples", statistics.sample_count},
      trace::TraceAttribute{"non_finite", statistics.non_finite_count},
      trace::TraceAttribute{"minimum", statistics.minimum},
      trace::TraceAttribute{"maximum", statistics.maximum},
      trace::TraceAttribute{"max_abs", statistics.max_abs},
      trace::TraceAttribute{"physical_volume", statistics.weight_sum},
      trace::TraceAttribute{"volume_mean", statistics.weighted_mean},
      trace::TraceAttribute{"net_charge", statistics.integral},
      trace::TraceAttribute{"absolute_charge", statistics.l1_integral},
      trace::TraceAttribute{"relative_imbalance", relative_imbalance},
      trace::TraceAttribute{"volume_semantics",
                            std::string_view{physical_volume.name()}},
  };
  sink({.kind = statisticsKind(severity),
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
  const std::array attributes{
      trace::TraceAttribute{"step", static_cast<std::uint64_t>(step)},
      trace::TraceAttribute{"time", time},
      trace::TraceAttribute{"samples", statistics.sample_count},
      trace::TraceAttribute{"non_finite", statistics.non_finite_count},
      trace::TraceAttribute{"minimum", statistics.minimum},
      trace::TraceAttribute{"maximum", statistics.maximum},
      trace::TraceAttribute{"max_abs", statistics.max_abs},
      trace::TraceAttribute{"physical_volume", statistics.weight_sum},
      trace::TraceAttribute{"volume_mean", statistics.weighted_mean},
      trace::TraceAttribute{"rms", statistics.weighted_rms},
      trace::TraceAttribute{"volume_semantics",
                            std::string_view{physical_volume.name()}},
  };
  sink({.kind = statisticsKind(severity),
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
  const std::array attributes{
      trace::TraceAttribute{"step", static_cast<std::uint64_t>(step)},
      trace::TraceAttribute{"time", time},
      trace::TraceAttribute{"samples", statistics.sample_count},
      trace::TraceAttribute{"non_finite", statistics.non_finite_count},
      trace::TraceAttribute{"minimum", statistics.minimum},
      trace::TraceAttribute{"maximum", statistics.maximum},
      trace::TraceAttribute{"max_abs", statistics.max_abs},
      trace::TraceAttribute{"physical_face_area", statistics.weight_sum},
      trace::TraceAttribute{"area_mean", statistics.weighted_mean},
      trace::TraceAttribute{"rms", statistics.weighted_rms},
      trace::TraceAttribute{"volume_semantics",
                            std::string_view{physical_volume.name()}},
  };
  sink({.kind = statisticsKind(severity),
        .domain = trace::DiagDomain::field,
        .category = "electric_field_normal",
        .name = "statistics",
        .severity = severity,
        .attributes = attributes});
}

template <trace::TraceSink Sink>
void emitPlasmaStatistics(
    Sink& sink, const physics::SpeciesSet& species,
    const physics::SpeciesCellFields& density,
    const field::CellField<double>& charge_density,
    const field::CellField<double>& potential,
    const field::FaceField<double>& electric_field_normal,
    const trace::StatisticsOptions& options, std::size_t step,
    double time) noexcept {
  if (!options.shouldSample(step)) {
    return;
  }

  const trace::PhysicalVolumeSemantics physical_volume{
      .mesh_dimension = density.mesh().dimension(),
      .planar_depth = options.planar_depth,
  };
  emitSpeciesStatistics(sink, species, density, physical_volume, step, time);
  emitChargeStatistics(sink, charge_density, physical_volume, step, time);
  emitPotentialStatistics(sink, potential, physical_volume, step, time);
  emitElectricFieldStatistics(sink, electric_field_normal, physical_volume,
                              step, time);
}

}  // namespace pemu::simulation::detail
