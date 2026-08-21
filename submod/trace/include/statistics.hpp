#pragma once

#include <pemu/unit/mp_units_bridge.hpp>

#include <mp-units/framework/quantity.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace pemu::trace {

// A 2D finite-volume mesh stores cell areas and face lengths.  Diagnostics
// turn them into physical measures by interpreting the mesh as a planar slab
// with a caller-supplied out-of-plane depth.  A 3D mesh already stores physical
// cell volumes and face areas, so the depth is ignored.
struct PhysicalVolumeSemantics {
  int mesh_dimension{};
  double planar_depth{1.0};
  units::precise_unit mesh_length_unit{units::precise::one};

  [[nodiscard]] constexpr bool valid() const noexcept {
    return numericSemanticsValid() &&
           mesh_length_unit.is_convertible(units::precise::m);
  }

  [[nodiscard]] double cellVolume(double mesh_cell_measure) const noexcept {
    return effectiveMeasure(mesh_cell_measure);
  }

  [[nodiscard]] double faceArea(double mesh_face_measure) const noexcept {
    return effectiveMeasure(mesh_face_measure);
  }

  [[nodiscard]] constexpr const char* name() const noexcept {
    return mesh_dimension == 2 ? "planar_extrusion" : "native_3d";
  }

  [[nodiscard]] constexpr units::precise_unit cellWeightUnit() const noexcept {
    return mesh_length_unit.pow(3);
  }

  [[nodiscard]] constexpr units::precise_unit faceWeightUnit() const noexcept {
    return mesh_length_unit.pow(2);
  }

 private:
  [[nodiscard]] constexpr bool numericSemanticsValid() const noexcept {
    return (mesh_dimension == 2 && planar_depth > 0.0) || mesh_dimension == 3;
  }

  [[nodiscard]] double effectiveMeasure(double mesh_measure) const noexcept {
    // Deliberately numeric-only: this method is called once per cell/face.
    // Runtime-unit validation and algebra stay outside the mesh scan.
    if (!numericSemanticsValid() || !std::isfinite(mesh_measure) ||
        mesh_measure <= 0.0) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return mesh_dimension == 2 ? mesh_measure * planar_depth : mesh_measure;
  }
};

class StatisticsOptions {
 public:
  constexpr StatisticsOptions() = default;

  template <mp_units::Quantity Depth, mp_units::Reference MeshLengthReference>
  constexpr StatisticsOptions(bool enabled, std::size_t sample_every_steps,
                              Depth planar_depth,
                              MeshLengthReference mesh_length_reference)
      : enabled_(enabled),
        sample_every_steps_(sample_every_steps),
        planar_depth_(static_cast<double>(planar_depth.numerical_value_in(
            mp_units::get_unit(mesh_length_reference)))),
        mesh_length_unit_(pemu::unit::bridgeUnit(mesh_length_reference)) {
    static_assert(pemu::unit::bridgeReference(MeshLengthReference{}).kind() ==
                  pemu::unit::QuantityKind::length);
  }

  [[nodiscard]] constexpr bool valid() const noexcept {
    return !enabled_ || (sample_every_steps_ > 0 && planar_depth_ > 0.0 &&
                         mesh_length_unit_.is_convertible(units::precise::m));
  }

  [[nodiscard]] constexpr bool shouldSample(std::size_t step) const noexcept {
    return enabled_ && valid() && step % sample_every_steps_ == 0;
  }

  [[nodiscard]] constexpr bool enabled() const noexcept { return enabled_; }

  [[nodiscard]] constexpr std::size_t sampleEverySteps() const noexcept {
    return sample_every_steps_;
  }

  [[nodiscard]] constexpr double planarDepth() const noexcept {
    return planar_depth_;
  }

  [[nodiscard]] constexpr units::precise_unit meshLengthUnit() const noexcept {
    return mesh_length_unit_;
  }

 private:
  bool enabled_{false};
  std::size_t sample_every_steps_{1};
  double planar_depth_{1.0};
  units::precise_unit mesh_length_unit_{units::precise::one};
};

struct ScalarFieldStatistics {
  std::uint64_t sample_count{};
  std::uint64_t finite_count{};
  std::uint64_t non_finite_count{};
  std::uint64_t negative_count{};
  std::uint64_t invalid_weight_count{};

  double minimum{};
  double maximum{};
  double max_abs{};
  double weight_sum{};
  double integral{};
  double l1_integral{};
  double weighted_mean{};
  double weighted_rms{};

  units::precise_unit value_unit{units::precise::one};
  units::precise_unit weight_unit{units::precise::one};

  [[nodiscard]] constexpr units::precise_unit integralUnit() const noexcept {
    return value_unit * weight_unit;
  }

  [[nodiscard]] bool valid() const noexcept {
    return sample_count != 0 && finite_count != 0 && non_finite_count == 0 &&
           invalid_weight_count == 0 && weight_sum > 0.0 &&
           std::isfinite(minimum) && std::isfinite(maximum) &&
           std::isfinite(max_abs) && std::isfinite(weight_sum) &&
           std::isfinite(integral) && std::isfinite(l1_integral) &&
           std::isfinite(weighted_mean) && std::isfinite(weighted_rms);
  }
};

class ScalarFieldStatisticsAccumulator {
 public:
  void add(double value, double physical_weight = 1.0) noexcept {
    ++statistics_.sample_count;

    if (!std::isfinite(value)) {
      ++statistics_.non_finite_count;
      return;
    }
    ++statistics_.finite_count;
    if (value < 0.0) {
      ++statistics_.negative_count;
    }

    if (!has_finite_value_) {
      statistics_.minimum = value;
      statistics_.maximum = value;
      statistics_.max_abs = std::abs(value);
      has_finite_value_ = true;
    } else {
      statistics_.minimum = std::min(statistics_.minimum, value);
      statistics_.maximum = std::max(statistics_.maximum, value);
      statistics_.max_abs = std::max(statistics_.max_abs, std::abs(value));
    }

    if (!std::isfinite(physical_weight) || physical_weight <= 0.0) {
      ++statistics_.invalid_weight_count;
      return;
    }

    weight_sum_.add(physical_weight);
    integral_.add(value * physical_weight);
    l1_integral_.add(std::abs(value) * physical_weight);
    square_integral_.add(value * value * physical_weight);
  }

  [[nodiscard]] ScalarFieldStatistics finish(
      units::precise_unit value_unit = units::precise::one,
      units::precise_unit weight_unit = units::precise::one) const noexcept {
    auto result = statistics_;
    result.value_unit = value_unit;
    result.weight_unit = weight_unit;
    result.weight_sum = weight_sum_.value();
    result.integral = integral_.value();
    result.l1_integral = l1_integral_.value();
    if (result.weight_sum > 0.0) {
      result.weighted_mean = result.integral / result.weight_sum;
      result.weighted_rms =
          std::sqrt(square_integral_.value() / result.weight_sum);
    }
    return result;
  }

 private:
  class CompensatedSum {
   public:
    void add(double value) noexcept {
      const double corrected = value - correction_;
      const double next = sum_ + corrected;
      correction_ = (next - sum_) - corrected;
      sum_ = next;
    }

    [[nodiscard]] double value() const noexcept { return sum_; }

   private:
    double sum_{};
    double correction_{};
  };

  ScalarFieldStatistics statistics_{};
  CompensatedSum weight_sum_{};
  CompensatedSum integral_{};
  CompensatedSum l1_integral_{};
  CompensatedSum square_integral_{};
  bool has_finite_value_{false};
};

}  // namespace pemu::trace
