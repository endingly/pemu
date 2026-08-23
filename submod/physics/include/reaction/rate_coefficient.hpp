#pragma once

#include <pemu/field/cell_field.hpp>
#include <pemu/field/utils.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pemu::physics::reaction {

/** @brief Selects interpolation in coefficient or logarithmic coefficient. */
enum class RateInterpolation {
  linear,
  log_linear,
};

/** @brief Selects behavior outside the tabulated coordinate interval. */
enum class RateTableBounds {
  error,
  clamp,
};

/**
 * @brief Immutable one-dimensional lookup table for a reaction coefficient.
 *
 * Coordinates are unit-agnostic numerical values chosen by the caller, most
 * commonly E/N in Td or Te in eV. Coefficients must use the same volume/time
 * unit expected by the density and reaction-rate fields. The table never
 * performs linear extrapolation because that can silently create negative or
 * rapidly unphysical chemistry rates.
 */
class TabulatedRateCoefficient {
 public:
  /** @brief Validates and stores one monotone coefficient table. */
  TabulatedRateCoefficient(
      std::vector<double> coordinates, std::vector<double> coefficients,
      RateInterpolation interpolation = RateInterpolation::log_linear,
      RateTableBounds bounds = RateTableBounds::error)
      : coordinates_(std::move(coordinates)),
        coefficients_(std::move(coefficients)),
        interpolation_(interpolation),
        bounds_(bounds) {
    validateTable();
  }

  /** @brief Returns the interpolated non-negative coefficient at coordinate. */
  [[nodiscard]] double operator()(double coordinate) const {
    validateCoordinate(coordinate);

    if (coordinate <= coordinates_.front()) {
      return coefficients_.front();
    }
    if (coordinate >= coordinates_.back()) {
      return coefficients_.back();
    }

    const auto upper =
        std::upper_bound(coordinates_.begin(), coordinates_.end(), coordinate);
    const auto upper_index =
        static_cast<std::size_t>(upper - coordinates_.begin());
    const std::size_t lower_index = upper_index - 1;
    const double fraction =
        (coordinate - coordinates_[lower_index]) /
        (coordinates_[upper_index] - coordinates_[lower_index]);

    if (interpolation_ == RateInterpolation::linear) {
      return coefficients_[lower_index] +
             fraction *
                 (coefficients_[upper_index] - coefficients_[lower_index]);
    }
    return std::exp(std::log(coefficients_[lower_index]) +
                    fraction * (std::log(coefficients_[upper_index]) -
                                std::log(coefficients_[lower_index])));
  }

  /** @brief Checks whether a coordinate can be evaluated under this policy. */
  void validateCoordinate(double coordinate) const {
    if (!std::isfinite(coordinate)) {
      throw std::invalid_argument("rate-table coordinate must be finite");
    }
    if (bounds_ == RateTableBounds::error &&
        (coordinate < coordinates_.front() ||
         coordinate > coordinates_.back())) {
      throw std::out_of_range("rate-table coordinate is outside table bounds");
    }
  }

  /** @brief Returns the first supported coordinate. */
  [[nodiscard]] double minimumCoordinate() const noexcept {
    return coordinates_.front();
  }

  /** @brief Returns the last supported coordinate. */
  [[nodiscard]] double maximumCoordinate() const noexcept {
    return coordinates_.back();
  }

  /** @brief Returns the interpolation rule fixed at construction. */
  [[nodiscard]] RateInterpolation interpolation() const noexcept {
    return interpolation_;
  }

  /** @brief Returns the out-of-range rule fixed at construction. */
  [[nodiscard]] RateTableBounds bounds() const noexcept { return bounds_; }

 private:
  /** @brief Rejects malformed axes and coefficients before first evaluation. */
  void validateTable() const {
    if (coordinates_.size() < 2 ||
        coordinates_.size() != coefficients_.size()) {
      throw std::invalid_argument(
          "rate table requires equally sized coordinate and coefficient "
          "arrays with at least two entries");
    }

    for (std::size_t index = 0; index < coordinates_.size(); ++index) {
      if (!std::isfinite(coordinates_[index])) {
        throw std::invalid_argument(
            "rate-table coordinates must contain only finite values");
      }
      if (index > 0 && coordinates_[index] <= coordinates_[index - 1]) {
        throw std::invalid_argument(
            "rate-table coordinates must be strictly increasing");
      }
      if (!std::isfinite(coefficients_[index]) || coefficients_[index] < 0.0) {
        throw std::invalid_argument(
            "rate coefficients must be finite and non-negative");
      }
      if (interpolation_ == RateInterpolation::log_linear &&
          coefficients_[index] == 0.0) {
        throw std::invalid_argument(
            "log-linear rate coefficients must be strictly positive");
      }
    }
  }

  std::vector<double> coordinates_;
  std::vector<double> coefficients_;
  RateInterpolation interpolation_;
  RateTableBounds bounds_;
};

/** @brief Evaluates one coefficient table over a cell coordinate field. */
inline void evaluateRateCoefficient(
    const field::CellField<double>& coordinate,
    const TabulatedRateCoefficient& table,
    field::CellField<double>& rate_coefficient) {
  field::ensureSameMesh(coordinate, rate_coefficient);

  for (const double value : coordinate) {
    table.validateCoordinate(value);
  }
  for (mesh::CellId cell = 0; cell < coordinate.mesh().numCells(); ++cell) {
    rate_coefficient[cell] = table(coordinate[cell]);
  }
}

}  // namespace pemu::physics::reaction
