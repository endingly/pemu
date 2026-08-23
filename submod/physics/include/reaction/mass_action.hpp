#pragma once

#include <pemu/field/cell_field.hpp>
#include <pemu/field/utils.hpp>

#include <llnl-units/units.hpp>

#include <array>
#include <cmath>
#include <limits>
#include <span>
#include <stdexcept>

namespace pemu::physics::reaction {

/**
 * @brief Read-only density factor n^alpha used by a mass-action rate law.
 *
 * A factor can reference a spatial field or hold a uniform background value.
 * This permits prescribed bath gases without introducing a separate kernel.
 */
class MassActionTerm {
 public:
  /** @brief Creates a spatially varying density factor. */
  explicit MassActionTerm(const field::CellField<double>& number_density,
                          double order = 1.0)
      : number_density_(&number_density), order_(order) {
    validateOrder();
  }

  /** @brief Creates a spatially uniform density factor. */
  explicit MassActionTerm(double uniform_number_density, double order = 1.0)
      : uniform_number_density_(uniform_number_density), order_(order) {
    validateOrder();
    if (!std::isfinite(uniform_number_density_) ||
        uniform_number_density_ < 0.0) {
      throw std::invalid_argument(
          "uniform mass-action density must be finite and non-negative");
    }
  }

  /** @brief Returns the kinetic exponent applied to this density. */
  [[nodiscard]] double order() const noexcept { return order_; }

  /** @brief Reports whether this factor references a cell field. */
  [[nodiscard]] bool isField() const noexcept {
    return number_density_ != nullptr;
  }

  /** @brief Returns the referenced field, or null for a uniform factor. */
  [[nodiscard]] const field::CellField<double>* field() const noexcept {
    return number_density_;
  }

  /** @brief Returns the raw density value for one output cell. */
  [[nodiscard]] double value(mesh::CellId cell) const noexcept {
    return isField() ? (*number_density_)[cell] : uniform_number_density_;
  }

 private:
  /** @brief Rejects non-positive or non-finite kinetic orders. */
  void validateOrder() const {
    if (!std::isfinite(order_) || order_ <= 0.0) {
      throw std::invalid_argument(
          "mass-action order must be finite and positive");
    }
  }

  const field::CellField<double>* number_density_{};
  double uniform_number_density_{};
  double order_{1.0};
};

/** @brief Read-only spatial or uniform mass-action rate coefficient. */
class RateCoefficientView {
 public:
  /** @brief References a spatially varying rate-coefficient field. */
  explicit RateCoefficientView(const field::CellField<double>& rate_coefficient)
      : rate_coefficient_(&rate_coefficient) {}

  /** @brief Stores a spatially uniform rate coefficient. */
  explicit RateCoefficientView(double uniform_rate_coefficient)
      : uniform_rate_coefficient_(uniform_rate_coefficient) {
    if (!std::isfinite(uniform_rate_coefficient_) ||
        uniform_rate_coefficient_ < 0.0) {
      throw std::invalid_argument(
          "uniform rate coefficient must be finite and non-negative");
    }
  }

  /** @brief Reports whether this coefficient references a cell field. */
  [[nodiscard]] bool isField() const noexcept {
    return rate_coefficient_ != nullptr;
  }

  /** @brief Returns the referenced field, or null for a uniform coefficient. */
  [[nodiscard]] const field::CellField<double>* field() const noexcept {
    return rate_coefficient_;
  }

  /** @brief Returns the raw coefficient value for one output cell. */
  [[nodiscard]] double value(mesh::CellId cell) const noexcept {
    return isField() ? (*rate_coefficient_)[cell] : uniform_rate_coefficient_;
  }

 private:
  const field::CellField<double>* rate_coefficient_{};
  double uniform_rate_coefficient_{};
};

/** @brief Returns the sum of all kinetic exponents in a mass-action law. */
[[nodiscard]] inline double totalReactionOrder(
    std::span<const MassActionTerm> terms) {
  if (terms.empty()) {
    throw std::invalid_argument(
        "mass-action reaction requires at least one density factor");
  }
  double result = 0.0;
  for (const auto& term : terms) {
    result += term.order();
    if (!std::isfinite(result)) {
      throw std::overflow_error("total reaction order overflowed");
    }
  }
  return result;
}

/**
 * @brief Returns the canonical cm-based coefficient unit for molecularity m.
 *
 * With every density in 1/cm^3 and the reaction rate in 1/(cm^3*s), an
 * elementary m-body coefficient has unit cm^(3*(m-1))/s.
 */
[[nodiscard]] inline units::precise_unit canonicalCentimetreRateCoefficientUnit(
    unsigned int molecularity) {
  if (molecularity == 0) {
    throw std::invalid_argument("reaction molecularity must be positive");
  }
  constexpr auto maximum_exponent =
      static_cast<unsigned int>(std::numeric_limits<int>::max() / 3);
  if (molecularity - 1 > maximum_exponent) {
    throw std::overflow_error("reaction molecularity is too large");
  }
  const int length_exponent = 3 * static_cast<int>(molecularity - 1);
  return units::precise::cm.pow(length_exponent) / units::precise::s;
}

namespace detail {

/** @brief Multiplies non-negative factors while detecting overflow. */
[[nodiscard]] inline double checkedNonNegativeProduct(double lhs, double rhs) {
  if (lhs == 0.0 || rhs == 0.0) {
    return 0.0;
  }
  if (lhs > std::numeric_limits<double>::max() / rhs) {
    throw std::overflow_error("mass-action reaction rate overflowed");
  }
  return lhs * rhs;
}

/** @brief Evaluates one already mesh-validated mass-action cell. */
[[nodiscard]] inline double evaluateMassActionCell(
    std::span<const MassActionTerm> terms,
    const RateCoefficientView& rate_coefficient, mesh::CellId cell) {
  double result = rate_coefficient.value(cell);
  if (!std::isfinite(result) || result < 0.0) {
    throw std::invalid_argument(
        "mass-action rate coefficient must be finite and non-negative");
  }
  if (result == 0.0) {
    return 0.0;
  }

  for (const auto& term : terms) {
    const double density = term.value(cell);
    if (!std::isfinite(density) || density < 0.0) {
      throw std::invalid_argument(
          "mass-action densities must be finite and non-negative");
    }
    if (density == 0.0) {
      return 0.0;
    }
    const double factor = std::pow(density, term.order());
    if (!std::isfinite(factor)) {
      throw std::overflow_error("mass-action density factor overflowed");
    }
    result = checkedNonNegativeProduct(result, factor);
  }
  return result;
}

/** @brief Writes rates after every input field and value has been validated. */
inline void writeValidatedMassActionReactionRate(
    std::span<const MassActionTerm> terms,
    const RateCoefficientView& rate_coefficient,
    field::CellField<double>& reaction_rate) {
  for (mesh::CellId cell = 0; cell < reaction_rate.mesh().numCells(); ++cell) {
    reaction_rate[cell] = evaluateMassActionCell(terms, rate_coefficient, cell);
  }
}

}  // namespace detail

/**
 * @brief Validates all inputs to one arbitrary-order mass-action evaluation.
 *
 * This operation does not modify the destination field. It is exposed so a
 * multi-reaction assembler can validate its complete batch before any write.
 */
inline void validateMassActionReactionRate(
    std::span<const MassActionTerm> terms,
    const RateCoefficientView& rate_coefficient,
    const field::CellField<double>& reaction_rate) {
  static_cast<void>(totalReactionOrder(terms));
  if (rate_coefficient.isField()) {
    field::ensureSameMesh(*rate_coefficient.field(), reaction_rate);
  }
  for (const auto& term : terms) {
    if (term.isField()) {
      field::ensureSameMesh(*term.field(), reaction_rate);
    }
  }

  for (mesh::CellId cell = 0; cell < reaction_rate.mesh().numCells(); ++cell) {
    static_cast<void>(
        detail::evaluateMassActionCell(terms, rate_coefficient, cell));
  }
}

/**
 * @brief Computes R = k product_j n_j^alpha_j for any number of factors.
 *
 * Inputs must use one consistent unit system. The complete candidate rate is
 * validated before the destination is overwritten, preserving atomic failure.
 */
inline void massActionReactionRate(std::span<const MassActionTerm> terms,
                                   const RateCoefficientView& rate_coefficient,
                                   field::CellField<double>& reaction_rate) {
  validateMassActionReactionRate(terms, rate_coefficient, reaction_rate);
  detail::writeValidatedMassActionReactionRate(terms, rate_coefficient,
                                               reaction_rate);
}

/** @brief Convenience wrapper for a spatial binary reaction. */
inline void binaryReactionRate(
    const field::CellField<double>& first_number_density,
    const field::CellField<double>& second_number_density,
    const field::CellField<double>& rate_coefficient,
    field::CellField<double>& reaction_rate) {
  const std::array terms{MassActionTerm{first_number_density},
                         MassActionTerm{second_number_density}};
  massActionReactionRate(terms, RateCoefficientView{rate_coefficient},
                         reaction_rate);
}

/** @brief Convenience wrapper for a binary reaction with uniform target and k. */
inline void binaryReactionRate(
    const field::CellField<double>& first_number_density,
    double second_number_density, double rate_coefficient,
    field::CellField<double>& reaction_rate) {
  const std::array terms{MassActionTerm{first_number_density},
                         MassActionTerm{second_number_density}};
  massActionReactionRate(terms, RateCoefficientView{rate_coefficient},
                         reaction_rate);
}

}  // namespace pemu::physics::reaction
