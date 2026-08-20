#include <pemu/equation/poisson_solver.hpp>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace pemu::equation {

PoissonSolver::PoissonSolver(
    discretization::PoissonFvm discretization,
    std::unique_ptr<linalg::ISolver> linear_solver,
    std::optional<PureNeumannOptions> pure_neumann_options)
    : discretization_(std::move(discretization)),
      linear_solver_(std::move(linear_solver)),
      pure_neumann_options_(std::move(pure_neumann_options)),
      pure_neumann_(!discretization_.hasDirichletBoundary()) {
  if (!linear_solver_) {
    throw std::invalid_argument("PoissonSolver requires a linear solver");
  }
  if (pure_neumann_ && !pure_neumann_options_) {
    throw std::invalid_argument(
        "pure Neumann Poisson problem requires a gauge condition");
  }
  if (!pure_neumann_ || !pure_neumann_options_) {
    return;
  }
  if (!std::isfinite(pure_neumann_options_->absolute_compatibility_tolerance) ||
      !std::isfinite(pure_neumann_options_->relative_compatibility_tolerance) ||
      pure_neumann_options_->absolute_compatibility_tolerance < 0.0 ||
      pure_neumann_options_->relative_compatibility_tolerance < 0.0) {
    throw std::invalid_argument(
        "pure Neumann compatibility tolerances must be finite and "
        "non-negative");
  }
  if (const auto* pin =
          std::get_if<PinCellGauge>(&pure_neumann_options_->gauge);
      pin != nullptr) {
    if (pin->cell >= discretization_.mesh().numCells()) {
      throw std::invalid_argument("pure Neumann pin cell is out of range");
    }
    if (!std::isfinite(pin->value)) {
      throw std::invalid_argument("pure Neumann pin value must be finite");
    }
  }
}

linalg::SolverResult PoissonSolver::solve(field::CellField<double>& solution) {
  if (&solution.mesh() != &discretization_.mesh()) {
    throw std::invalid_argument(
        "Poisson solution field "
        "belongs to another mesh");
  }
  if (!initialized_) {
    const auto status = initialize();
    if (status != linalg::SolverStatus::Success) {
      return {.status = status};
    }
  }
  // --------------------------------------------------------
  // Only RHS is rebuilt.
  // --------------------------------------------------------
  discretization_.assembleRhs(b_);

  if (pure_neumann_) {
    if (!hasCompatibleRhs()) {
      return {.status = linalg::SolverStatus::IncompatibleRhs};
    }
    applyGaugeToRhs();
  }
  // --------------------------------------------------------
  // Solve directly into CellField storage unless the zero-mean gauge adds a
  // Lagrange multiplier to the system.
  // --------------------------------------------------------
  if (usesZeroMeanGauge()) {
    augmented_solution_.resize(b_.size());
    const auto result = linear_solver_->solve(b_, augmented_solution_);
    if (result.success()) {
      for (mesh::CellId cell = 0; cell < solution.size(); ++cell) {
        solution[cell] = augmented_solution_[static_cast<linalg::Index>(cell)];
      }
    }
    return result;
  }

  Eigen::Map<linalg::Vector> x(solution.data(),
                               static_cast<linalg::Index>(solution.size()));
  return linear_solver_->solve(b_, x);
}

void PoissonSolver::reset() {
  linear_solver_->reset();
  A_.resize(0, 0);
  b_.resize(0);
  pin_column_.resize(0);
  augmented_solution_.resize(0);
  initialized_ = false;
}

linalg::SolverStatus PoissonSolver::initialize() {
  discretization_.assembleMatrix(A_);

  if (pure_neumann_) {
    applyGaugeToMatrix();
  }

  const auto analyze_status = linear_solver_->analyzePattern(A_);

  if (analyze_status != linalg::SolverStatus::Success) {

    initialized_ = false;

    return analyze_status;
  }

  const auto factor_status = linear_solver_->factorize(A_);

  if (factor_status != linalg::SolverStatus::Success) {

    initialized_ = false;

    return factor_status;
  }

  initialized_ = true;

  return linalg::SolverStatus::Success;
}

void PoissonSolver::applyGaugeToMatrix() {
  std::visit(
      [this](const auto& gauge) {
        using Gauge = std::remove_cvref_t<decltype(gauge)>;

        if constexpr (std::same_as<Gauge, PinCellGauge>) {
          using Triplet = Eigen::Triplet<linalg::Scalar, linalg::Index>;

          const auto pin = static_cast<linalg::Index>(gauge.cell);
          pin_column_ = linalg::Vector::Zero(A_.rows());

          for (linalg::SparseMatrix::InnerIterator entry(A_, pin); entry;
               ++entry) {
            pin_column_[entry.row()] = entry.value();
          }

          std::vector<Triplet> entries;
          entries.reserve(static_cast<std::size_t>(A_.nonZeros()) + 1);
          for (linalg::Index column = 0; column < A_.outerSize(); ++column) {
            for (linalg::SparseMatrix::InnerIterator entry(A_, column); entry;
                 ++entry) {
              if (entry.row() != pin && entry.col() != pin) {
                entries.emplace_back(entry.row(), entry.col(), entry.value());
              }
            }
          }
          entries.emplace_back(pin, pin, 1.0);

          linalg::SparseMatrix constrained(A_.rows(), A_.cols());
          constrained.setFromTriplets(entries.begin(), entries.end());
          constrained.makeCompressed();
          A_.swap(constrained);
        } else if constexpr (std::same_as<Gauge, ZeroMeanGauge>) {
          using Triplet = Eigen::Triplet<linalg::Scalar, linalg::Index>;

          const auto cells =
              static_cast<linalg::Index>(discretization_.mesh().numCells());
          double total_volume = 0.0;
          for (mesh::CellId cell = 0; cell < discretization_.mesh().numCells();
               ++cell) {
            total_volume += discretization_.mesh().cellVolume(cell);
          }
          if (total_volume <= 0.0) {
            throw std::runtime_error("mesh has non-positive total volume");
          }

          std::vector<Triplet> entries;
          entries.reserve(static_cast<std::size_t>(A_.nonZeros()) +
                          2 * discretization_.mesh().numCells());
          for (linalg::Index column = 0; column < A_.outerSize(); ++column) {
            for (linalg::SparseMatrix::InnerIterator entry(A_, column); entry;
                 ++entry) {
              entries.emplace_back(entry.row(), entry.col(), entry.value());
            }
          }
          for (mesh::CellId cell = 0; cell < discretization_.mesh().numCells();
               ++cell) {
            const auto index = static_cast<linalg::Index>(cell);
            const double weight =
                discretization_.mesh().cellVolume(cell) / total_volume;
            entries.emplace_back(index, cells, weight);
            entries.emplace_back(cells, index, weight);
          }

          linalg::SparseMatrix augmented(cells + 1, cells + 1);
          augmented.setFromTriplets(entries.begin(), entries.end());
          augmented.makeCompressed();
          A_.swap(augmented);
        }
      },
      pure_neumann_options_->gauge);
}

void PoissonSolver::applyGaugeToRhs() {
  std::visit(
      [this](const auto& gauge) {
        using Gauge = std::remove_cvref_t<decltype(gauge)>;

        if constexpr (std::same_as<Gauge, PinCellGauge>) {
          const auto pin = static_cast<linalg::Index>(gauge.cell);
          b_ -= pin_column_ * gauge.value;
          b_[pin] = gauge.value;
        } else if constexpr (std::same_as<Gauge, ZeroMeanGauge>) {
          const auto cells = b_.size();
          b_.conservativeResize(cells + 1);
          b_[cells] = 0.0;
        }
      },
      pure_neumann_options_->gauge);
}

bool PoissonSolver::hasCompatibleRhs() const {
  const double residual = std::abs(b_.sum());
  const double scale = b_.lpNorm<1>();
  return residual <=
         pure_neumann_options_->absolute_compatibility_tolerance +
             pure_neumann_options_->relative_compatibility_tolerance * scale;
}

bool PoissonSolver::usesZeroMeanGauge() const noexcept {
  return pure_neumann_ && pure_neumann_options_ &&
         std::holds_alternative<ZeroMeanGauge>(pure_neumann_options_->gauge);
}

}  // namespace pemu::equation
