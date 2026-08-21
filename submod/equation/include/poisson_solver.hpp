#pragma once

#include <pemu/discretization/poisson_fvm.hpp>
#include <pemu/field/cell_field.hpp>
#include <pemu/linalg/i_solver.hpp>

#include <memory>
#include <optional>
#include <variant>

namespace pemu::equation {

/**
 * @brief Fixes one cell value to remove the constant null space of a pure
 * Neumann Poisson problem.
 */
struct PinCellGauge {
  mesh::CellId cell{0};
  double value{0.0};
};

/**
 * @brief Constrains the volume-weighted mean solution to zero for a pure
 * Neumann Poisson problem.
 */
struct ZeroMeanGauge {};

/** @brief Gauge strategies supported for pure Neumann Poisson systems. */
using PureNeumannGauge = std::variant<PinCellGauge, ZeroMeanGauge>;

/**
 * @brief Configuration for solving a pure Neumann Poisson problem.
 */
struct PureNeumannOptions {
  PureNeumannGauge gauge;
  double absolute_compatibility_tolerance{1e-12};
  double relative_compatibility_tolerance{1e-10};
};

class PoissonSolver {
 public:
  /**
   * @brief Creates a reusable Poisson solver.
   *
   * @param discretization Finite-volume Poisson discretization.
   * @param linear_solver Backend used to factorize and solve the system.
   * @param pure_neumann_options Required when every boundary is Neumann and
   * ignored when at least one Dirichlet boundary is present.
   * @throws std::invalid_argument for a null backend, a pure Neumann problem
   * without gauge options, an invalid pin configuration, or invalid
   * compatibility tolerances.
   */
  PoissonSolver(
      discretization::PoissonFvm discretization,
      std::unique_ptr<linalg::ISolver> linear_solver,
      std::optional<PureNeumannOptions> pure_neumann_options = std::nullopt);

  PoissonSolver(const PoissonSolver&) = delete;

  PoissonSolver& operator=(const PoissonSolver&) = delete;

  PoissonSolver(PoissonSolver&&) noexcept = default;

  PoissonSolver& operator=(PoissonSolver&&) noexcept = default;

  linalg::SolverResult solve(field::CellField<double>& solution);

  void reset();

  linalg::SolverResult initialize();

  bool is_initialized() const noexcept { return initialized_; }

 private:
  /** @brief Applies the configured gauge to the assembled matrix. */
  void applyGaugeToMatrix();

  /** @brief Applies the configured gauge to an already assembled RHS. */
  void applyGaugeToRhs();

  /**
   * @brief Checks the discrete pure Neumann compatibility condition.
   * @return `true` when the RHS is orthogonal to the constant null space.
   */
  [[nodiscard]]
  bool hasCompatibleRhs() const;

  /** @brief Returns whether the configured gauge is the zero-mean strategy. */
  [[nodiscard]]
  bool usesZeroMeanGauge() const noexcept;

  discretization::PoissonFvm discretization_;
  std::unique_ptr<linalg::ISolver> linear_solver_;
  std::optional<PureNeumannOptions> pure_neumann_options_;
  linalg::SparseMatrix A_;
  linalg::Vector b_;
  linalg::Vector pin_column_;
  linalg::Vector augmented_solution_;
  bool pure_neumann_{false};
  bool initialized_{false};
};

}  // namespace pemu::equation
