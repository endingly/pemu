#pragma once

#include <pemu/discretization/poisson_fvm.hpp>
#include <pemu/field/cell_field.hpp>
#include <pemu/linalg/i_solver.hpp>

#include <memory>

namespace pemu::equation {

class PoissonSolver {
 public:
  PoissonSolver(discretization::PoissonFvm discretization,
                std::unique_ptr<linalg::ISolver> linear_solver);

  PoissonSolver(const PoissonSolver&) = delete;

  PoissonSolver& operator=(const PoissonSolver&) = delete;

  PoissonSolver(PoissonSolver&&) noexcept = default;

  PoissonSolver& operator=(PoissonSolver&&) noexcept = default;

  linalg::SolverResult solve(field::CellField<double>& solution);

  void reset();

 private:
  discretization::PoissonFvm discretization_;

  std::unique_ptr<linalg::ISolver> linear_solver_;
};

}  // namespace pemu::equation