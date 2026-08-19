#pragma once

#include <pemu/mesh/types.hpp>
#include <variant>

namespace pemu::boundary {

struct Dirichlet {
  double value;
};

//
// Prescribed outward normal diffusive flux
//
//     -k grad(u) · n = value
//
// where k is the diffusion/conductivity coefficient
// of the equation using this boundary condition.
//
struct Neumann {
  double value;
};

using BoundaryCondition = std::variant<Dirichlet, Neumann>;

template <typename T>
concept BoundaryConditionType =
    std::same_as<std::remove_cvref_t<T>, Dirichlet> ||
    std::same_as<std::remove_cvref_t<T>, Neumann>;

}  // namespace pemu::boundary