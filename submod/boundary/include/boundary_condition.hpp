#pragma once

#include <pemu/mesh/types.hpp>
#include <variant>

namespace pemu::boundary {

struct Dirichlet {
  double value;
};

//
// Specifies outward flux:
//
//   -epsilon grad(phi) · n = value
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