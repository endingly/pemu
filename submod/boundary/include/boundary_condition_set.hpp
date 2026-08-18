#pragma once

#include <pemu/boundary/boundary_condition.hpp>
#include <pemu/mesh/types.hpp>

#include <concepts>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <utility>

namespace pemu::boundary {

class BoundaryConditionSet {
 public:
  template <BoundaryConditionType T>
  void set(mesh::BoundaryId id, T&& condition) {
    conditions_.insert_or_assign(id,
                                 BoundaryCondition{std::forward<T>(condition)});
  }

  void setDirichlet(mesh::BoundaryId id, double value) {
    set(id, Dirichlet{value});
  }

  void setNeumann(mesh::BoundaryId id, double value) {
    set(id, Neumann{value});
  }

  bool contains(mesh::BoundaryId id) const noexcept {
    return conditions_.contains(id);
  }

  const BoundaryCondition& at(mesh::BoundaryId id) const {
    const auto it = conditions_.find(id);

    if (it == conditions_.end()) {

      throw std::out_of_range(
          "boundary condition "
          "is not defined");
    }

    return it->second;
  }

  std::size_t size() const noexcept { return conditions_.size(); }

  void clear() noexcept { conditions_.clear(); }

 private:
  std::unordered_map<mesh::BoundaryId, BoundaryCondition> conditions_;
};

}  // namespace pemu::boundary