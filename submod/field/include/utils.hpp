#pragma once

#include <stdexcept>
#include <pemu/field/types.hpp>

namespace pemu::field {

template <FieldLike FieldA, FieldLike FieldB>
void ensureSameMesh(const FieldA& a, const FieldB& b) {
  if (&a.mesh() != &b.mesh()) {
    throw std::invalid_argument("fields belong to different meshes");
  }
}

}  // namespace pemu::field