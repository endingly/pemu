#pragma once

#include <pemu/mesh/types.hpp>

namespace pemu::mesh {

[[nodiscard]]
constexpr double dot(const Vec3& a, const Vec3& b) noexcept {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

[[nodiscard]]
constexpr double normSquared(const Vec3& v) noexcept {
  return dot(v, v);
}

}  // namespace pemu::mesh