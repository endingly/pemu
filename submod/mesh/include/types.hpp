#pragma once

#include <cstdint>
#include <limits>

namespace pemu::mesh {

using CellId = std::uint32_t;
using FaceId = std::uint32_t;
using VertexId = std::uint32_t;
using BoundaryId = std::uint32_t;

inline constexpr CellId invalid_cell = std::numeric_limits<CellId>::max();

inline constexpr FaceId invalid_face = std::numeric_limits<FaceId>::max();

inline constexpr BoundaryId invalid_boundary =
    std::numeric_limits<BoundaryId>::max();

struct Vec3 {
  double x{};
  double y{};
  double z{};

  constexpr Vec3 operator+(const Vec3& rhs) const noexcept {
    return {x + rhs.x, y + rhs.y, z + rhs.z};
  }

  constexpr Vec3 operator-(const Vec3& rhs) const noexcept {
    return {x - rhs.x, y - rhs.y, z - rhs.z};
  }

  constexpr Vec3 operator*(double s) const noexcept {
    return {x * s, y * s, z * s};
  }

  constexpr Vec3 operator/(double s) const noexcept {
    return {x / s, y / s, z / s};
  }
};

}  // namespace pemu::mesh