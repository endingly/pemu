#pragma once

#include <pemu/mesh/i_mesh.hpp>

#include <concepts>
#include <cstddef>

namespace pemu::field {

enum class Location { Cell, Face };

template <typename T>
concept FieldLike = requires(T field, const T const_field) {
  typename T::value_type;
  typename T::size_type;

  { const_field.size() } -> std::convertible_to<std::size_t>;

  { const_field.mesh() } -> std::same_as<const mesh::IMesh&>;

  { field.data() } -> std::same_as<typename T::value_type*>;

  { const_field.data() } -> std::same_as<const typename T::value_type*>;
};

template <typename T>
concept CellFieldLike = FieldLike<T> && (T::location == Location::Cell);

template <typename T>
concept FaceFieldLike = FieldLike<T> && (T::location == Location::Face);

};  // namespace pemu::field
