#pragma once

#include <pemu/field/types.hpp>

#include <algorithm>
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

namespace pemu::field {

template <typename T>
class CellField {
 public:
  using value_type = T;
  using size_type = std::size_t;
  static constexpr Location location = Location::Cell;

 public:
  explicit CellField(const mesh::IMesh& mesh)
      : mesh_(&mesh), data_(mesh.numCells()) {}

  CellField(const mesh::IMesh& mesh, const T& initial_value)
      : mesh_(&mesh), data_(mesh.numCells(), initial_value) {}

  [[nodiscard]]
  size_type size() const noexcept {
    return data_.size();
  }

  [[nodiscard]]
  bool empty() const noexcept {
    return data_.empty();
  }

  [[nodiscard]]
  const mesh::IMesh& mesh() const noexcept {
    return *mesh_;
  }

  T& operator[](const mesh::CellId cell) noexcept {
    return data_[static_cast<size_type>(cell)];
  }

  const T& operator[](const mesh::CellId cell) const noexcept {
    return data_[static_cast<size_type>(cell)];
  }

  T& at(const mesh::CellId cell) {
    return data_.at(static_cast<size_type>(cell));
  }

  const T& at(const mesh::CellId cell) const {
    return data_.at(static_cast<size_type>(cell));
  }

  [[nodiscard]]
  T* data() noexcept {
    return data_.data();
  }

  [[nodiscard]]
  const T* data() const noexcept {
    return data_.data();
  }

  [[nodiscard]]
  std::span<T> span() noexcept {
    return data_;
  }

  [[nodiscard]]
  std::span<const T> span() const noexcept {
    return data_;
  }

  auto begin() noexcept { return data_.begin(); }

  auto end() noexcept { return data_.end(); }

  auto begin() const noexcept { return data_.begin(); }

  auto end() const noexcept { return data_.end(); }

  void fill(const T& value) { std::fill(data_.begin(), data_.end(), value); }

 private:
  const mesh::IMesh* mesh_;
  std::vector<T> data_;
};

}  // namespace pemu::field