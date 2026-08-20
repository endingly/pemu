#pragma once

#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>
#include <pemu/field/field_metadata.hpp>
#include <pemu/field/types.hpp>

#include <concepts>
#include <cstddef>
#include <span>
#include <vector>

namespace pemu::field {

template <typename Id>
concept DenseFieldId = requires(Id id) {
  { id.value } -> std::convertible_to<std::size_t>;
};

// Owns a homogeneous collection of fields indexed by a domain-specific dense
// ID. Physics may provide aliases with SpeciesId or ReactionId without owning
// another field-container implementation.
template <FieldLike Field, DenseFieldId Id>
class FieldSet {
 public:
  using field_type = Field;
  using value_type = typename Field::value_type;
  using id_type = Id;
  using size_type = std::size_t;

  FieldSet(const mesh::IMesh& mesh, size_type field_count,
           const value_type& initial_value = value_type{},
           FieldMetadata metadata = {})
      : mesh_(&mesh) {
    fields_.reserve(field_count);

    for (size_type i = 0; i < field_count; ++i) {
      fields_.emplace_back(mesh, initial_value, metadata);
    }
  }

  [[nodiscard]]
  size_type size() const noexcept {
    return fields_.size();
  }

  [[nodiscard]]
  bool empty() const noexcept {
    return fields_.empty();
  }

  [[nodiscard]]
  const mesh::IMesh& mesh() const noexcept {
    return *mesh_;
  }

  Field& operator[](Id id) { return fields_.at(index(id)); }

  const Field& operator[](Id id) const { return fields_.at(index(id)); }

  Field& at(Id id) { return fields_.at(index(id)); }

  const Field& at(Id id) const { return fields_.at(index(id)); }

  void fill(const value_type& value) {
    for (auto& field : fields_) {
      field.fill(value);
    }
  }

  [[nodiscard]]
  std::span<Field> span() noexcept {
    return fields_;
  }

  [[nodiscard]]
  std::span<const Field> span() const noexcept {
    return fields_;
  }

  auto begin() noexcept { return fields_.begin(); }
  auto end() noexcept { return fields_.end(); }
  auto begin() const noexcept { return fields_.begin(); }
  auto end() const noexcept { return fields_.end(); }

 private:
  [[nodiscard]]
  static size_type index(Id id) noexcept {
    return static_cast<size_type>(id.value);
  }

  const mesh::IMesh* mesh_;
  std::vector<Field> fields_;
};

template <typename T, DenseFieldId Id>
using CellFieldSet = FieldSet<CellField<T>, Id>;

template <typename T, DenseFieldId Id>
using FaceFieldSet = FieldSet<FaceField<T>, Id>;

}  // namespace pemu::field
