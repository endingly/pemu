#pragma once

#include <pemu/field/mp_units_metadata.hpp>
#include <pemu/field/types.hpp>

#include <mp-units/framework/quantity.h>

#include <concepts>
#include <type_traits>

namespace pemu::field {

// Unit-aware boundary adapters. They validate metadata and then immediately
// convert to/from the raw representation used by numerical kernels.
template <FieldLike Field, typename Index, mp_units::Reference R>
  requires std::is_arithmetic_v<typename Field::value_type>
[[nodiscard]] auto quantityAt(const Field& field, Index index, R reference) {
  ensureReference(field.metadata(), reference);
  return static_cast<typename Field::value_type>(field[index]) * reference;
}

template <FieldLike Field, typename Index, mp_units::Quantity Q,
          mp_units::Reference R>
  requires(std::is_arithmetic_v<typename Field::value_type> &&
           mp_units::implicitly_convertible(Q::quantity_spec,
                                            mp_units::get_quantity_spec(R{})) &&
           requires(Q quantity, R reference) {
             quantity.numerical_value_in(mp_units::get_unit(reference));
           })
void setQuantity(Field& field, Index index, Q quantity, R reference) {
  ensureReference(field.metadata(), reference);
  field[index] = static_cast<typename Field::value_type>(
      quantity.numerical_value_in(mp_units::get_unit(reference)));
}

template <FieldLike Field, mp_units::Quantity Q, mp_units::Reference R>
  requires(std::is_arithmetic_v<typename Field::value_type> &&
           mp_units::implicitly_convertible(Q::quantity_spec,
                                            mp_units::get_quantity_spec(R{})) &&
           requires(Q quantity, R reference) {
             quantity.numerical_value_in(mp_units::get_unit(reference));
           })
void fillQuantity(Field& field, Q quantity, R reference) {
  ensureReference(field.metadata(), reference);
  field.fill(static_cast<typename Field::value_type>(
      quantity.numerical_value_in(mp_units::get_unit(reference))));
}

}  // namespace pemu::field
