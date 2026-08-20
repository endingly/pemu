#pragma once

#include <pemu/field/field_metadata.hpp>

#include <mp-units/framework/reference.h>
#include <mp-units/framework/unit.h>

#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace pemu::field {

template <mp_units::Reference R>
[[nodiscard]] FieldMetadata makeFieldMetadata(std::string name, R reference) {
  const auto quantity_spec = mp_units::get_quantity_spec(reference);
  const auto unit = mp_units::get_unit(reference);

  return {
      .name = std::move(name),
      .physical_quantity = pemu::unit::PhysicalQuantityMetadata(
          std::string(mp_units::unit_symbol(unit)),
          typeid(std::remove_cvref_t<decltype(quantity_spec)>),
          typeid(std::remove_cvref_t<decltype(unit)>)),
  };
}

template <mp_units::Reference R>
void ensureReference(const FieldMetadata& metadata, R expected_reference) {
  if (!metadata.physical_quantity.has_value()) {
    throw std::invalid_argument("field has no physical quantity metadata");
  }

  if (!metadata.physical_quantity->represents(
          mp_units::get_quantity_spec(expected_reference),
          mp_units::get_unit(expected_reference))) {
    throw std::invalid_argument(
        "field physical quantity metadata does not match requested reference");
  }
}

}  // namespace pemu::field
