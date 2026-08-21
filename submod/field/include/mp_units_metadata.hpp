#pragma once

#include <pemu/field/field_metadata.hpp>
#include <pemu/unit/mp_units_bridge.hpp>

#include <mp-units/framework/reference.h>

#include <stdexcept>
#include <string>
#include <utility>

namespace pemu::field {

template <mp_units::Reference R>
[[nodiscard]] FieldMetadata makeFieldMetadata(std::string name, R reference) {
  return {
      .name = std::move(name),
      .physical_quantity = pemu::unit::bridgeReference(reference),
  };
}

template <mp_units::Reference R>
void ensureReference(const FieldMetadata& metadata, R expected_reference) {
  if (!metadata.physical_quantity.has_value()) {
    throw std::invalid_argument("field has no physical quantity metadata");
  }

  if (*metadata.physical_quantity !=
      pemu::unit::bridgeReference(expected_reference)) {
    throw std::invalid_argument(
        "field physical quantity metadata does not match requested reference");
  }
}

}  // namespace pemu::field
