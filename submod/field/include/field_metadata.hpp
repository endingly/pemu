#pragma once

#include <pemu/unit/quantity_metadata.hpp>

#include <optional>
#include <string>

namespace pemu::field {

struct FieldMetadata {
  std::string name;
  std::optional<pemu::unit::PhysicalQuantityMetadata> physical_quantity;

  [[nodiscard]]
  bool hasPhysicalQuantity() const noexcept {
    return physical_quantity.has_value();
  }
};

}  // namespace pemu::field
