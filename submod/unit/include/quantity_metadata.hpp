#pragma once

#include <string>
#include <type_traits>
#include <typeindex>
#include <utility>

namespace pemu::unit {

// Runtime, type-erased metadata for an mp-units quantity specification and
// unit. The mp-units types never become part of a field's scalar type.
class PhysicalQuantityMetadata {
 public:
  PhysicalQuantityMetadata(std::string unit_symbol,
                           std::type_index quantity_spec_type,
                           std::type_index unit_type)
      : unit_symbol_(std::move(unit_symbol)),
        quantity_spec_type_(quantity_spec_type),
        unit_type_(unit_type) {}

  [[nodiscard]]
  const std::string& unitSymbol() const noexcept {
    return unit_symbol_;
  }

  template <typename QS, typename U>
  [[nodiscard]] bool represents(QS, U) const noexcept {
    return quantity_spec_type_ ==
               std::type_index(typeid(std::remove_cvref_t<QS>)) &&
           unit_type_ == std::type_index(typeid(std::remove_cvref_t<U>));
  }

  friend bool operator==(const PhysicalQuantityMetadata&,
                         const PhysicalQuantityMetadata&) = default;

 private:
  std::string unit_symbol_;
  std::type_index quantity_spec_type_;
  std::type_index unit_type_;
};

}  // namespace pemu::unit
