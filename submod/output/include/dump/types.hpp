#pragma once

#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>
#include <pemu/output/common/types.hpp>

#include <filesystem>
#include <functional>
#include <span>
#include <string>

namespace pemu::output::dump {

/** @brief Describes one visualization-dump series. */
struct SeriesRequest {
  const mesh::IMesh* mesh{};
  std::filesystem::path path;
  bool overwrite{};
};

/** @brief Metadata to serialize for a selected visualization field. */
struct FieldMetadata {
  std::string meta_name;
  unit::QuantityKind quantity_kind{unit::QuantityKind::dimensionless};
  units::precise_unit unit{units::precise::one};
};

/** @brief Selects a cell field and optionally overrides its dump name. */
struct CellFieldSelection {
  const field::CellField<double>* field{};
  std::string name;
  FieldMetadata meta_data;
};

/** @brief Selects a face field and its optional visualization projection. */
struct FaceFieldSelection {
  const field::FaceField<double>* field{};
  std::string name;
  bool include_cell_centered_visualization{};
  std::string cell_centered_name;
  FieldMetadata meta_data;
};

/** @brief Describes one visualization snapshot. */
struct Request {
  const mesh::IMesh* mesh{};
  std::filesystem::path path;
  OutputStamp stamp;
  std::span<const std::reference_wrapper<const field::CellField<double>>>
      cell_fields;
  std::span<const CellFieldSelection> cell_field_selections;
  std::span<const FaceFieldSelection> face_fields;
  bool overwrite{};
};

}  // namespace pemu::output::dump
