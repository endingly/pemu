#pragma once

#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>

#include <cstdint>
#include <filesystem>
#include <functional>
#include <span>
#include <string>

namespace pemu::output {

struct OutputStamp {
  std::uint64_t step{};
  double time{};
};

struct OutputRecord {
  std::filesystem::path path;
  OutputStamp stamp;
};

/** @brief Describes one output series before any temporal snapshots arrive. */
struct FieldSeriesRequest {
  const mesh::IMesh* mesh{};
  std::filesystem::path path;
  bool overwrite{};
};

/**
 * @brief Selects a cell field for output and optionally overrides its output
 * name without copying its data.
 */
struct CellFieldSelection {
  const field::CellField<double>* field{};
  std::string name;
};

struct FaceFieldSelection {
  const field::FaceField<double>* field{};
  std::string name;
  bool include_cell_centered_visualization{};
  std::string cell_centered_name;
};

struct FieldDumpRequest {
  const mesh::IMesh* mesh{};
  std::filesystem::path path;
  OutputStamp stamp;
  std::span<const std::reference_wrapper<const field::CellField<double>>>
      cell_fields;
  std::span<const CellFieldSelection> cell_field_selections;
  std::span<const FaceFieldSelection> face_fields;
  bool overwrite{};
};

}  // namespace pemu::output
