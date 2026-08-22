#pragma once

#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>
#include <pemu/output/common/types.hpp>

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <span>
#include <string>
#include <vector>

namespace pemu::output::checkpoint {

enum class FieldAssociation : std::uint8_t { cell, face };

/** @brief Selects an immutable cell field for a checkpoint. */
struct CellFieldSource {
  const field::CellField<double>* field{};
  std::string key;
};

/** @brief Selects an immutable face field for a checkpoint. */
struct FaceFieldSource {
  const field::FaceField<double>* field{};
  std::string key;
};

/** @brief Complete request for one restartable checkpoint artifact. */
struct WriteRequest {
  const mesh::IMesh* mesh{};
  std::filesystem::path path;
  OutputStamp stamp;
  std::span<const CellFieldSource> cell_fields;
  std::span<const FaceFieldSource> face_fields;
  bool overwrite{};
};

/** @brief Selects a mutable cell field restored from a checkpoint key. */
struct CellFieldTarget {
  field::CellField<double>* field{};
  std::string key;
};

/** @brief Selects a mutable face field restored from a checkpoint key. */
struct FaceFieldTarget {
  field::FaceField<double>* field{};
  std::string key;
};

/** @brief Request for an all-or-nothing restore into an existing mesh state. */
struct RestoreRequest {
  const mesh::IMesh* mesh{};
  std::span<const CellFieldTarget> cell_fields;
  std::span<const FaceFieldTarget> face_fields;
  bool require_all_fields{true};
};

/** @brief Metadata reported without loading values into simulation fields. */
struct FieldDescriptor {
  FieldAssociation association{};
  std::string key;
  field::FieldMetadata metadata;
  std::size_t value_count{};
};

/** @brief Versioned checkpoint manifest. */
struct Manifest {
  std::uint32_t format_version{};
  OutputStamp stamp;
  std::size_t num_cells{};
  std::size_t num_faces{};
  std::vector<FieldDescriptor> fields;
};

}  // namespace pemu::output::checkpoint
