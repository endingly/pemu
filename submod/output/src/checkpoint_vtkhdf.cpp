#include <pemu/output/checkpoint/vtkhdf.hpp>

#include <pemu/output/common/mesh_adapter.hpp>

#include <llnl-units/units.hpp>

#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkErrorCode.h>
#include <vtkFieldData.h>
#include <vtkHDFReader.h>
#include <vtkHDFWriter.h>
#include <vtkIdList.h>
#include <vtkIdTypeArray.h>
#include <vtkNew.h>
#include <vtkStringArray.h>
#include <vtkTypeUInt32Array.h>
#include <vtkTypeUInt64Array.h>
#include <vtkUnstructuredGrid.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_set>
#include <utility>

namespace pemu::output::checkpoint {

namespace {

constexpr std::string_view kFormatName = "pemu.vtkhdf.checkpoint";
constexpr std::string_view kFormatArray = "pemu_checkpoint_format";
constexpr std::string_view kVersionArray = "pemu_checkpoint_version";
constexpr std::string_view kStepArray = "pemu_checkpoint_step";
constexpr std::string_view kTimeArray = "pemu_checkpoint_time";
constexpr std::string_view kNumFacesArray = "pemu_checkpoint_num_faces";
constexpr std::string_view kAssociationArray =
    "pemu_checkpoint_field_association";
constexpr std::string_view kKeyArray = "pemu_checkpoint_field_key";
constexpr std::string_view kMetadataNameArray =
    "pemu_checkpoint_field_metadata_name";
constexpr std::string_view kQuantityKindArray =
    "pemu_checkpoint_field_quantity_kind";
constexpr std::string_view kUnitArray = "pemu_checkpoint_field_unit";
constexpr std::string_view kValueCountArray =
    "pemu_checkpoint_field_value_count";
constexpr std::string_view kFaceOwnerArray = "pemu_checkpoint_face_owner";
constexpr std::string_view kFaceNeighborArray = "pemu_checkpoint_face_neighbor";
constexpr std::string_view kFaceCenterArray = "pemu_checkpoint_face_center";
constexpr std::string_view kFaceAreaArray = "pemu_checkpoint_face_area";
constexpr std::string_view kFaceNormalArray = "pemu_checkpoint_face_normal";
constexpr std::string_view kFaceBoundaryArray =
    "pemu_checkpoint_face_boundary_id";

constexpr std::array<std::string_view, 17> kReservedNames{
    kFormatArray,       kVersionArray,     kStepArray,       kTimeArray,
    kNumFacesArray,     kAssociationArray, kKeyArray,        kMetadataNameArray,
    kQuantityKindArray, kUnitArray,        kValueCountArray, kFaceOwnerArray,
    kFaceNeighborArray, kFaceCenterArray,  kFaceAreaArray,   kFaceNormalArray,
    kFaceBoundaryArray,
};

struct LoadedCheckpoint {
  vtkSmartPointer<vtkUnstructuredGrid> grid;
  Manifest manifest;
};

std::atomic_uint64_t temporary_file_sequence{};

[[nodiscard]] std::filesystem::path temporaryCheckpointPath(
    const std::filesystem::path& path) {
  const auto nonce =
      std::chrono::steady_clock::now().time_since_epoch().count();
  const auto sequence =
      temporary_file_sequence.fetch_add(1, std::memory_order_relaxed);
  return path.parent_path() /
         (path.stem().string() + ".tmp-" + std::to_string(nonce) + "-" +
          std::to_string(sequence) + path.extension().string());
}

class TemporaryCheckpointFile {
 public:
  explicit TemporaryCheckpointFile(std::filesystem::path path)
      : path_(std::move(path)) {}

  ~TemporaryCheckpointFile() {
    std::error_code error;
    std::filesystem::remove(path_, error);
  }

  TemporaryCheckpointFile(const TemporaryCheckpointFile&) = delete;
  TemporaryCheckpointFile& operator=(const TemporaryCheckpointFile&) = delete;

 private:
  std::filesystem::path path_;
};

[[nodiscard]] std::filesystem::path normalizedPath(
    const std::filesystem::path& requested_path) {
  if (requested_path.empty()) {
    throw std::invalid_argument("checkpoint path must not be empty");
  }
  auto path = std::filesystem::absolute(requested_path).lexically_normal();
  if (path.extension().empty()) {
    path.replace_extension(".vtkhdf");
  } else if (path.extension() != ".vtkhdf" && path.extension() != ".hdf") {
    throw std::invalid_argument(
        "VTKHDF checkpoint path must use .vtkhdf or .hdf");
  }
  return path;
}

[[nodiscard]] std::string unitName(
    const unit::PhysicalQuantityMetadata& metadata) {
  return metadata.unit() == units::precise::one
             ? std::string{"1"}
             : units::to_string(metadata.unit());
}

[[nodiscard]] bool isReserved(std::string_view name) {
  return std::ranges::find(kReservedNames, name) != kReservedNames.end();
}

void validateKey(std::string_view key) {
  if (key.empty()) {
    throw std::invalid_argument("checkpoint field key must not be empty");
  }
  if (key.find_first_of("./") != std::string_view::npos) {
    throw std::invalid_argument(
        "checkpoint field key must not contain '.' or '/'");
  }
  if (isReserved(key)) {
    throw std::invalid_argument("checkpoint field key is reserved");
  }
}

template <typename Array, typename Value>
void addScalar(vtkFieldData& field_data, std::string_view name, Value value) {
  vtkNew<Array> array;
  array->SetName(name.data());
  array->InsertNextValue(value);
  field_data.AddArray(array);
}

void addStringScalar(vtkFieldData& field_data, std::string_view name,
                     std::string_view value) {
  vtkNew<vtkStringArray> array;
  array->SetName(name.data());
  array->InsertNextValue(std::string{value});
  field_data.AddArray(array);
}

vtkSmartPointer<vtkDoubleArray> copyValues(std::string_view name,
                                           std::span<const double> values) {
  vtkNew<vtkDoubleArray> array;
  array->SetName(std::string{name}.c_str());
  array->SetNumberOfValues(static_cast<vtkIdType>(values.size()));
  for (std::size_t i = 0; i < values.size(); ++i) {
    array->SetValue(static_cast<vtkIdType>(i), values[i]);
  }
  return array;
}

void addFaceTopology(vtkUnstructuredGrid& grid, const mesh::IMesh& mesh) {
  const auto count = static_cast<vtkIdType>(mesh.numFaces());
  vtkNew<vtkIdTypeArray> owners;
  owners->SetName(kFaceOwnerArray.data());
  owners->SetNumberOfValues(count);
  vtkNew<vtkIdTypeArray> neighbors;
  neighbors->SetName(kFaceNeighborArray.data());
  neighbors->SetNumberOfValues(count);
  vtkNew<vtkDoubleArray> centers;
  centers->SetName(kFaceCenterArray.data());
  centers->SetNumberOfComponents(3);
  centers->SetNumberOfTuples(count);
  vtkNew<vtkDoubleArray> areas;
  areas->SetName(kFaceAreaArray.data());
  areas->SetNumberOfValues(count);
  vtkNew<vtkDoubleArray> normals;
  normals->SetName(kFaceNormalArray.data());
  normals->SetNumberOfComponents(3);
  normals->SetNumberOfTuples(count);
  vtkNew<vtkIdTypeArray> boundaries;
  boundaries->SetName(kFaceBoundaryArray.data());
  boundaries->SetNumberOfValues(count);

  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {
    const auto index = static_cast<vtkIdType>(face);
    const auto center = mesh.faceCenter(face);
    const auto normal = mesh.faceNormal(face);
    owners->SetValue(index, static_cast<vtkIdType>(mesh.owner(face)));
    neighbors->SetValue(index,
                        mesh.isBoundary(face)
                            ? vtkIdType{-1}
                            : static_cast<vtkIdType>(mesh.neighbor(face)));
    centers->SetTuple3(index, center.x, center.y, center.z);
    areas->SetValue(index, mesh.faceArea(face));
    normals->SetTuple3(index, normal.x, normal.y, normal.z);
    boundaries->SetValue(index,
                         mesh.boundaryId(face) == mesh::invalid_boundary
                             ? vtkIdType{-1}
                             : static_cast<vtkIdType>(mesh.boundaryId(face)));
  }
  grid.GetFieldData()->AddArray(owners);
  grid.GetFieldData()->AddArray(neighbors);
  grid.GetFieldData()->AddArray(centers);
  grid.GetFieldData()->AddArray(areas);
  grid.GetFieldData()->AddArray(normals);
  grid.GetFieldData()->AddArray(boundaries);
}

void addManifest(vtkUnstructuredGrid& grid, const WriteRequest& request) {
  auto& field_data = *grid.GetFieldData();
  addStringScalar(field_data, kFormatArray, kFormatName);
  addScalar<vtkTypeUInt32Array>(field_data, kVersionArray,
                                vtkhdf_format_version);
  addScalar<vtkTypeUInt64Array>(field_data, kStepArray, request.stamp.step);
  addScalar<vtkDoubleArray>(field_data, kTimeArray, request.stamp.time);
  addScalar<vtkTypeUInt64Array>(
      field_data, kNumFacesArray,
      static_cast<std::uint64_t>(request.mesh->numFaces()));

  vtkNew<vtkStringArray> associations;
  associations->SetName(kAssociationArray.data());
  vtkNew<vtkStringArray> keys;
  keys->SetName(kKeyArray.data());
  vtkNew<vtkStringArray> metadata_names;
  metadata_names->SetName(kMetadataNameArray.data());
  vtkNew<vtkStringArray> quantity_kinds;
  quantity_kinds->SetName(kQuantityKindArray.data());
  vtkNew<vtkStringArray> field_units;
  field_units->SetName(kUnitArray.data());
  vtkNew<vtkTypeUInt64Array> value_counts;
  value_counts->SetName(kValueCountArray.data());

  const auto add_field = [&](std::string_view association, std::string_view key,
                             const field::FieldMetadata& metadata,
                             std::size_t value_count) {
    associations->InsertNextValue(std::string{association});
    keys->InsertNextValue(std::string{key});
    metadata_names->InsertNextValue(metadata.name);
    if (metadata.physical_quantity.has_value()) {
      quantity_kinds->InsertNextValue(
          std::string{pemu::to_string(metadata.physical_quantity->kind())});
      field_units->InsertNextValue(unitName(*metadata.physical_quantity));
    } else {
      quantity_kinds->InsertNextValue("");
      field_units->InsertNextValue("");
    }
    value_counts->InsertNextValue(static_cast<std::uint64_t>(value_count));
  };

  for (const auto& source : request.cell_fields) {
    add_field("cell", source.key, source.field->metadata(),
              source.field->size());
  }
  for (const auto& source : request.face_fields) {
    add_field("face", source.key, source.field->metadata(),
              source.field->size());
  }
  for (const auto& source : request.scalars) {
    add_field("scalar", source.key, source.metadata, 1);
  }

  field_data.AddArray(associations);
  field_data.AddArray(keys);
  field_data.AddArray(metadata_names);
  field_data.AddArray(quantity_kinds);
  field_data.AddArray(field_units);
  field_data.AddArray(value_counts);
}

template <typename Source>
void validateSources(const mesh::IMesh& mesh, std::span<const Source> sources,
                     std::unordered_set<std::string>& keys) {
  for (const auto& source : sources) {
    if (source.field == nullptr) {
      throw std::invalid_argument("checkpoint field source must not be null");
    }
    if (&source.field->mesh() != &mesh) {
      throw std::invalid_argument(
          "checkpoint field belongs to a different mesh");
    }
    validateKey(source.key);
    if (!keys.insert(source.key).second) {
      throw std::invalid_argument("duplicate checkpoint field key");
    }
  }
}

void validateScalarSources(std::span<const ScalarSource> sources,
                           std::unordered_set<std::string>& keys) {
  for (const auto& source : sources) {
    if (source.value == nullptr) {
      throw std::invalid_argument("checkpoint scalar source must not be null");
    }
    if (!std::isfinite(*source.value)) {
      throw std::invalid_argument("checkpoint scalar must be finite");
    }
    validateKey(source.key);
    if (!keys.insert(source.key).second) {
      throw std::invalid_argument("duplicate checkpoint scalar key");
    }
  }
}

template <typename Value>
[[nodiscard]] Value readNumericScalar(vtkFieldData& field_data,
                                      std::string_view name) {
  auto* array = field_data.GetArray(std::string{name}.c_str());
  if (array == nullptr || array->GetNumberOfValues() != 1) {
    throw std::runtime_error("invalid checkpoint scalar: " + std::string{name});
  }
  if constexpr (std::is_floating_point_v<Value>) {
    return static_cast<Value>(array->GetTuple1(0));
  } else {
    const auto value = array->GetVariantValue(0).ToUnsignedLongLong();
    if (value > std::numeric_limits<Value>::max()) {
      throw std::overflow_error("checkpoint scalar is out of range: " +
                                std::string{name});
    }
    return static_cast<Value>(value);
  }
}

[[nodiscard]] std::string readStringScalar(vtkFieldData& field_data,
                                           std::string_view name) {
  auto* array = vtkStringArray::SafeDownCast(
      field_data.GetAbstractArray(std::string{name}.c_str()));
  if (array == nullptr || array->GetNumberOfValues() != 1) {
    throw std::runtime_error("invalid checkpoint string: " + std::string{name});
  }
  return array->GetValue(0);
}

[[nodiscard]] unit::QuantityKind parseQuantityKind(std::string_view value) {
  for (unsigned int raw = 0;
       raw <=
       static_cast<unsigned int>(unit::QuantityKind::normal_drift_velocity);
       ++raw) {
    const auto kind = static_cast<unit::QuantityKind>(raw);
    if (pemu::to_string(kind) == value) {
      return kind;
    }
  }
  throw std::runtime_error("invalid checkpoint quantity kind");
}

[[nodiscard]] units::precise_unit parseUnit(std::string_view value) {
  if (value == "1") {
    return units::precise::one;
  }
  auto parsed = units::unit_from_string(std::string{value});
  if (!units::is_valid(parsed) || units::is_error(parsed)) {
    throw std::runtime_error("invalid checkpoint unit");
  }
  return parsed;
}

[[nodiscard]] LoadedCheckpoint loadCheckpoint(
    const std::filesystem::path& requested_path) {
  const auto path = normalizedPath(requested_path);
  if (!std::filesystem::is_regular_file(path)) {
    throw std::filesystem::filesystem_error(
        "checkpoint file does not exist", path,
        std::make_error_code(std::errc::no_such_file_or_directory));
  }

  vtkNew<vtkHDFReader> reader;
  reader->SetFileName(path.string().c_str());
  reader->Update();
  auto* output =
      vtkUnstructuredGrid::SafeDownCast(reader->GetOutputDataObject(0));
  if (output == nullptr || output->GetFieldData() == nullptr) {
    throw std::runtime_error("file is not an unstructured-grid checkpoint");
  }

  LoadedCheckpoint loaded;
  loaded.grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  loaded.grid->DeepCopy(output);
  auto& field_data = *loaded.grid->GetFieldData();
  if (readStringScalar(field_data, kFormatArray) != kFormatName) {
    throw std::runtime_error("file is not a pemu VTKHDF checkpoint");
  }

  loaded.manifest.format_version =
      readNumericScalar<std::uint32_t>(field_data, kVersionArray);
  if (loaded.manifest.format_version != vtkhdf_format_version) {
    throw std::runtime_error("unsupported checkpoint format version");
  }
  loaded.manifest.stamp.step =
      readNumericScalar<std::uint64_t>(field_data, kStepArray);
  loaded.manifest.stamp.time =
      readNumericScalar<double>(field_data, kTimeArray);
  if (!std::isfinite(loaded.manifest.stamp.time)) {
    throw std::runtime_error("checkpoint time is not finite");
  }
  loaded.manifest.num_cells =
      static_cast<std::size_t>(loaded.grid->GetNumberOfCells());
  loaded.manifest.num_faces = static_cast<std::size_t>(
      readNumericScalar<std::uint64_t>(field_data, kNumFacesArray));

  auto* associations = vtkStringArray::SafeDownCast(
      field_data.GetAbstractArray(kAssociationArray.data()));
  auto* keys = vtkStringArray::SafeDownCast(
      field_data.GetAbstractArray(kKeyArray.data()));
  auto* metadata_names = vtkStringArray::SafeDownCast(
      field_data.GetAbstractArray(kMetadataNameArray.data()));
  auto* quantity_kinds = vtkStringArray::SafeDownCast(
      field_data.GetAbstractArray(kQuantityKindArray.data()));
  auto* field_units = vtkStringArray::SafeDownCast(
      field_data.GetAbstractArray(kUnitArray.data()));
  auto* value_counts = field_data.GetArray(kValueCountArray.data());
  if (associations == nullptr || keys == nullptr || metadata_names == nullptr ||
      quantity_kinds == nullptr || field_units == nullptr ||
      value_counts == nullptr) {
    throw std::runtime_error("checkpoint field manifest is incomplete");
  }

  const vtkIdType count = associations->GetNumberOfValues();
  if (keys->GetNumberOfValues() != count ||
      metadata_names->GetNumberOfValues() != count ||
      quantity_kinds->GetNumberOfValues() != count ||
      field_units->GetNumberOfValues() != count ||
      value_counts->GetNumberOfValues() != count) {
    throw std::runtime_error("checkpoint field manifest lengths differ");
  }

  std::unordered_set<std::string> unique_fields;
  std::unordered_set<std::string> unique_field_data_keys;
  loaded.manifest.fields.reserve(static_cast<std::size_t>(count));
  for (vtkIdType i = 0; i < count; ++i) {
    const std::string association = associations->GetValue(i);
    const std::string key = keys->GetValue(i);
    validateKey(key);
    FieldDescriptor descriptor;
    if (association == "cell") {
      descriptor.association = FieldAssociation::cell;
    } else if (association == "face") {
      descriptor.association = FieldAssociation::face;
    } else if (association == "scalar") {
      descriptor.association = FieldAssociation::scalar;
    } else {
      throw std::runtime_error("invalid checkpoint field association");
    }
    const std::string unique_key = association + ":" + key;
    if (!unique_fields.insert(unique_key).second) {
      throw std::runtime_error("duplicate checkpoint field manifest entry");
    }
    if (descriptor.association != FieldAssociation::cell &&
        !unique_field_data_keys.insert(key).second) {
      throw std::runtime_error(
          "checkpoint face and scalar storage keys collide");
    }
    descriptor.key = key;
    descriptor.metadata.name = metadata_names->GetValue(i);
    const std::string quantity = quantity_kinds->GetValue(i);
    const std::string unit_name = field_units->GetValue(i);
    if (quantity.empty() != unit_name.empty()) {
      throw std::runtime_error("incomplete checkpoint physical metadata");
    }
    if (!quantity.empty()) {
      descriptor.metadata.physical_quantity = unit::PhysicalQuantityMetadata{
          parseQuantityKind(quantity), parseUnit(unit_name)};
    }
    const auto raw_count =
        value_counts->GetVariantValue(i).ToUnsignedLongLong();
    if (raw_count > std::numeric_limits<std::size_t>::max()) {
      throw std::overflow_error("checkpoint field is too large");
    }
    descriptor.value_count = static_cast<std::size_t>(raw_count);

    vtkDataArray* values =
        descriptor.association == FieldAssociation::cell
            ? loaded.grid->GetCellData()->GetArray(key.c_str())
            : field_data.GetArray(key.c_str());
    if (vtkDoubleArray::SafeDownCast(values) == nullptr ||
        values->GetNumberOfValues() !=
            static_cast<vtkIdType>(descriptor.value_count)) {
      throw std::runtime_error("checkpoint field data does not match manifest");
    }
    if (descriptor.association == FieldAssociation::scalar &&
        !std::isfinite(values->GetTuple1(0))) {
      throw std::runtime_error("checkpoint scalar is not finite");
    }
    const std::size_t expected =
        descriptor.association == FieldAssociation::cell
            ? loaded.manifest.num_cells
        : descriptor.association == FieldAssociation::face
            ? loaded.manifest.num_faces
            : 1;
    if (descriptor.value_count != expected) {
      throw std::runtime_error("checkpoint field has invalid association size");
    }
    loaded.manifest.fields.push_back(std::move(descriptor));
  }
  return loaded;
}

[[nodiscard]] bool topologyMatches(const mesh::IMesh& mesh,
                                   vtkUnstructuredGrid& grid) {
  if (grid.GetNumberOfPoints() != static_cast<vtkIdType>(mesh.numVertices()) ||
      grid.GetNumberOfCells() != static_cast<vtkIdType>(mesh.numCells())) {
    return false;
  }
  for (mesh::VertexId vertex = 0; vertex < mesh.numVertices(); ++vertex) {
    double point[3]{};
    grid.GetPoint(static_cast<vtkIdType>(vertex), point);
    const auto expected = mesh.vertex(vertex);
    if (point[0] != expected.x || point[1] != expected.y ||
        point[2] != expected.z) {
      return false;
    }
  }
  for (mesh::CellId cell = 0; cell < mesh.numCells(); ++cell) {
    const auto expected = mesh.cellVertices(cell);
    auto* ids = grid.GetCell(static_cast<vtkIdType>(cell))->GetPointIds();
    if (ids->GetNumberOfIds() != static_cast<vtkIdType>(expected.size())) {
      return false;
    }
    for (std::size_t i = 0; i < expected.size(); ++i) {
      if (ids->GetId(static_cast<vtkIdType>(i)) !=
          static_cast<vtkIdType>(expected[i])) {
        return false;
      }
    }
  }
  auto* owners = grid.GetFieldData()->GetArray(kFaceOwnerArray.data());
  auto* neighbors = grid.GetFieldData()->GetArray(kFaceNeighborArray.data());
  auto* centers = grid.GetFieldData()->GetArray(kFaceCenterArray.data());
  auto* areas = grid.GetFieldData()->GetArray(kFaceAreaArray.data());
  auto* normals = grid.GetFieldData()->GetArray(kFaceNormalArray.data());
  auto* boundaries = grid.GetFieldData()->GetArray(kFaceBoundaryArray.data());
  const auto face_count = static_cast<vtkIdType>(mesh.numFaces());
  if (owners == nullptr || neighbors == nullptr || centers == nullptr ||
      areas == nullptr || normals == nullptr || boundaries == nullptr ||
      owners->GetNumberOfValues() != face_count ||
      neighbors->GetNumberOfValues() != face_count ||
      centers->GetNumberOfTuples() != face_count ||
      centers->GetNumberOfComponents() != 3 ||
      areas->GetNumberOfValues() != face_count ||
      normals->GetNumberOfTuples() != face_count ||
      normals->GetNumberOfComponents() != 3 ||
      boundaries->GetNumberOfValues() != face_count) {
    return false;
  }
  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {
    const auto index = static_cast<vtkIdType>(face);
    const auto center = mesh.faceCenter(face);
    const auto normal = mesh.faceNormal(face);
    double stored_center[3]{};
    double stored_normal[3]{};
    centers->GetTuple(index, stored_center);
    normals->GetTuple(index, stored_normal);
    const auto expected_neighbor =
        mesh.isBoundary(face) ? vtkIdType{-1}
                              : static_cast<vtkIdType>(mesh.neighbor(face));
    const auto expected_boundary =
        mesh.boundaryId(face) == mesh::invalid_boundary
            ? vtkIdType{-1}
            : static_cast<vtkIdType>(mesh.boundaryId(face));
    if (static_cast<vtkIdType>(owners->GetTuple1(index)) != mesh.owner(face) ||
        static_cast<vtkIdType>(neighbors->GetTuple1(index)) !=
            expected_neighbor ||
        stored_center[0] != center.x || stored_center[1] != center.y ||
        stored_center[2] != center.z ||
        areas->GetTuple1(index) != mesh.faceArea(face) ||
        stored_normal[0] != normal.x || stored_normal[1] != normal.y ||
        stored_normal[2] != normal.z ||
        static_cast<vtkIdType>(boundaries->GetTuple1(index)) !=
            expected_boundary) {
      return false;
    }
  }
  return true;
}

[[nodiscard]] const FieldDescriptor& findDescriptor(
    const Manifest& manifest, FieldAssociation association,
    std::string_view key) {
  const auto iterator = std::ranges::find_if(
      manifest.fields, [association, key](const FieldDescriptor& descriptor) {
        return descriptor.association == association && descriptor.key == key;
      });
  if (iterator == manifest.fields.end()) {
    throw std::invalid_argument("checkpoint does not contain requested field");
  }
  return *iterator;
}

[[nodiscard]] bool metadataMatches(const field::FieldMetadata& lhs,
                                   const field::FieldMetadata& rhs) {
  return lhs.name == rhs.name && lhs.physical_quantity == rhs.physical_quantity;
}

template <typename Target>
void validateTargets(const RestoreRequest& request,
                     std::span<const Target> targets,
                     FieldAssociation association,
                     const LoadedCheckpoint& loaded,
                     std::unordered_set<std::string>& keys) {
  for (const auto& target : targets) {
    if (target.field == nullptr) {
      throw std::invalid_argument("checkpoint field target must not be null");
    }
    validateKey(target.key);
    const std::string unique_key =
        std::to_string(static_cast<unsigned int>(association)) + ":" +
        target.key;
    if (!keys.insert(unique_key).second) {
      throw std::invalid_argument("duplicate checkpoint restore target");
    }
    if (&target.field->mesh() != request.mesh) {
      throw std::invalid_argument(
          "checkpoint restore target belongs to a different mesh");
    }
    const auto& descriptor =
        findDescriptor(loaded.manifest, association, target.key);
    if (target.field->size() != descriptor.value_count) {
      throw std::invalid_argument("checkpoint restore target size differs");
    }
    if (!metadataMatches(target.field->metadata(), descriptor.metadata)) {
      throw std::invalid_argument("checkpoint restore metadata differs");
    }
  }
}

void copyCellTargets(std::span<const CellFieldTarget> targets,
                     vtkUnstructuredGrid& grid) {
  for (const auto& target : targets) {
    auto* values = vtkDoubleArray::SafeDownCast(
        grid.GetCellData()->GetArray(target.key.c_str()));
    for (mesh::CellId cell = 0; cell < target.field->size(); ++cell) {
      (*target.field)[cell] = values->GetValue(static_cast<vtkIdType>(cell));
    }
  }
}

void copyFaceTargets(std::span<const FaceFieldTarget> targets,
                     vtkUnstructuredGrid& grid) {
  for (const auto& target : targets) {
    auto* values = vtkDoubleArray::SafeDownCast(
        grid.GetFieldData()->GetArray(target.key.c_str()));
    for (mesh::FaceId face = 0; face < target.field->size(); ++face) {
      (*target.field)[face] = values->GetValue(static_cast<vtkIdType>(face));
    }
  }
}

void validateScalarTargets(const RestoreRequest& request,
                           const LoadedCheckpoint& loaded,
                           std::unordered_set<std::string>& keys) {
  for (const auto& target : request.scalars) {
    if (target.value == nullptr) {
      throw std::invalid_argument("checkpoint scalar target must not be null");
    }
    validateKey(target.key);
    const std::string unique_key =
        std::to_string(static_cast<unsigned int>(FieldAssociation::scalar)) +
        ":" + target.key;
    if (!keys.insert(unique_key).second) {
      throw std::invalid_argument("duplicate checkpoint restore target");
    }
    const auto& descriptor =
        findDescriptor(loaded.manifest, FieldAssociation::scalar, target.key);
    if (descriptor.value_count != 1) {
      throw std::invalid_argument("checkpoint scalar size differs");
    }
    if (!metadataMatches(target.metadata, descriptor.metadata)) {
      throw std::invalid_argument("checkpoint scalar metadata differs");
    }
  }
}

void copyScalarTargets(std::span<const ScalarTarget> targets,
                       vtkUnstructuredGrid& grid) {
  for (const auto& target : targets) {
    auto* values = vtkDoubleArray::SafeDownCast(
        grid.GetFieldData()->GetArray(target.key.c_str()));
    *target.value = values->GetValue(0);
  }
}

}  // namespace

OutputRecord VtkHdfWriter::write(const WriteRequest& request) const {
  if (request.mesh == nullptr) {
    throw std::invalid_argument("checkpoint requires a mesh");
  }
  if (!std::isfinite(request.stamp.time)) {
    throw std::invalid_argument("checkpoint time must be finite");
  }
  const auto path = normalizedPath(request.path);
  if (std::filesystem::exists(path) && !request.overwrite) {
    throw std::filesystem::filesystem_error(
        "refusing to overwrite checkpoint", path,
        std::make_error_code(std::errc::file_exists));
  }
  if (!path.parent_path().empty()) {
    std::filesystem::create_directories(path.parent_path());
  }

  std::unordered_set<std::string> cell_keys;
  std::unordered_set<std::string> field_data_keys;
  validateSources(*request.mesh, request.cell_fields, cell_keys);
  validateSources(*request.mesh, request.face_fields, field_data_keys);
  validateScalarSources(request.scalars, field_data_keys);
  auto grid = common::toVtkUnstructuredGrid(*request.mesh);
  addFaceTopology(*grid, *request.mesh);
  for (const auto& source : request.cell_fields) {
    grid->GetCellData()->AddArray(copyValues(source.key, source.field->span()));
  }
  for (const auto& source : request.face_fields) {
    grid->GetFieldData()->AddArray(
        copyValues(source.key, source.field->span()));
  }
  for (const auto& source : request.scalars) {
    grid->GetFieldData()->AddArray(
        copyValues(source.key, std::span<const double>{source.value, 1}));
  }
  addManifest(*grid, request);

  const auto temporary_path = temporaryCheckpointPath(path);
  TemporaryCheckpointFile temporary_file{temporary_path};
  vtkNew<vtkHDFWriter> writer;
  writer->SetFileName(temporary_path.string().c_str());
  writer->SetOverwrite(false);
  writer->SetWriteAllTimeSteps(false);
  writer->SetInputData(grid);
  const auto success = writer->Write();
  if (success == 0 || writer->GetErrorCode() != vtkErrorCode::NoError ||
      !std::filesystem::is_regular_file(temporary_path)) {
    throw std::runtime_error("vtkHDFWriter failed to create checkpoint");
  }
  if (request.overwrite) {
    std::filesystem::rename(temporary_path, path);
  } else {
    std::filesystem::create_hard_link(temporary_path, path);
  }
  if (!std::filesystem::is_regular_file(path)) {
    throw std::runtime_error("failed to commit checkpoint atomically");
  }
  return {.path = path, .stamp = request.stamp};
}

Manifest VtkHdfReader::inspect(const std::filesystem::path& path) const {
  return loadCheckpoint(path).manifest;
}

OutputRecord VtkHdfReader::restore(const std::filesystem::path& requested_path,
                                   const RestoreRequest& request) const {
  if (request.mesh == nullptr) {
    throw std::invalid_argument("checkpoint restore requires a mesh");
  }
  auto loaded = loadCheckpoint(requested_path);
  if (loaded.manifest.num_faces != request.mesh->numFaces() ||
      !topologyMatches(*request.mesh, *loaded.grid)) {
    throw std::invalid_argument("checkpoint mesh topology differs");
  }
  if (request.require_all_fields && request.cell_fields.size() +
                                            request.face_fields.size() +
                                            request.scalars.size() !=
                                        loaded.manifest.fields.size()) {
    throw std::invalid_argument(
        "checkpoint restore does not target every stored field");
  }

  std::unordered_set<std::string> keys;
  validateTargets(request, request.cell_fields, FieldAssociation::cell, loaded,
                  keys);
  validateTargets(request, request.face_fields, FieldAssociation::face, loaded,
                  keys);
  validateScalarTargets(request, loaded, keys);
  copyCellTargets(request.cell_fields, *loaded.grid);
  copyFaceTargets(request.face_fields, *loaded.grid);
  copyScalarTargets(request.scalars, *loaded.grid);
  return {.path = normalizedPath(requested_path),
          .stamp = loaded.manifest.stamp};
}

}  // namespace pemu::output::checkpoint
