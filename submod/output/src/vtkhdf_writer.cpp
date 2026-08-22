#include <pemu/output/dump/vtkhdf_writer.hpp>

#include <pemu/output/common/mesh_adapter.hpp>
#include <pemu/unit/quantity_metadata.hpp>

#include <llnl-units/units.hpp>

#include <vtkCellData.h>
#include <vtkDataObject.h>
#include <vtkDoubleArray.h>
#include <vtkErrorCode.h>
#include <vtkFieldData.h>
#include <vtkHDFWriter.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkStringArray.h>
#include <vtkTypeUInt64Array.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridAlgorithm.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

namespace pemu::output::dump {

namespace {

constexpr std::string_view kFieldAssociationArrayName =
    "pemu_field_association";
constexpr std::string_view kFieldNameArrayName = "pemu_field_name";
constexpr std::string_view kFieldQuantityKindArrayName =
    "pemu_field_quantity_kind";
constexpr std::string_view kFieldUnitArrayName = "pemu_field_unit";
constexpr std::string_view kFaceOwnerArrayName = "pemu_face_owner";
constexpr std::string_view kFaceNeighborArrayName = "pemu_face_neighbor";
constexpr std::string_view kFaceCenterArrayName = "pemu_face_center";
constexpr std::string_view kFaceAreaArrayName = "pemu_face_area";
constexpr std::string_view kFaceBoundaryArrayName = "pemu_face_boundary_id";
constexpr std::string_view kStepArrayName = "pemu_step";
constexpr std::string_view kTimeArrayName = "pemu_time";

struct MetadataRow {
  std::string association;
  std::string name;
  std::string quantity_kind;
  std::string unit;

  /** @brief Compares serialized field-schema rows component by component. */
  [[nodiscard]] bool operator==(const MetadataRow&) const = default;
};

struct NamedDoubleBuffer {
  std::string name;
  std::vector<double> values;
};

struct TemporalSnapshot {
  OutputStamp stamp;
  std::vector<NamedDoubleBuffer> cell_arrays;
  std::vector<NamedDoubleBuffer> field_arrays;
};

struct CapturedSnapshot {
  TemporalSnapshot snapshot;
  std::vector<MetadataRow> schema;
  bool has_face_fields{};
};

/** @brief Converts runtime unit metadata to its stable serialized name. */
std::string unitName(const unit::PhysicalQuantityMetadata& metadata) {
  return metadata.unit() == units::precise::one
             ? std::string{"1"}
             : units::to_string(metadata.unit());
}

/** @brief Reports whether a field name collides with pemu VTKHDF metadata. */
bool isReservedFieldDataName(std::string_view name) {
  return name == kFieldAssociationArrayName || name == kFieldNameArrayName ||
         name == kFieldQuantityKindArrayName || name == kFieldUnitArrayName ||
         name == kFaceOwnerArrayName || name == kFaceNeighborArrayName ||
         name == kFaceCenterArrayName || name == kFaceAreaArrayName ||
         name == kFaceBoundaryArrayName || name == kStepArrayName ||
         name == kTimeArrayName;
}

/** @brief Creates a stable schema row from field metadata. */
MetadataRow makeMetadataRow(std::string_view association,
                            const field::FieldMetadata& metadata) {
  MetadataRow row{.association = std::string{association},
                  .name = metadata.name};
  if (metadata.physical_quantity.has_value()) {
    row.quantity_kind =
        std::string{pemu::to_string(metadata.physical_quantity->kind())};
    row.unit = unitName(*metadata.physical_quantity);
  }
  return row;
}

/**
 * @brief Resolves output metadata without changing the source field.
 *
 * The VTK array name is an output concern, so it is reflected in fallback
 * metadata.  An explicit selection override instead preserves its semantic
 * metadata name while supplying its own physical-quantity description.
 */
field::FieldMetadata resolveOutputMetadata(
    const field::FieldMetadata& source_metadata, std::string_view output_name,
    const FieldMetadata& output_metadata) {
  auto metadata = source_metadata;
  if (output_metadata.meta_name.empty()) {
    metadata.name = std::string{output_name};
    return metadata;
  }

  metadata.name = output_metadata.meta_name;
  metadata.physical_quantity = unit::PhysicalQuantityMetadata{
      output_metadata.quantity_kind, output_metadata.unit};
  return metadata;
}

/** @brief Validates a cell field against its output mesh. */
void validateField(const mesh::IMesh& mesh,
                   const field::CellField<double>& values) {
  if (&values.mesh() != &mesh) {
    throw std::invalid_argument("cell field belongs to a different mesh");
  }
  if (values.metadata().name.empty()) {
    throw std::invalid_argument("cell field name must not be empty");
  }
}

/** @brief Validates a face field against its output mesh. */
void validateField(const mesh::IMesh& mesh,
                   const field::FaceField<double>& values) {
  if (&values.mesh() != &mesh) {
    throw std::invalid_argument("face field belongs to a different mesh");
  }
  if (values.metadata().name.empty()) {
    throw std::invalid_argument("face field name must not be empty");
  }
}

/** @brief Copies one mutable simulation field into temporal snapshot storage. */
NamedDoubleBuffer captureValues(std::string name,
                                std::span<const double> values) {
  return {.name = std::move(name),
          .values = std::vector<double>{values.begin(), values.end()}};
}

/** @brief Adds mesh-owned face topology once to the reusable VTK grid. */
void addFaceTopology(vtkUnstructuredGrid& grid, const mesh::IMesh& mesh) {
  const auto num_faces = static_cast<vtkIdType>(mesh.numFaces());

  vtkNew<vtkIdTypeArray> owners;
  owners->SetName(kFaceOwnerArrayName.data());
  owners->SetNumberOfValues(num_faces);
  vtkNew<vtkIdTypeArray> neighbors;
  neighbors->SetName(kFaceNeighborArrayName.data());
  neighbors->SetNumberOfValues(num_faces);
  vtkNew<vtkDoubleArray> centers;
  centers->SetName(kFaceCenterArrayName.data());
  centers->SetNumberOfComponents(3);
  centers->SetNumberOfTuples(num_faces);
  vtkNew<vtkDoubleArray> areas;
  areas->SetName(kFaceAreaArrayName.data());
  areas->SetNumberOfValues(num_faces);
  vtkNew<vtkIdTypeArray> boundary_ids;
  boundary_ids->SetName(kFaceBoundaryArrayName.data());
  boundary_ids->SetNumberOfValues(num_faces);

  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {
    const auto index = static_cast<vtkIdType>(face);
    const auto center = mesh.faceCenter(face);
    owners->SetValue(index, static_cast<vtkIdType>(mesh.owner(face)));
    neighbors->SetValue(index,
                        mesh.isBoundary(face)
                            ? vtkIdType{-1}
                            : static_cast<vtkIdType>(mesh.neighbor(face)));
    centers->SetTuple3(index, center.x, center.y, center.z);
    areas->SetValue(index, mesh.faceArea(face));
    boundary_ids->SetValue(index,
                           mesh.boundaryId(face) == mesh::invalid_boundary
                               ? vtkIdType{-1}
                               : static_cast<vtkIdType>(mesh.boundaryId(face)));
  }

  grid.GetFieldData()->AddArray(owners);
  grid.GetFieldData()->AddArray(neighbors);
  grid.GetFieldData()->AddArray(centers);
  grid.GetFieldData()->AddArray(areas);
  grid.GetFieldData()->AddArray(boundary_ids);
}

/** @brief Computes a cell visualization buffer from one face field. */
std::vector<double> makeCellCenteredFaceValues(
    const mesh::IMesh& mesh, const field::FaceField<double>& face_field) {
  std::vector<double> result(mesh.numCells());
  for (mesh::CellId cell = 0; cell < mesh.numCells(); ++cell) {
    double weighted_sum = 0.0;
    double area_sum = 0.0;
    for (const mesh::FaceId face : mesh.cellFaces(cell)) {
      const double area = mesh.faceArea(face);
      weighted_sum += face_field[face] * area;
      area_sum += area;
    }
    if (!(area_sum > 0.0) || !std::isfinite(area_sum)) {
      throw std::runtime_error("cannot center face field on invalid cell");
    }
    result[cell] = weighted_sum / area_sum;
  }
  return result;
}

/** @brief Normalizes and validates a VTKHDF destination path. */
std::filesystem::path normalizedOutputPath(
    const std::filesystem::path& requested_path) {
  if (requested_path.empty()) {
    throw std::invalid_argument("VTKHDF output path must not be empty");
  }
  auto path = std::filesystem::absolute(requested_path).lexically_normal();
  if (path.extension().empty()) {
    path.replace_extension(".vtkhdf");
  } else if (path.extension() != ".vtkhdf" && path.extension() != ".hdf") {
    throw std::invalid_argument("VTKHDF output path must use .vtkhdf or .hdf");
  }
  return path;
}

/** @brief Adds immutable field-schema metadata arrays to a reusable grid. */
void addSchemaMetadata(vtkUnstructuredGrid& grid,
                       const std::vector<MetadataRow>& schema) {
  vtkNew<vtkStringArray> associations;
  associations->SetName(kFieldAssociationArrayName.data());
  vtkNew<vtkStringArray> names;
  names->SetName(kFieldNameArrayName.data());
  vtkNew<vtkStringArray> quantity_kinds;
  quantity_kinds->SetName(kFieldQuantityKindArrayName.data());
  vtkNew<vtkStringArray> field_units;
  field_units->SetName(kFieldUnitArrayName.data());

  for (const auto& row : schema) {
    associations->InsertNextValue(row.association);
    names->InsertNextValue(row.name);
    quantity_kinds->InsertNextValue(row.quantity_kind);
    field_units->InsertNextValue(row.unit);
  }

  grid.GetFieldData()->AddArray(associations);
  grid.GetFieldData()->AddArray(names);
  grid.GetFieldData()->AddArray(quantity_kinds);
  grid.GetFieldData()->AddArray(field_units);
}

/** @brief Captures one request while preserving its temporal field schema. */
CapturedSnapshot captureSnapshot(const mesh::IMesh& mesh,
                                 const Request& request) {
  CapturedSnapshot captured;
  captured.snapshot.stamp = request.stamp;
  std::unordered_set<std::string> cell_names;

  const auto add_cell_field = [&](const field::CellField<double>& values,
                                  std::string output_name,
                                  field::FieldMetadata metadata) {
    validateField(mesh, values);
    if (output_name.empty()) {
      throw std::invalid_argument("cell output field name must not be empty");
    }
    if (!cell_names.insert(output_name).second) {
      throw std::invalid_argument("duplicate cell field name");
    }
    captured.snapshot.cell_arrays.push_back(
        captureValues(std::move(output_name), values.span()));
    captured.schema.push_back(makeMetadataRow("cell", metadata));
  };

  for (const auto& field_reference : request.cell_fields) {
    const auto& values = field_reference.get();
    add_cell_field(values, values.metadata().name, values.metadata());
  }
  for (const auto& selection : request.cell_field_selections) {
    if (selection.field == nullptr) {
      throw std::invalid_argument("cell field selection must not be null");
    }
    const std::string output_name = selection.name.empty()
                                        ? selection.field->metadata().name
                                        : selection.name;
    add_cell_field(*selection.field, output_name,
                   resolveOutputMetadata(selection.field->metadata(),
                                         output_name, selection.meta_data));
  }

  std::unordered_set<std::string> face_names;
  for (const auto& selection : request.face_fields) {
    if (selection.field == nullptr) {
      throw std::invalid_argument("face field selection must not be null");
    }
    const auto& values = *selection.field;
    validateField(mesh, values);
    const std::string output_name =
        selection.name.empty() ? values.metadata().name : selection.name;
    const auto metadata = resolveOutputMetadata(values.metadata(), output_name,
                                                selection.meta_data);
    if (isReservedFieldDataName(output_name)) {
      throw std::invalid_argument(
          "face field name is reserved by VTKHDF output");
    }
    if (!face_names.insert(output_name).second) {
      throw std::invalid_argument("duplicate face field name");
    }

    captured.has_face_fields = true;
    captured.snapshot.field_arrays.push_back(
        captureValues(output_name, values.span()));
    captured.schema.push_back(makeMetadataRow("face", metadata));

    if (selection.include_cell_centered_visualization) {
      const std::string centered_name = selection.cell_centered_name.empty()
                                            ? output_name + "_cell_centered"
                                            : selection.cell_centered_name;
      if (!cell_names.insert(centered_name).second) {
        throw std::invalid_argument("duplicate cell visualization field name");
      }
      captured.snapshot.cell_arrays.push_back(
          {.name = centered_name,
           .values = makeCellCenteredFaceValues(mesh, values)});
      auto centered_metadata = metadata;
      if (selection.meta_data.meta_name.empty()) {
        centered_metadata.name = centered_name;
      }
      captured.schema.push_back(
          makeMetadataRow("cell_visualization", centered_metadata));
    }
  }
  return captured;
}

/** @brief Wraps snapshot memory in a VTK array without copying its values. */
vtkSmartPointer<vtkDoubleArray> wrapDoubleBuffer(
    const NamedDoubleBuffer& buffer) {
  vtkNew<vtkDoubleArray> array;
  array->SetName(buffer.name.c_str());
  array->SetArray(const_cast<double*>(buffer.values.data()),
                  static_cast<vtkIdType>(buffer.values.size()), 1);
  return array;
}

class TemporalUnstructuredGridSource final
    : public vtkUnstructuredGridAlgorithm {
 public:
  /** @brief Allocates a temporal unstructured-grid source. */
  static TemporalUnstructuredGridSource* New();
  vtkTypeMacro(TemporalUnstructuredGridSource, vtkUnstructuredGridAlgorithm);

  /**
   * @brief Binds reusable topology and immutable snapshot buffers.
   * @param topology Grid whose points and connectivity stay unchanged.
   * @param snapshots Ordered temporal snapshot storage.
   */
  void setSeries(vtkUnstructuredGrid* topology,
                 const std::vector<TemporalSnapshot>* snapshots) {
    topology_ = topology;
    snapshots_ = snapshots;
    Modified();
  }

 protected:
  /** @brief Creates a source with no input ports. */
  TemporalUnstructuredGridSource() { SetNumberOfInputPorts(0); }

  /** @brief Publishes the available simulation times to the VTK pipeline. */
  int RequestInformation(vtkInformation*, vtkInformationVector**,
                         vtkInformationVector* output_vector) override {
    if (snapshots_ == nullptr || snapshots_->empty()) {
      return 0;
    }
    times_.clear();
    times_.reserve(snapshots_->size());
    for (const auto& snapshot : *snapshots_) {
      times_.push_back(snapshot.stamp.time);
    }
    auto* output_info = output_vector->GetInformationObject(0);
    output_info->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
                     times_.data(), static_cast<int>(times_.size()));
    const double range[2] = {times_.front(), times_.back()};
    output_info->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), range, 2);
    return 1;
  }

  /** @brief Presents one requested snapshot while sharing the static mesh. */
  int RequestData(vtkInformation*, vtkInformationVector**,
                  vtkInformationVector* output_vector) override {
    if (topology_ == nullptr || snapshots_ == nullptr || snapshots_->empty()) {
      return 0;
    }
    auto* output_info = output_vector->GetInformationObject(0);
    const double requested_time =
        output_info->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())
            ? output_info->Get(
                  vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())
            : snapshots_->front().stamp.time;
    const auto iterator =
        std::lower_bound(snapshots_->begin(), snapshots_->end(), requested_time,
                         [](const TemporalSnapshot& snapshot, double time) {
                           return snapshot.stamp.time < time;
                         });
    const auto& snapshot =
        iterator == snapshots_->end() ? snapshots_->back() : *iterator;

    auto* output = vtkUnstructuredGrid::GetData(output_vector, 0);
    output->ShallowCopy(topology_);
    for (const auto& buffer : snapshot.cell_arrays) {
      output->GetCellData()->AddArray(wrapDoubleBuffer(buffer));
    }
    for (const auto& buffer : snapshot.field_arrays) {
      output->GetFieldData()->AddArray(wrapDoubleBuffer(buffer));
    }

    vtkNew<vtkTypeUInt64Array> step;
    step->SetName(kStepArrayName.data());
    step->InsertNextValue(snapshot.stamp.step);
    output->GetFieldData()->AddArray(step);
    vtkNew<vtkDoubleArray> time;
    time->SetName(kTimeArrayName.data());
    time->InsertNextValue(snapshot.stamp.time);
    output->GetFieldData()->AddArray(time);
    output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(),
                                  snapshot.stamp.time);
    return 1;
  }

 private:
  vtkSmartPointer<vtkUnstructuredGrid> topology_;
  const std::vector<TemporalSnapshot>* snapshots_{};
  std::vector<double> times_;
};

vtkStandardNewMacro(TemporalUnstructuredGridSource);

class VtkHdfOutputSeries final : public ISeries {
 public:
  /**
   * @brief Initializes one series and converts its mesh topology exactly once.
   * @param request Series mesh, output path, and overwrite policy.
   */
  explicit VtkHdfOutputSeries(const SeriesRequest& request)
      : mesh_(request.mesh),
        path_(normalizedOutputPath(request.path)),
        overwrite_(request.overwrite) {
    if (mesh_ == nullptr) {
      throw std::invalid_argument("field output series requires a mesh");
    }
    if (std::filesystem::exists(path_) && !overwrite_) {
      throw std::filesystem::filesystem_error(
          "refusing to overwrite VTKHDF output", path_,
          std::make_error_code(std::errc::file_exists));
    }
    if (!path_.parent_path().empty()) {
      std::filesystem::create_directories(path_.parent_path());
    }
    topology_ = common::toVtkUnstructuredGrid(*mesh_);
  }

  /** @copydoc ISeries::append */
  [[nodiscard]] OutputRecord append(const Request& request) override {
    if (finished_) {
      throw std::logic_error("cannot append to a finished output series");
    }
    if (request.mesh != nullptr && request.mesh != mesh_) {
      throw std::invalid_argument("snapshot belongs to a different mesh");
    }
    if (!request.path.empty() && normalizedOutputPath(request.path) != path_) {
      throw std::invalid_argument("snapshot belongs to a different path");
    }
    if (!std::isfinite(request.stamp.time)) {
      throw std::invalid_argument("output time must be finite");
    }
    if (!snapshots_.empty() &&
        (request.stamp.step <= snapshots_.back().stamp.step ||
         request.stamp.time <= snapshots_.back().stamp.time)) {
      throw std::invalid_argument(
          "output snapshots must have increasing step and time");
    }

    auto captured = captureSnapshot(*mesh_, request);
    if (snapshots_.empty()) {
      schema_ = std::move(captured.schema);
      if (captured.has_face_fields) {
        addFaceTopology(*topology_, *mesh_);
      }
      addSchemaMetadata(*topology_, schema_);
    } else if (captured.schema != schema_) {
      throw std::invalid_argument(
          "output field schema changed within a time series");
    }
    snapshots_.push_back(std::move(captured.snapshot));
    return {.path = path_, .stamp = snapshots_.back().stamp};
  }

  /** @copydoc ISeries::finish */
  [[nodiscard]] OutputRecord finish() override {
    if (finished_) {
      throw std::logic_error("output series has already been finished");
    }
    if (snapshots_.empty()) {
      throw std::logic_error("cannot finish an empty output series");
    }

    vtkNew<TemporalUnstructuredGridSource> source;
    source->setSeries(topology_, &snapshots_);
    vtkNew<vtkHDFWriter> writer;
    writer->SetFileName(path_.string().c_str());
    writer->SetOverwrite(overwrite_);
    writer->SetWriteAllTimeSteps(true);
    writer->SetUseExternalTimeSteps(false);
    writer->SetInputConnection(source->GetOutputPort());
    const auto success = writer->Write();

    if (success == 0 || writer->GetErrorCode() != vtkErrorCode::NoError ||
        !std::filesystem::is_regular_file(path_)) {
      throw std::runtime_error("vtkHDFWriter failed to create output series");
    }
    finished_ = true;
    return {.path = path_, .stamp = snapshots_.back().stamp};
  }

 private:
  const mesh::IMesh* mesh_{};
  std::filesystem::path path_;
  bool overwrite_{};
  bool finished_{};
  vtkSmartPointer<vtkUnstructuredGrid> topology_;
  std::vector<MetadataRow> schema_;
  std::vector<TemporalSnapshot> snapshots_;
};

}  // namespace

std::unique_ptr<ISeries> VtkHdfWriter::openSeries(
    const SeriesRequest& request) const {
  return std::make_unique<VtkHdfOutputSeries>(request);
}

}  // namespace pemu::output::dump
