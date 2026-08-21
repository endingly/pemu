#include <pemu/output/vtkhdf_writer.hpp>

#include <pemu/output/mesh_adapter.hpp>
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
#include <vtkNew.h>
#include <vtkStringArray.h>
#include <vtkTypeUInt64Array.h>
#include <vtkUnstructuredGrid.h>

#include <cmath>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>

namespace pemu::output {

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

std::string unitName(const unit::PhysicalQuantityMetadata& metadata) {
  return metadata.unit() == units::precise::one
             ? std::string{"1"}
             : units::to_string(metadata.unit());
}

bool isReservedFieldDataName(std::string_view name) {
  return name == kFieldAssociationArrayName || name == kFieldNameArrayName ||
         name == kFieldQuantityKindArrayName || name == kFieldUnitArrayName ||
         name == kFaceOwnerArrayName || name == kFaceNeighborArrayName ||
         name == kFaceCenterArrayName || name == kFaceAreaArrayName ||
         name == kFaceBoundaryArrayName || name == kStepArrayName ||
         name == kTimeArrayName;
}

void appendMetadata(vtkStringArray& associations, vtkStringArray& names,
                    vtkStringArray& quantity_kinds, vtkStringArray& units,
                    std::string_view association,
                    const field::FieldMetadata& field_metadata) {
  associations.InsertNextValue(std::string{association});
  names.InsertNextValue(field_metadata.name);

  if (field_metadata.physical_quantity.has_value()) {
    quantity_kinds.InsertNextValue(std::string{
        pemu::to_string(field_metadata.physical_quantity->kind())});
    units.InsertNextValue(unitName(*field_metadata.physical_quantity));
  } else {
    quantity_kinds.InsertNextValue("");
    units.InsertNextValue("");
  }
}

void validateField(const mesh::IMesh& mesh,
                   const field::CellField<double>& values) {
  if (&values.mesh() != &mesh) {
    throw std::invalid_argument("cell field belongs to a different mesh");
  }
  if (values.metadata().name.empty()) {
    throw std::invalid_argument("cell field name must not be empty");
  }
}

void validateField(const mesh::IMesh& mesh,
                   const field::FaceField<double>& values) {
  if (&values.mesh() != &mesh) {
    throw std::invalid_argument("face field belongs to a different mesh");
  }
  if (values.metadata().name.empty()) {
    throw std::invalid_argument("face field name must not be empty");
  }
}

vtkSmartPointer<vtkDoubleArray> makeDoubleArray(
    std::string_view name, std::span<const double> values) {
  vtkNew<vtkDoubleArray> array;
  array->SetName(std::string{name}.c_str());
  array->SetNumberOfValues(static_cast<vtkIdType>(values.size()));
  for (std::size_t i = 0; i < values.size(); ++i) {
    array->SetValue(static_cast<vtkIdType>(i), values[i]);
  }
  return array;
}

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
    neighbors->SetValue(
        index, mesh.isBoundary(face)
                   ? vtkIdType{-1}
                   : static_cast<vtkIdType>(mesh.neighbor(face)));
    centers->SetTuple3(index, center.x, center.y, center.z);
    areas->SetValue(index, mesh.faceArea(face));
    boundary_ids->SetValue(
        index, mesh.boundaryId(face) == mesh::invalid_boundary
                   ? vtkIdType{-1}
                   : static_cast<vtkIdType>(mesh.boundaryId(face)));
  }

  grid.GetFieldData()->AddArray(owners);
  grid.GetFieldData()->AddArray(neighbors);
  grid.GetFieldData()->AddArray(centers);
  grid.GetFieldData()->AddArray(areas);
  grid.GetFieldData()->AddArray(boundary_ids);
}

vtkSmartPointer<vtkDoubleArray> makeCellCenteredFaceField(
    const mesh::IMesh& mesh, const field::FaceField<double>& face_field,
    std::string_view name) {
  vtkNew<vtkDoubleArray> result;
  result->SetName(std::string{name}.c_str());
  result->SetNumberOfValues(static_cast<vtkIdType>(mesh.numCells()));

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

    result->SetValue(static_cast<vtkIdType>(cell),
                     weighted_sum / area_sum);
  }

  return result;
}

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

}  // namespace

OutputRecord VtkHdfWriter::write(const FieldDumpRequest& request) const {
  if (request.mesh == nullptr) {
    throw std::invalid_argument("field dump requires a mesh");
  }
  if (!std::isfinite(request.stamp.time)) {
    throw std::invalid_argument("output time must be finite");
  }

  const auto path = normalizedOutputPath(request.path);
  if (std::filesystem::exists(path) && !request.overwrite) {
    throw std::filesystem::filesystem_error(
        "refusing to overwrite VTKHDF output", path,
        std::make_error_code(std::errc::file_exists));
  }
  if (!path.parent_path().empty()) {
    std::filesystem::create_directories(path.parent_path());
  }

  const mesh::IMesh& mesh = *request.mesh;
  auto grid = toVtkUnstructuredGrid(mesh);
  grid->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(),
                              request.stamp.time);

  vtkNew<vtkTypeUInt64Array> step;
  step->SetName(kStepArrayName.data());
  step->InsertNextValue(request.stamp.step);
  grid->GetFieldData()->AddArray(step);

  vtkNew<vtkDoubleArray> time;
  time->SetName(kTimeArrayName.data());
  time->InsertNextValue(request.stamp.time);
  grid->GetFieldData()->AddArray(time);

  vtkNew<vtkStringArray> field_associations;
  field_associations->SetName(kFieldAssociationArrayName.data());
  vtkNew<vtkStringArray> field_names;
  field_names->SetName(kFieldNameArrayName.data());
  vtkNew<vtkStringArray> field_quantity_kinds;
  field_quantity_kinds->SetName(kFieldQuantityKindArrayName.data());
  vtkNew<vtkStringArray> field_units;
  field_units->SetName(kFieldUnitArrayName.data());

  std::unordered_set<std::string> cell_names;
  for (const auto& field_reference : request.cell_fields) {
    const auto& values = field_reference.get();
    validateField(mesh, values);
    if (!cell_names.insert(values.metadata().name).second) {
      throw std::invalid_argument("duplicate cell field name");
    }

    grid->GetCellData()->AddArray(
        makeDoubleArray(values.metadata().name, values.span()));
    appendMetadata(*field_associations, *field_names, *field_quantity_kinds,
                   *field_units, "cell", values.metadata());
  }

  if (!request.face_fields.empty()) {
    addFaceTopology(*grid, mesh);
  }

  std::unordered_set<std::string> face_names;
  for (const auto& selection : request.face_fields) {
    if (selection.field == nullptr) {
      throw std::invalid_argument("face field selection must not be null");
    }
    const auto& values = *selection.field;
    validateField(mesh, values);
    if (isReservedFieldDataName(values.metadata().name)) {
      throw std::invalid_argument(
          "face field name is reserved by VTKHDF output");
    }
    if (!face_names.insert(values.metadata().name).second) {
      throw std::invalid_argument("duplicate face field name");
    }

    // VTK has no face-data association for an unstructured grid. Keep the
    // raw face-order values in FieldData, alongside explicit face topology.
    grid->GetFieldData()->AddArray(
        makeDoubleArray(values.metadata().name, values.span()));
    appendMetadata(*field_associations, *field_names, *field_quantity_kinds,
                   *field_units, "face", values.metadata());

    if (selection.include_cell_centered_visualization) {
      const std::string centered_name =
          selection.cell_centered_name.empty()
          ? values.metadata().name + "_cell_centered"
          : selection.cell_centered_name;
      if (!cell_names.insert(centered_name).second) {
        throw std::invalid_argument("duplicate cell visualization field name");
      }

      grid->GetCellData()->AddArray(
          makeCellCenteredFaceField(mesh, values, centered_name));
      auto centered_metadata = values.metadata();
      centered_metadata.name = centered_name;
      appendMetadata(*field_associations, *field_names, *field_quantity_kinds,
                     *field_units, "cell_visualization", centered_metadata);
    }
  }

  grid->GetFieldData()->AddArray(field_associations);
  grid->GetFieldData()->AddArray(field_names);
  grid->GetFieldData()->AddArray(field_quantity_kinds);
  grid->GetFieldData()->AddArray(field_units);

  vtkNew<vtkHDFWriter> writer;
  writer->SetFileName(path.string().c_str());
  writer->SetOverwrite(request.overwrite);
  writer->SetInputData(grid);
  writer->Write();

  if (writer->GetErrorCode() != vtkErrorCode::NoError ||
      !std::filesystem::is_regular_file(path)) {
    throw std::runtime_error("vtkHDFWriter failed to create output");
  }

  return {.path = path, .stamp = request.stamp};
}

}  // namespace pemu::output
