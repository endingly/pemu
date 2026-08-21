#include <gtest/gtest.h>

#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>
#include <pemu/mesh/moab_mesh.hpp>
#include <pemu/output/i_checkpoint_writer.hpp>
#include <pemu/output/mesh_adapter.hpp>
#include <pemu/output/vtkhdf_writer.hpp>
#include <pemu/unit/quantity_metadata.hpp>

#include <llnl-units/units.hpp>

#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkDataArray.h>
#include <vtkFieldData.h>
#include <vtkHDFReader.h>
#include <vtkNew.h>
#include <vtkStringArray.h>
#include <vtkUnstructuredGrid.h>

#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <string>
#include <system_error>
#include <type_traits>

namespace pemu::output::test {

namespace {

std::filesystem::path testMeshPath() {
#ifdef PEMU_MESH_TEST_DATA_DIR
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
#else
  return std::filesystem::path{"data/meshfiles/two_quads.msh"};
#endif
}

class TemporaryDirectory {
 public:
  TemporaryDirectory() {
    const auto nonce = std::chrono::steady_clock::now()
                           .time_since_epoch()
                           .count();
    path_ = std::filesystem::temp_directory_path() /
            ("pemu-vtkhdf-test-" + std::to_string(nonce));
    std::filesystem::create_directories(path_);
  }

  ~TemporaryDirectory() {
    std::error_code error;
    std::filesystem::remove_all(path_, error);
  }

  const std::filesystem::path& path() const noexcept { return path_; }

 private:
  std::filesystem::path path_;
};

field::FieldMetadata potentialMetadata() {
  return {
      .name = "potential",
      .physical_quantity = unit::PhysicalQuantityMetadata{
          unit::QuantityKind::electric_potential, units::precise::V},
  };
}

field::FieldMetadata faceFluxMetadata() {
  return {
      .name = "face_flux",
      .physical_quantity = unit::PhysicalQuantityMetadata{
          unit::QuantityKind::normal_electric_field_strength,
          units::precise::V / units::precise::cm},
  };
}

vtkSmartPointer<vtkUnstructuredGrid> readGrid(
    const std::filesystem::path& path) {
  vtkNew<vtkHDFReader> reader;
  reader->SetFileName(path.string().c_str());
  reader->Update();

  auto* output =
      vtkUnstructuredGrid::SafeDownCast(reader->GetOutputDataObject(0));
  if (output == nullptr) {
    return {};
  }

  vtkSmartPointer<vtkUnstructuredGrid> grid =
      vtkSmartPointer<vtkUnstructuredGrid>::New();
  grid->DeepCopy(output);
  return grid;
}

bool hasMetadataRow(vtkFieldData& field_data, std::string_view association,
                    std::string_view name, std::string_view quantity_kind,
                    std::string_view unit_name) {
  auto* associations = vtkStringArray::SafeDownCast(
      field_data.GetAbstractArray("pemu_field_association"));
  auto* names = vtkStringArray::SafeDownCast(
      field_data.GetAbstractArray("pemu_field_name"));
  auto* quantity_kinds = vtkStringArray::SafeDownCast(
      field_data.GetAbstractArray("pemu_field_quantity_kind"));
  auto* units = vtkStringArray::SafeDownCast(
      field_data.GetAbstractArray("pemu_field_unit"));
  if (associations == nullptr || names == nullptr ||
      quantity_kinds == nullptr || units == nullptr ||
      associations->GetNumberOfValues() != names->GetNumberOfValues() ||
      associations->GetNumberOfValues() !=
          quantity_kinds->GetNumberOfValues() ||
      associations->GetNumberOfValues() != units->GetNumberOfValues()) {
    return false;
  }

  for (vtkIdType i = 0; i < associations->GetNumberOfValues(); ++i) {
    if (associations->GetValue(i) == association &&
        names->GetValue(i) == name &&
        quantity_kinds->GetValue(i) == quantity_kind &&
        units->GetValue(i) == unit_name) {
      return true;
    }
  }
  return false;
}

struct CapturedEvent {
  trace::DiagDomain domain{trace::DiagDomain::trace};
  std::string category;
  std::string name;
  std::string path;
  std::uint64_t step{};
  double time{};
  std::size_t attribute_count{};
};

class CapturingSink {
 public:
  void operator()(const trace::TraceEvent& event) noexcept {
    captured.domain = event.domain;
    captured.category = event.category;
    captured.name = event.name;
    captured.attribute_count = event.attributes.size();

    for (const auto& attribute : event.attributes) {
      if (attribute.name == "path") {
        captured.path =
            std::string{std::get<std::string_view>(attribute.value)};
      } else if (attribute.name == "step") {
        captured.step = std::get<std::uint64_t>(attribute.value);
      } else if (attribute.name == "time") {
        captured.time = std::get<double>(attribute.value);
      }
    }
  }

  CapturedEvent captured;
};

static_assert(trace::TraceSink<CapturingSink>);
static_assert(
    std::has_virtual_destructor_v<ICheckpointWriter<CapturedEvent>>);

}  // namespace

TEST(MeshAdapterTest, ConvertsIMeshToPolygonalUnstructuredGrid) {
  const mesh::MoabMesh mesh(testMeshPath().string());
  const auto grid = toVtkUnstructuredGrid(mesh);

  ASSERT_NE(grid, nullptr);
  EXPECT_EQ(grid->GetNumberOfPoints(), 6);
  EXPECT_EQ(grid->GetNumberOfCells(), 2);
  for (vtkIdType cell = 0; cell < grid->GetNumberOfCells(); ++cell) {
    EXPECT_EQ(grid->GetCellType(cell), VTK_POLYGON);
    EXPECT_EQ(grid->GetCell(cell)->GetNumberOfPoints(), 4);
  }
}

TEST(ParaViewReadabilityTest,
     OfficialVtkHdfReaderRoundTripsFieldsMetadataAndTime) {
  const mesh::MoabMesh mesh(testMeshPath().string());
  field::CellField<double> potential(mesh, potentialMetadata());
  potential[0] = 10.0;
  potential[1] = 20.0;

  field::FaceField<double> face_flux(mesh, faceFluxMetadata());
  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {
    face_flux[face] = 100.0 + static_cast<double>(face);
  }

  const std::array cell_fields{
      std::cref(static_cast<const field::CellField<double>&>(potential))};
  const std::array face_fields{FaceFieldSelection{
      .field = &face_flux,
      .include_cell_centered_visualization = true,
      .cell_centered_name = "face_flux_viz",
  }};
  TemporaryDirectory temporary;
  const auto path = temporary.path() / "step-00000042.vtkhdf";

  const VtkHdfWriter writer;
  const auto record = writer.write({
      .mesh = &mesh,
      .path = path,
      .stamp = {.step = 42, .time = 1.25e-6},
      .cell_fields = cell_fields,
      .face_fields = face_fields,
  });

  EXPECT_EQ(record.path, std::filesystem::absolute(path));
  const auto grid = readGrid(record.path);
  ASSERT_NE(grid, nullptr);
  EXPECT_EQ(grid->GetNumberOfPoints(), 6);
  EXPECT_EQ(grid->GetNumberOfCells(), 2);

  auto* potential_values = grid->GetCellData()->GetArray("potential");
  ASSERT_NE(potential_values, nullptr);
  EXPECT_DOUBLE_EQ(potential_values->GetTuple1(0), 10.0);
  EXPECT_DOUBLE_EQ(potential_values->GetTuple1(1), 20.0);

  // Raw face values remain FieldData and are never presented as CellData.
  EXPECT_EQ(grid->GetCellData()->GetArray("face_flux"), nullptr);
  auto* raw_face_values = grid->GetFieldData()->GetArray("face_flux");
  ASSERT_NE(raw_face_values, nullptr);
  ASSERT_EQ(raw_face_values->GetNumberOfTuples(),
            static_cast<vtkIdType>(mesh.numFaces()));
  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {
    EXPECT_DOUBLE_EQ(raw_face_values->GetTuple1(face), face_flux[face]);
  }

  auto* owners = grid->GetFieldData()->GetArray("pemu_face_owner");
  auto* neighbors = grid->GetFieldData()->GetArray("pemu_face_neighbor");
  ASSERT_NE(owners, nullptr);
  ASSERT_NE(neighbors, nullptr);
  ASSERT_EQ(owners->GetNumberOfTuples(),
            static_cast<vtkIdType>(mesh.numFaces()));
  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {
    EXPECT_EQ(static_cast<vtkIdType>(owners->GetTuple1(face)),
              mesh.owner(face));
    EXPECT_EQ(static_cast<vtkIdType>(neighbors->GetTuple1(face)),
              mesh.isBoundary(face)
                  ? -1
                  : static_cast<vtkIdType>(mesh.neighbor(face)));
  }

  auto* centered = grid->GetCellData()->GetArray("face_flux_viz");
  ASSERT_NE(centered, nullptr);
  for (mesh::CellId cell = 0; cell < mesh.numCells(); ++cell) {
    double expected = 0.0;
    double total_area = 0.0;
    for (const mesh::FaceId face : mesh.cellFaces(cell)) {
      expected += face_flux[face] * mesh.faceArea(face);
      total_area += mesh.faceArea(face);
    }
    EXPECT_DOUBLE_EQ(centered->GetTuple1(cell), expected / total_area);
  }

  auto* step = grid->GetFieldData()->GetArray("pemu_step");
  auto* time = grid->GetFieldData()->GetArray("pemu_time");
  ASSERT_NE(step, nullptr);
  ASSERT_NE(time, nullptr);
  EXPECT_EQ(static_cast<std::uint64_t>(step->GetTuple1(0)), 42u);
  EXPECT_DOUBLE_EQ(time->GetTuple1(0), 1.25e-6);

  EXPECT_TRUE(hasMetadataRow(*grid->GetFieldData(), "cell", "potential",
                             "electric_potential", "V"));
  EXPECT_TRUE(hasMetadataRow(*grid->GetFieldData(), "face", "face_flux",
                             "normal_electric_field_strength", "V/cm"));
  EXPECT_TRUE(hasMetadataRow(*grid->GetFieldData(), "cell_visualization",
                             "face_flux_viz",
                             "normal_electric_field_strength", "V/cm"));
}

TEST(OutputTraceTest, EmitsOnlyLightweightCompletionContextAfterWrite) {
  const mesh::MoabMesh mesh(testMeshPath().string());
  field::CellField<double> potential(mesh, 3.0, potentialMetadata());
  const std::array cell_fields{
      std::cref(static_cast<const field::CellField<double>&>(potential))};
  TemporaryDirectory temporary;
  CapturingSink sink;

  const VtkHdfWriter writer;
  const auto record = writer.writeAndTrace(
      {
          .mesh = &mesh,
          .path = temporary.path() / "trace.vtkhdf",
          .stamp = {.step = 7, .time = 0.5},
          .cell_fields = cell_fields,
      },
      sink);

  EXPECT_EQ(sink.captured.domain, trace::DiagDomain::output);
  EXPECT_TRUE(sink.captured.category.empty());
  EXPECT_EQ(sink.captured.name, "completed");
  EXPECT_EQ(sink.captured.attribute_count, 3u);
  EXPECT_EQ(sink.captured.path, record.path.string());
  EXPECT_EQ(sink.captured.step, 7u);
  EXPECT_DOUBLE_EQ(sink.captured.time, 0.5);
}

TEST(VtkHdfWriterTest, RefusesOverwriteUnlessExplicitlyEnabled) {
  const mesh::MoabMesh mesh(testMeshPath().string());
  field::CellField<double> potential(mesh, 3.0, potentialMetadata());
  const std::array cell_fields{
      std::cref(static_cast<const field::CellField<double>&>(potential))};
  TemporaryDirectory temporary;
  const auto path = temporary.path() / "existing.vtkhdf";
  const VtkHdfWriter writer;
  CapturingSink sink;

  (void)writer.write(
      {.mesh = &mesh, .path = path, .cell_fields = cell_fields});
  EXPECT_THROW(
      writer.writeAndTrace(
          {.mesh = &mesh, .path = path, .cell_fields = cell_fields}, sink),
      std::filesystem::filesystem_error);
  EXPECT_EQ(sink.captured.attribute_count, 0u);
  EXPECT_NO_THROW(writer.write({.mesh = &mesh,
                                .path = path,
                                .cell_fields = cell_fields,
                                .overwrite = true}));
}

}  // namespace pemu::output::test
