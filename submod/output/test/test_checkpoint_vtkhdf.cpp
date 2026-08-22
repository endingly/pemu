#include <gtest/gtest.h>

#include <pemu/mesh/moab_mesh.hpp>
#include <pemu/output/checkpoint/vtkhdf.hpp>
#include <pemu/output/dump/vtkhdf_writer.hpp>
#include <pemu/unit/quantity_metadata.hpp>

#include <llnl-units/units.hpp>

#include <array>
#include <chrono>
#include <filesystem>
#include <string>
#include <system_error>
#include <type_traits>

namespace pemu::output::checkpoint::test {

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
    const auto nonce =
        std::chrono::steady_clock::now().time_since_epoch().count();
    path_ = std::filesystem::temp_directory_path() /
            ("pemu-checkpoint-test-" + std::to_string(nonce));
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
      .name = "electric potential",
      .physical_quantity =
          unit::PhysicalQuantityMetadata{unit::QuantityKind::electric_potential,
                                         units::precise::V},
  };
}

field::FieldMetadata electricFieldMetadata() {
  return {
      .name = "normal electric field",
      .physical_quantity =
          unit::PhysicalQuantityMetadata{
              unit::QuantityKind::normal_electric_field_strength,
              units::precise::V / units::precise::cm},
  };
}

static_assert(std::has_virtual_destructor_v<IWriter>);
static_assert(std::has_virtual_destructor_v<IReader>);

}  // namespace

TEST(VtkHdfCheckpointTest, RoundTripsManifestCellAndFaceState) {
  const mesh::MoabMesh mesh(testMeshPath().string());
  field::CellField<double> potential(mesh, potentialMetadata());
  potential[0] = 10.25;
  potential[1] = -3.5;
  field::FaceField<double> electric_field(mesh, electricFieldMetadata());
  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {
    electric_field[face] = 100.0 + static_cast<double>(face);
  }

  const std::array cell_sources{
      CellFieldSource{.field = &potential, .key = "potential"}};
  const std::array face_sources{FaceFieldSource{
      .field = &electric_field, .key = "normal_electric_field"}};
  double previous_time_step = 1.25e-9;
  const field::FieldMetadata previous_time_step_metadata{
      .name = "adaptive_previous_time_step"};
  const std::array scalar_sources{
      ScalarSource{.value = &previous_time_step,
                   .key = "adaptive_previous_time_step",
                   .metadata = previous_time_step_metadata}};
  TemporaryDirectory temporary;
  const auto path = temporary.path() / "state.vtkhdf";

  const VtkHdfWriter writer;
  const auto written = writer.write({
      .mesh = &mesh,
      .path = path,
      .stamp = {.step = 17, .time = 4.5e-8},
      .cell_fields = cell_sources,
      .face_fields = face_sources,
      .scalars = scalar_sources,
  });
  EXPECT_EQ(written.path, std::filesystem::absolute(path));

  const VtkHdfReader reader;
  const auto manifest = reader.inspect(path);
  EXPECT_EQ(manifest.format_version, vtkhdf_format_version);
  EXPECT_EQ(manifest.stamp.step, 17u);
  EXPECT_DOUBLE_EQ(manifest.stamp.time, 4.5e-8);
  EXPECT_EQ(manifest.num_cells, mesh.numCells());
  EXPECT_EQ(manifest.num_faces, mesh.numFaces());
  ASSERT_EQ(manifest.fields.size(), 3u);
  EXPECT_EQ(manifest.fields[0].association, FieldAssociation::cell);
  EXPECT_EQ(manifest.fields[0].key, "potential");
  EXPECT_EQ(manifest.fields[0].metadata.name, "electric potential");
  ASSERT_TRUE(manifest.fields[0].metadata.physical_quantity.has_value());
  EXPECT_EQ(manifest.fields[0].metadata.physical_quantity->kind(),
            unit::QuantityKind::electric_potential);
  EXPECT_EQ(manifest.fields[0].metadata.physical_quantity->unit(),
            units::precise::V);
  EXPECT_EQ(manifest.fields[2].association, FieldAssociation::scalar);
  EXPECT_EQ(manifest.fields[2].key, "adaptive_previous_time_step");
  EXPECT_EQ(manifest.fields[2].value_count, 1u);

  potential.fill(0.0);
  electric_field.fill(0.0);
  previous_time_step = 0.0;
  const std::array cell_targets{
      CellFieldTarget{.field = &potential, .key = "potential"}};
  const std::array face_targets{FaceFieldTarget{
      .field = &electric_field, .key = "normal_electric_field"}};
  const std::array scalar_targets{
      ScalarTarget{.value = &previous_time_step,
                   .key = "adaptive_previous_time_step",
                   .metadata = previous_time_step_metadata}};
  const auto restored = reader.restore(path, {.mesh = &mesh,
                                              .cell_fields = cell_targets,
                                              .face_fields = face_targets,
                                              .scalars = scalar_targets});

  EXPECT_EQ(restored.stamp.step, 17u);
  EXPECT_DOUBLE_EQ(restored.stamp.time, 4.5e-8);
  EXPECT_DOUBLE_EQ(potential[0], 10.25);
  EXPECT_DOUBLE_EQ(potential[1], -3.5);
  EXPECT_DOUBLE_EQ(previous_time_step, 1.25e-9);
  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {
    EXPECT_DOUBLE_EQ(electric_field[face], 100.0 + static_cast<double>(face));
  }
}

TEST(VtkHdfCheckpointTest, AtomicallyReplacesExistingCheckpoint) {
  const mesh::MoabMesh mesh(testMeshPath().string());
  field::CellField<double> potential(mesh, 1.0, potentialMetadata());
  const std::array sources{
      CellFieldSource{.field = &potential, .key = "potential"}};
  TemporaryDirectory temporary;
  const auto path = temporary.path() / "rolling.vtkhdf";
  const VtkHdfWriter writer;

  (void)writer.write(
      {.mesh = &mesh, .path = path, .cell_fields = sources, .overwrite = true});
  potential.fill(9.0);
  (void)writer.write(
      {.mesh = &mesh, .path = path, .cell_fields = sources, .overwrite = true});

  potential.fill(0.0);
  const std::array targets{
      CellFieldTarget{.field = &potential, .key = "potential"}};
  const VtkHdfReader reader;
  (void)reader.restore(path, {.mesh = &mesh, .cell_fields = targets});
  EXPECT_DOUBLE_EQ(potential[0], 9.0);
  EXPECT_DOUBLE_EQ(potential[1], 9.0);

  std::size_t artifact_count{};
  for (const auto& entry :
       std::filesystem::directory_iterator(temporary.path())) {
    ++artifact_count;
    EXPECT_EQ(entry.path(), path);
  }
  EXPECT_EQ(artifact_count, 1u);
}

TEST(VtkHdfCheckpointTest, ValidatesAllTargetsBeforeChangingAnyField) {
  const mesh::MoabMesh mesh(testMeshPath().string());
  field::CellField<double> first(mesh, 1.0, potentialMetadata());
  field::CellField<double> second(mesh, 2.0, potentialMetadata());
  const std::array sources{
      CellFieldSource{.field = &first, .key = "first"},
      CellFieldSource{.field = &second, .key = "second"},
  };
  TemporaryDirectory temporary;
  const auto path = temporary.path() / "atomic-restore.hdf";
  const VtkHdfWriter writer;
  (void)writer.write({.mesh = &mesh, .path = path, .cell_fields = sources});

  first.fill(9.0);
  second.fill(8.0);
  second.setMetadata({.name = "wrong metadata"});
  const std::array targets{
      CellFieldTarget{.field = &first, .key = "first"},
      CellFieldTarget{.field = &second, .key = "second"},
  };

  const VtkHdfReader reader;
  EXPECT_THROW(
      auto _ = reader.restore(path, {.mesh = &mesh, .cell_fields = targets}),
      std::invalid_argument);
  EXPECT_DOUBLE_EQ(first[0], 9.0);
  EXPECT_DOUBLE_EQ(second[0], 8.0);
}

TEST(VtkHdfCheckpointTest, RejectsKeysThatVtkHdfWouldRewrite) {
  const mesh::MoabMesh mesh(testMeshPath().string());
  field::CellField<double> potential(mesh, 1.0, potentialMetadata());
  const std::array sources{
      CellFieldSource{.field = &potential, .key = "potential/raw"}};
  TemporaryDirectory temporary;
  const VtkHdfWriter writer;

  EXPECT_THROW(auto _ = writer.write({.mesh = &mesh,
                                      .path = temporary.path() / "bad.vtkhdf",
                                      .cell_fields = sources}),
               std::invalid_argument);
}

TEST(VtkHdfCheckpointTest, RejectsVisualizationDumpWithoutCheckpointManifest) {
  const mesh::MoabMesh mesh(testMeshPath().string());
  field::CellField<double> potential(mesh, 1.0, potentialMetadata());
  const std::array fields{
      std::cref(static_cast<const field::CellField<double>&>(potential))};
  TemporaryDirectory temporary;
  const auto path = temporary.path() / "visualization.vtkhdf";
  const dump::VtkHdfWriter dump_writer;
  (void)dump_writer.write({.mesh = &mesh, .path = path, .cell_fields = fields});

  const VtkHdfReader checkpoint_reader;
  EXPECT_THROW(auto _ = checkpoint_reader.inspect(path), std::runtime_error);
}

}  // namespace pemu::output::checkpoint::test
