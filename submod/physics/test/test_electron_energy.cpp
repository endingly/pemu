#include <gtest/gtest.h>

#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>
#include <pemu/mesh/moab_mesh.hpp>
#include <pemu/physics/electron_energy.hpp>

#include <filesystem>
#include <limits>

namespace pemu::physics::test {
namespace {

/** @brief Returns the orthogonal two-cell fixture used by field-power tests. */
[[nodiscard]] std::filesystem::path testMeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
}

class ElectronFieldPowerTest : public ::testing::Test {
 protected:
  /** @brief Loads a pair of unit square finite-volume cells. */
  ElectronFieldPowerTest() : mesh_(testMeshPath().string()) {}

  mesh::MoabMesh mesh_;
};

}  // namespace

TEST(ElectronTemperatureTest, ConvertsMaxwellianMeanEnergyToTeInEv) {
  EXPECT_DOUBLE_EQ(electronTemperatureEv(0.0), 0.0);
  EXPECT_DOUBLE_EQ(electronTemperatureEv(3.0), 2.0);
  EXPECT_THROW(static_cast<void>(electronTemperatureEv(-1.0)),
               std::invalid_argument);
  EXPECT_THROW(static_cast<void>(
                   electronTemperatureEv(std::numeric_limits<double>::infinity())),
               std::invalid_argument);
}

TEST_F(ElectronFieldPowerTest, ConvertsMeanEnergyFieldAtomically) {
  field::CellField<double> mean_energy(mesh_, 0.0);
  field::CellField<double> temperature(mesh_, 91.0);
  mean_energy[mesh::CellId{0}] = 1.5;
  mean_energy[mesh::CellId{1}] = 6.0;

  computeElectronTemperatureEv(mean_energy, temperature);

  EXPECT_DOUBLE_EQ(temperature[mesh::CellId{0}], 1.0);
  EXPECT_DOUBLE_EQ(temperature[mesh::CellId{1}], 4.0);

  mean_energy[mesh::CellId{1}] = -1.0;
  temperature.fill(91.0);
  EXPECT_THROW(computeElectronTemperatureEv(mean_energy, temperature),
               std::invalid_argument);
  EXPECT_DOUBLE_EQ(temperature[mesh::CellId{0}], 91.0);
  EXPECT_DOUBLE_EQ(temperature[mesh::CellId{1}], 91.0);
}

TEST_F(ElectronFieldPowerTest,
       ReconstructsPositiveHeatingForElectronFluxAgainstField) {
  field::FaceField<double> electron_flux_normal(mesh_, 0.0);
  field::FaceField<double> electric_field_normal(mesh_, 0.0);
  field::CellField<double> field_power_density(mesh_, 0.0);

  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {
    const auto normal = mesh_.faceNormal(face);
    electron_flux_normal[face] = -2.0 * normal.x;
    electric_field_normal[face] = 3.0 * normal.x;
  }

  computeElectronFieldPowerDensity(electron_flux_normal,
                                   electric_field_normal,
                                   field_power_density);

  for (const double value : field_power_density) {
    EXPECT_NEAR(value, 6.0, 1e-12);
  }
}

TEST_F(ElectronFieldPowerTest, RejectsFieldsFromDifferentMeshes) {
  mesh::MoabMesh other_mesh(testMeshPath().string());
  field::FaceField<double> electron_flux_normal(mesh_, 0.0);
  field::FaceField<double> electric_field_normal(other_mesh, 0.0);
  field::CellField<double> field_power_density(mesh_, 0.0);

  EXPECT_THROW(computeElectronFieldPowerDensity(
                   electron_flux_normal, electric_field_normal,
                   field_power_density),
               std::invalid_argument);
}

TEST_F(ElectronFieldPowerTest, RejectsNonFiniteFaceData) {
  field::FaceField<double> electron_flux_normal(mesh_, 0.0);
  field::FaceField<double> electric_field_normal(mesh_, 0.0);
  field::CellField<double> field_power_density(mesh_, 0.0);
  electron_flux_normal[mesh::FaceId{0}] =
      std::numeric_limits<double>::infinity();

  EXPECT_THROW(computeElectronFieldPowerDensity(
                   electron_flux_normal, electric_field_normal,
                   field_power_density),
               std::invalid_argument);
}

}  // namespace pemu::physics::test
