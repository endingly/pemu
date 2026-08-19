#include <gtest/gtest.h>
#include <filesystem>
#include <pemu/mesh/moab_mesh.hpp>
#include <pemu/physics/charge_density.hpp>

namespace pemu::physics::test {

namespace {

std::filesystem::path testMeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
}

class ChargeDensityTest : public ::testing::Test {
 protected:
  ChargeDensityTest() : mesh_(testMeshPath().string()) {}

  mesh::MoabMesh mesh_;
};

};  // namespace

TEST_F(ChargeDensityTest, EqualOppositeSpeciesAreNeutral) {
  field::CellField<double> electron_density(mesh_, 3.0);

  field::CellField<double> ion_density(mesh_, 3.0);

  field::CellField<double> rho(mesh_, 0.0);

  physics::addSpeciesChargeDensity(electron_density, -2.0, rho);

  physics::addSpeciesChargeDensity(ion_density, +2.0, rho);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(rho[cell], 0.0, 1e-12);
  }
}

};  // namespace pemu::physics::test