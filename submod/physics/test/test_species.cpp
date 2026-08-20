#include <gtest/gtest.h>

#include <pemu/mesh/moab_mesh.hpp>
#include <pemu/physics/species.hpp>

#include <filesystem>

namespace pemu::physics::test {

namespace {

std::filesystem::path testMeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
}

class SpeciesFieldsTest : public ::testing::Test {
 protected:
  SpeciesFieldsTest() : mesh_(testMeshPath().string()) {}
  mesh::MoabMesh mesh_;
};

};  // namespace

TEST(SpeciesSetTest, AssignsDenseStableIds) {
  physics::SpeciesSet species;
  const auto electron = species.add({.name = "e", .charge = -1.0});
  const auto ion = species.add({.name = "Ar+", .charge = 1.0});
  EXPECT_EQ(electron.value, 0u);
  EXPECT_EQ(ion.value, 1u);
  EXPECT_EQ(species.size(), 2u);
}

TEST_F(SpeciesFieldsTest, StoresIndependentFieldsPerSpecies) {
  physics::SpeciesSet species;
  const auto electron = species.add({.name = "e", .charge = -1.0});
  const auto ion = species.add({.name = "ion", .charge = 1.0});
  physics::SpeciesCellFields density(mesh_, species.size(), 0.0);
  density[electron][0] = 2.0;
  density[ion][0] = 5.0;
  EXPECT_DOUBLE_EQ(density[electron][0], 2.0);
  EXPECT_DOUBLE_EQ(density[ion][0], 5.0);
}

};  // namespace pemu::physics::test
