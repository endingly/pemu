#include <gtest/gtest.h>

#include <filesystem>
#include <pemu/boundary/boundary_condition.hpp>
#include <pemu/boundary/boundary_condition_set.hpp>
#include <pemu/mesh/moab_mesh.hpp>
#include <variant>

namespace pemu::boundary::test {

TEST(BoundaryConditionTest, StoresDirichlet) {
  BoundaryConditionSet bc;

  bc.setDirichlet(mesh::BoundaryId{1}, 3.5);

  ASSERT_TRUE(bc.contains(mesh::BoundaryId{1}));

  const auto& condition = bc.at(mesh::BoundaryId{1});

  ASSERT_TRUE(std::holds_alternative<Dirichlet>(condition));

  EXPECT_DOUBLE_EQ(std::get<Dirichlet>(condition).value, 3.5);
}

TEST(BoundaryConditionTest, StoresNeumann) {
  BoundaryConditionSet bc;

  bc.setNeumann(mesh::BoundaryId{2}, -1.25);

  const auto& condition = bc.at(mesh::BoundaryId{2});

  ASSERT_TRUE(std::holds_alternative<Neumann>(condition));

  EXPECT_DOUBLE_EQ(std::get<Neumann>(condition).value, -1.25);
}

TEST(BoundaryConditionTest, ReplacesExistingCondition) {
  BoundaryConditionSet bc;

  bc.setDirichlet(mesh::BoundaryId{1}, 1.0);

  bc.setNeumann(mesh::BoundaryId{1}, 2.0);

  const auto& condition = bc.at(mesh::BoundaryId{1});

  EXPECT_TRUE(std::holds_alternative<Neumann>(condition));

  EXPECT_DOUBLE_EQ(std::get<Neumann>(condition).value, 2.0);
}

TEST(BoundaryConditionTest, MissingBoundaryThrows) {
  BoundaryConditionSet bc;

  EXPECT_THROW(bc.at(mesh::BoundaryId{999}), std::out_of_range);
}

TEST(BoundaryConditionTest, ReportsSize) {
  BoundaryConditionSet bc;

  bc.setDirichlet(mesh::BoundaryId{1}, 0.0);

  bc.setNeumann(mesh::BoundaryId{2}, 0.0);

  EXPECT_EQ(bc.size(), 2u);
}

namespace {

std::filesystem::path testMeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
}

class BoundaryIntegrationTest : public ::testing::Test {
 protected:
  BoundaryIntegrationTest() : mesh_(testMeshPath().string()) {}

  mesh::MoabMesh mesh_;
};

}  // namespace

TEST_F(BoundaryIntegrationTest, EveryBoundaryFaceHasCondition) {
  BoundaryConditionSet bc;

  bc.setDirichlet(mesh::BoundaryId{1}, 0.0);

  bc.setDirichlet(mesh::BoundaryId{2}, 0.0);

  bc.setDirichlet(mesh::BoundaryId{3}, 0.0);

  bc.setDirichlet(mesh::BoundaryId{4}, 0.0);

  std::size_t boundary_count = 0;

  for (mesh::FaceId f = 0; f < mesh_.numFaces(); ++f) {

    if (!mesh_.isBoundary(f)) {
      continue;
    }

    ++boundary_count;

    const auto id = mesh_.boundaryId(f);

    ASSERT_NE(id, mesh::invalid_boundary);

    EXPECT_TRUE(bc.contains(id));
  }

  EXPECT_EQ(boundary_count, 6u);
}

TEST_F(BoundaryIntegrationTest, DetectsMissingBoundaryCondition) {
  BoundaryConditionSet bc;

  bc.setDirichlet(mesh::BoundaryId{1}, 0.0);

  bc.setDirichlet(mesh::BoundaryId{2}, 0.0);

  bc.setDirichlet(mesh::BoundaryId{3}, 0.0);

  // BoundaryId{4} deliberately missing.

  bool found_missing = false;

  for (mesh::FaceId f = 0; f < mesh_.numFaces(); ++f) {

    if (!mesh_.isBoundary(f)) {
      continue;
    }

    const auto id = mesh_.boundaryId(f);

    if (!bc.contains(id)) {
      found_missing = true;
      EXPECT_EQ(id, mesh::BoundaryId{4});
    }
  }

  EXPECT_TRUE(found_missing);
}

}  // namespace pemu::boundary::test