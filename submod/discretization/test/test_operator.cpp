#include <gtest/gtest.h>
#include <filesystem>
#include <pemu/boundary/boundary_condition_set.hpp>
#include <pemu/discretization/operators/bernoulli.hpp>
#include <pemu/discretization/operators/diffusion_flux.hpp>
#include <pemu/discretization/operators/divergence.hpp>
#include <pemu/discretization/operators/scharfetter_gummel_flux.hpp>
#include <pemu/discretization/operators/upwind_advection_flux.hpp>
#include <pemu/mesh/moab_mesh.hpp>

#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>

namespace pemu::discretization::test {

namespace {
mesh::FaceId findInternalFace(const mesh::IMesh& mesh) {
  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {

    if (!mesh.isBoundary(face)) {
      return face;
    }
  }
  throw std::runtime_error("internal face not found");
}

std::filesystem::path testMeshPath() {
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
}

class OperatorTest : public ::testing::Test {
 protected:
  OperatorTest() : mesh_(testMeshPath().string()) {}

  mesh::MoabMesh mesh_;
};

mesh::FaceId findBoundaryFace(const mesh::IMesh& mesh, mesh::BoundaryId id) {
  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {

    if (!mesh.isBoundary(face)) {
      continue;
    }

    if (mesh.boundaryId(face) == id) {
      return face;
    }
  }

  throw std::runtime_error("boundary face not found");
}

};  // namespace

TEST_F(OperatorTest, DivergenceConservesInternalFlux) {
  field::FaceField<double> flux(mesh_, 0.0);

  field::CellField<double> div(mesh_, 0.0);

  const auto internal_face = findInternalFace(mesh_);

  flux[internal_face] = 2.5;

  discretization::operators::divergence(flux, div);

  double integral = 0.0;

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    integral += div[cell] * mesh_.cellVolume(cell);
  }

  EXPECT_NEAR(integral, 0.0, 1e-12);
}

TEST_F(OperatorTest, DivergenceUsesOwnerNeighborOrientation) {
  field::FaceField<double> flux(mesh_, 0.0);

  field::CellField<double> div(mesh_, 0.0);

  const auto face = findInternalFace(mesh_);

  const auto owner = mesh_.owner(face);

  const auto neighbor = mesh_.neighbor(face);

  flux[face] = 2.0;

  discretization::operators::divergence(flux, div);

  EXPECT_NEAR(div[owner], +2.0, 1e-12);

  EXPECT_NEAR(div[neighbor], -2.0, 1e-12);
}

TEST_F(OperatorTest, DivergenceSatisfiesDiscreteGaussTheorem) {
  field::FaceField<double> flux(mesh_, 0.0);

  field::CellField<double> div(mesh_, 0.0);

  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {

    if (mesh_.isBoundary(face)) {
      flux[face] = 1.0;
    }
  }

  discretization::operators::divergence(flux, div);

  double volume_integral = 0.0;

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    volume_integral += div[cell] * mesh_.cellVolume(cell);
  }

  double boundary_integral = 0.0;

  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {

    if (!mesh_.isBoundary(face)) {
      continue;
    }

    boundary_integral += flux[face] * mesh_.faceArea(face);
  }

  EXPECT_NEAR(volume_integral, boundary_integral, 1e-12);

  EXPECT_NEAR(volume_integral, 6.0, 1e-12);
}

TEST_F(OperatorTest, DiffusionFluxIsExactForLinearField) {
  field::CellField<double> u(mesh_);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    const auto center = mesh_.cellCenter(cell);

    u[cell] = center.x;
  }

  boundary::BoundaryConditionSet bc;

  // u = x

  bc.setDirichlet(1, 0.0);  // left
  bc.setDirichlet(2, 2.0);  // right

  //
  // y-boundaries:
  //
  // grad(u) dot n = 0
  //
  bc.setNeumann(3, 0.0);
  bc.setNeumann(4, 0.0);

  field::FaceField<double> flux(mesh_, 0.0);

  discretization::operators::diffusionFlux(u, 2.0, bc, flux);

  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {

    const auto n = mesh_.faceNormal(face);

    const double expected = -2.0 * n.x;

    EXPECT_NEAR(flux[face], expected, 1e-12);
  }
}

TEST_F(OperatorTest, DivergenceOfLinearDiffusionFluxIsZero) {
  field::CellField<double> u(mesh_);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    const auto center = mesh_.cellCenter(cell);

    u[cell] = center.x;
  }

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(1, 0.0);
  bc.setDirichlet(2, 2.0);

  bc.setNeumann(3, 0.0);
  bc.setNeumann(4, 0.0);

  field::FaceField<double> flux(mesh_, 0.0);

  field::CellField<double> div(mesh_, 0.0);

  discretization::operators::diffusionFlux(u, 1.0, bc, flux);

  discretization::operators::divergence(flux, div);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(div[cell], 0.0, 1e-12);
  }
}

TEST_F(OperatorTest, DivergenceRejectsDifferentMeshes) {
  mesh::MoabMesh another_mesh(testMeshPath().string());

  field::FaceField<double> flux(mesh_, 0.0);

  field::CellField<double> div(another_mesh, 0.0);

  EXPECT_THROW(discretization::operators::divergence(flux, div),
               std::invalid_argument);
}

TEST_F(OperatorTest, UpwindFluxUsesOwnerForPositiveVelocity) {
  field::CellField<double> u(mesh_, 0.0);

  u[0] = 2.0;
  u[1] = 5.0;

  field::FaceField<double> velocity(mesh_, 0.0);

  field::FaceField<double> flux(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  const auto face = findInternalFace(mesh_);

  const auto owner = mesh_.owner(face);

  velocity[face] = 3.0;

  discretization::operators::upwindAdvectionFlux(u, velocity, bc, flux);

  EXPECT_NEAR(flux[face], 3.0 * u[owner], 1e-12);
}

TEST_F(OperatorTest, UpwindFluxUsesNeighborForNegativeVelocity) {
  field::CellField<double> u(mesh_, 0.0);

  u[0] = 2.0;
  u[1] = 5.0;

  field::FaceField<double> velocity(mesh_, 0.0);

  field::FaceField<double> flux(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  const auto face = findInternalFace(mesh_);

  const auto neighbor = mesh_.neighbor(face);

  velocity[face] = -3.0;

  discretization::operators::upwindAdvectionFlux(u, velocity, bc, flux);

  EXPECT_NEAR(flux[face], -3.0 * u[neighbor], 1e-12);
}

TEST_F(OperatorTest, UpwindFluxUsesInteriorStateForOutflow) {
  field::CellField<double> u(mesh_, 4.0);

  field::FaceField<double> velocity(mesh_, 0.0);

  field::FaceField<double> flux(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  const auto face = findBoundaryFace(mesh_, mesh::BoundaryId{2});

  const auto owner = mesh_.owner(face);

  velocity[face] = +2.0;

  //
  // Intentionally no BC.
  //
  discretization::operators::upwindAdvectionFlux(u, velocity, bc, flux);

  EXPECT_NEAR(flux[face], 2.0 * u[owner], 1e-12);
}

TEST_F(OperatorTest, UpwindFluxUsesBoundaryStateForInflow) {
  field::CellField<double> u(mesh_, 4.0);

  field::FaceField<double> velocity(mesh_, 0.0);

  field::FaceField<double> flux(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(mesh::BoundaryId{1}, 7.0);

  const auto face = findBoundaryFace(mesh_, mesh::BoundaryId{1});

  velocity[face] = -2.0;

  discretization::operators::upwindAdvectionFlux(u, velocity, bc, flux);

  EXPECT_NEAR(flux[face], -2.0 * 7.0, 1e-12);
}

TEST_F(OperatorTest, UpwindFluxRejectsMissingInflowBoundaryValue) {
  field::CellField<double> u(mesh_, 1.0);

  field::FaceField<double> velocity(mesh_, 0.0);

  field::FaceField<double> flux(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  const auto face = findBoundaryFace(mesh_, mesh::BoundaryId{1});

  velocity[face] = -1.0;

  EXPECT_THROW(
      discretization::operators::upwindAdvectionFlux(u, velocity, bc, flux),
      std::runtime_error);
}

TEST_F(OperatorTest, UpwindFluxRejectsNeumannOnInflow) {
  field::CellField<double> u(mesh_, 1.0);

  field::FaceField<double> velocity(mesh_, 0.0);

  field::FaceField<double> flux(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  bc.setNeumann(mesh::BoundaryId{1}, 0.0);

  const auto face = findBoundaryFace(mesh_, mesh::BoundaryId{1});

  velocity[face] = -1.0;

  EXPECT_THROW(
      discretization::operators::upwindAdvectionFlux(u, velocity, bc, flux),
      std::runtime_error);
}

TEST_F(OperatorTest, ConstantStateHasZeroAdvectionDivergence) {
  field::CellField<double> u(mesh_, 3.0);

  field::FaceField<double> velocity(mesh_, 0.0);

  field::FaceField<double> flux(mesh_, 0.0);

  field::CellField<double> div(mesh_, 0.0);

  // --------------------------------------------------------
  // Constant velocity:
  //
  //     v = (1, 0)
  //
  // therefore:
  //
  //     v_n = v · n = n_x
  // --------------------------------------------------------

  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {

    const auto normal = mesh_.faceNormal(face);

    velocity[face] = normal.x;
  }

  // --------------------------------------------------------
  // Left boundary is inflow:
  //
  //     u_in = 3
  //
  // same as interior constant solution.
  //
  // Other boundaries do not need values.
  // --------------------------------------------------------

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(mesh::BoundaryId{1}, 3.0);

  discretization::operators::upwindAdvectionFlux(u, velocity, bc, flux);

  discretization::operators::divergence(flux, div);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(div[cell], 0.0, 1e-12);
  }
}

TEST_F(OperatorTest, AdvectionFluxIsGloballyConservative) {
  field::CellField<double> u(mesh_, 0.0);

  u[0] = 2.0;
  u[1] = 5.0;

  field::FaceField<double> velocity(mesh_, 0.0);

  field::FaceField<double> flux(mesh_, 0.0);

  field::CellField<double> div(mesh_, 0.0);

  //
  // Only internal transport.
  //
  const auto face = findInternalFace(mesh_);

  velocity[face] = 1.0;

  boundary::BoundaryConditionSet bc;

  discretization::operators::upwindAdvectionFlux(u, velocity, bc, flux);

  discretization::operators::divergence(flux, div);

  double total = 0.0;

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    total += div[cell] * mesh_.cellVolume(cell);
  }

  EXPECT_NEAR(total, 0.0, 1e-12);
}

TEST_F(OperatorTest, UpwindFluxRejectsDifferentMeshes) {
  mesh::MoabMesh other_mesh(testMeshPath().string());

  field::CellField<double> u(mesh_, 1.0);

  field::FaceField<double> velocity(other_mesh, 0.0);

  field::FaceField<double> flux(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  EXPECT_THROW(
      discretization::operators::upwindAdvectionFlux(u, velocity, bc, flux),
      std::invalid_argument);
}

TEST(BernoulliTest, ValueAtZeroIsOne) {
  EXPECT_DOUBLE_EQ(discretization::operators::bernoulli(0.0), 1.0);
}

TEST(BernoulliTest, MatchesKnownValues) {
  EXPECT_NEAR(discretization::operators::bernoulli(1.0), 0.5819767068693265,
              1e-14);
  EXPECT_NEAR(discretization::operators::bernoulli(-1.0), 1.5819767068693265,
              1e-14);
}

TEST(BernoulliTest, SatisfiesDifferenceIdentity) {
  for (const double x : {-20.0, -5.0, -1.0, -0.01, 0.01, 1.0, 5.0, 20.0}) {
    EXPECT_NEAR(discretization::operators::bernoulli(-x) -
                    discretization::operators::bernoulli(x),
                x, 1e-12);
  }
}

TEST(BernoulliTest, IsStableNearZero) {
  const double x = 1e-12;

  const double expected = 1.0 - x / 2.0;

  EXPECT_NEAR(discretization::operators::bernoulli(x), expected, 1e-15);
}

TEST(BernoulliTest, HasCorrectLargeArgumentLimits) {
  EXPECT_LT(discretization::operators::bernoulli(100.0), 1e-40);
  EXPECT_NEAR(discretization::operators::bernoulli(-100.0), 100.0, 1e-12);
}

TEST_F(OperatorTest, ScharfetterGummelReducesToDiffusionAtZeroVelocity) {
  field::CellField<double> u(mesh_, 0.0);

  u[0] = 2.0;
  u[1] = 5.0;

  field::FaceField<double> velocity(mesh_, 0.0);

  field::FaceField<double> sg_flux(mesh_, 0.0);

  field::FaceField<double> diffusion_flux(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(1, 0.0);
  bc.setDirichlet(2, 0.0);
  bc.setDirichlet(3, 0.0);
  bc.setDirichlet(4, 0.0);

  constexpr double D = 2.0;

  discretization::operators::scharfetterGummelFlux(u, velocity, D, bc, sg_flux);

  discretization::operators::diffusionFlux(u, D, bc, diffusion_flux);

  const auto face = findInternalFace(mesh_);

  EXPECT_NEAR(sg_flux[face], diffusion_flux[face], 1e-12);
}

TEST_F(OperatorTest,
       ScharfetterGummelApproachesOwnerUpwindForStrongPositiveDrift) {
  field::CellField<double> u(mesh_, 0.0);

  const auto face = findInternalFace(mesh_);

  const auto owner = mesh_.owner(face);

  const auto neighbor = mesh_.neighbor(face);

  u[owner] = 2.0;

  u[neighbor] = 5.0;

  field::FaceField<double> velocity(mesh_, 0.0);

  velocity[face] = 1.0;

  field::FaceField<double> flux(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(1, 0.0);
  bc.setDirichlet(2, 0.0);
  bc.setDirichlet(3, 0.0);
  bc.setDirichlet(4, 0.0);

  discretization::operators::scharfetterGummelFlux(u, velocity, 1e-3, bc, flux);

  EXPECT_NEAR(flux[face], 1.0 * u[owner], 1e-10);
}

TEST_F(OperatorTest,
       ScharfetterGummelApproachesNeighborUpwindForStrongNegativeDrift) {
  field::CellField<double> u(mesh_, 0.0);

  const auto face = findInternalFace(mesh_);

  const auto owner = mesh_.owner(face);

  const auto neighbor = mesh_.neighbor(face);

  u[owner] = 2.0;

  u[neighbor] = 5.0;

  field::FaceField<double> velocity(mesh_, 0.0);

  velocity[face] = -1.0;

  field::FaceField<double> flux(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(1, 0.0);
  bc.setDirichlet(2, 0.0);
  bc.setDirichlet(3, 0.0);
  bc.setDirichlet(4, 0.0);

  discretization::operators::scharfetterGummelFlux(u, velocity, 1e-3, bc, flux);

  EXPECT_NEAR(flux[face], -1.0 * u[neighbor], 1e-10);
}

TEST_F(OperatorTest, ScharfetterGummelPreservesConstantStateFlux) {
  constexpr double value = 3.0;

  constexpr double D = 0.2;

  field::CellField<double> u(mesh_, value);

  field::FaceField<double> velocity(mesh_, 0.0);

  field::FaceField<double> flux(mesh_, 0.0);

  //
  // v = (1, 0)
  //
  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {

    const auto normal = mesh_.faceNormal(face);

    velocity[face] = normal.x;
  }

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(1, value);
  bc.setDirichlet(2, value);
  bc.setDirichlet(3, value);
  bc.setDirichlet(4, value);

  discretization::operators::scharfetterGummelFlux(u, velocity, D, bc, flux);

  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {

    EXPECT_NEAR(flux[face], velocity[face] * value, 1e-12);
  }
}

TEST_F(OperatorTest, ConstantSgFluxHasZeroDivergence) {
  constexpr double value = 3.0;

  field::CellField<double> u(mesh_, value);

  field::FaceField<double> velocity(mesh_, 0.0);

  field::FaceField<double> flux(mesh_, 0.0);

  field::CellField<double> div(mesh_, 0.0);

  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {

    velocity[face] = mesh_.faceNormal(face).x;
  }

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(1, value);
  bc.setDirichlet(2, value);
  bc.setDirichlet(3, value);
  bc.setDirichlet(4, value);

  discretization::operators::scharfetterGummelFlux(u, velocity, 0.2, bc, flux);

  discretization::operators::divergence(flux, div);

  for (mesh::CellId cell = 0; cell < mesh_.numCells(); ++cell) {

    EXPECT_NEAR(div[cell], 0.0, 1e-12);
  }
}

TEST_F(OperatorTest, ScharfetterGummelExactlyPreservesExponentialEquilibrium) {
  const auto face = findInternalFace(mesh_);

  const auto owner = mesh_.owner(face);

  const auto neighbor = mesh_.neighbor(face);

  constexpr double D = 1.0;

  constexpr double vn = 0.5;

  const auto delta = mesh_.cellCenter(neighbor) - mesh_.cellCenter(owner);

  const double distance = mesh::dot(delta, mesh_.faceNormal(face));

  const double pe = vn * distance / D;

  field::CellField<double> u(mesh_, 0.0);

  u[owner] = 1.0;

  u[neighbor] = std::exp(pe);

  field::FaceField<double> velocity(mesh_, 0.0);

  velocity[face] = vn;

  field::FaceField<double> flux(mesh_, 0.0);

  boundary::BoundaryConditionSet bc;

  bc.setDirichlet(1, 0.0);
  bc.setDirichlet(2, 0.0);
  bc.setDirichlet(3, 0.0);
  bc.setDirichlet(4, 0.0);

  discretization::operators::scharfetterGummelFlux(u, velocity, D, bc, flux);

  EXPECT_NEAR(flux[face], 0.0, 1e-12);
}

TEST_F(OperatorTest, ScharfetterGummelNeumannAddsPrescribedDiffusiveFlux) {
  constexpr double state_value = 3.0;
  constexpr double normal_velocity = 2.0;
  constexpr double diffusive_flux = -0.5;
  field::CellField<double> state(mesh_, state_value);
  field::FaceField<double> velocity(mesh_, 0.0);
  field::FaceField<double> flux(mesh_, 0.0);
  const auto face = findBoundaryFace(mesh_, mesh::BoundaryId{1});
  velocity[face] = normal_velocity;
  boundary::BoundaryConditionSet conditions;
  conditions.setNeumann(mesh::BoundaryId{1}, diffusive_flux);
  conditions.setDirichlet(mesh::BoundaryId{2}, state_value);
  conditions.setDirichlet(mesh::BoundaryId{3}, state_value);
  conditions.setDirichlet(mesh::BoundaryId{4}, state_value);

  discretization::operators::scharfetterGummelFlux(state, velocity, 0.2,
                                                   conditions, flux);

  EXPECT_DOUBLE_EQ(flux[face], normal_velocity * state_value + diffusive_flux);
}

TEST_F(OperatorTest, ScharfetterGummelHomogeneousNeumannKeepsDriftFlux) {
  constexpr double state_value = 3.0;
  field::CellField<double> state(mesh_, state_value);
  field::FaceField<double> velocity(mesh_, 0.0);
  field::FaceField<double> flux(mesh_, 0.0);
  boundary::BoundaryConditionSet conditions;
  for (mesh::BoundaryId boundary = 1; boundary <= 4; ++boundary) {
    conditions.setNeumann(boundary, 0.0);
  }
  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {
    velocity[face] = mesh_.faceNormal(face).x;
  }

  discretization::operators::scharfetterGummelFlux(state, velocity, 0.2,
                                                   conditions, flux);

  for (mesh::FaceId face = 0; face < mesh_.numFaces(); ++face) {
    if (mesh_.isBoundary(face)) {
      EXPECT_DOUBLE_EQ(flux[face], velocity[face] * state_value);
    }
  }
}

};  // namespace pemu::discretization::test
