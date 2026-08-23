#include <gtest/gtest.h>

#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>
#include <pemu/field/field_set.hpp>
#include <pemu/field/quantity_io.hpp>
#include <pemu/mesh/moab_mesh.hpp>
#include <pemu/unit/plasma_quantities.hpp>

#include <mp-units/systems/si.h>

#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <type_traits>

namespace pemu::field::test {

namespace {

std::filesystem::path testMeshPath() {
#ifdef PEMU_MESH_TEST_DATA_DIR
  return std::filesystem::path{PEMU_MESH_TEST_DATA_DIR} / "two_quads.msh";
#else
  return std::filesystem::path{"two_quads.msh"};
#endif
}

class FieldTest : public ::testing::Test {
 protected:
  FieldTest() : mesh_(testMeshPath().string()) {}

  mesh::MoabMesh mesh_;
};

struct TestFieldId {
  std::uint32_t value{};
};

template <typename QS, typename U>
concept FormsQuantityReference =
    requires(std::remove_cvref_t<QS> quantity_spec,
             std::remove_cvref_t<U> unit) { quantity_spec[unit]; };

}  // namespace

// ------------------------------------------------------------
// CellField
// ------------------------------------------------------------

TEST_F(FieldTest, CellFieldSizeMatchesMesh) {
  CellField<double> field(mesh_);

  EXPECT_EQ(field.size(), mesh_.numCells());

  EXPECT_EQ(field.size(), 2u);
}

TEST_F(FieldTest, CellFieldSupportsInitialValue) {
  CellField<double> field(mesh_, 3.5);

  for (mesh::CellId c = 0; c < mesh_.numCells(); ++c) {

    EXPECT_DOUBLE_EQ(field[c], 3.5);
  }
}

TEST_F(FieldTest, CellFieldSupportsDenseIdAccess) {
  CellField<double> field(mesh_);

  field[0] = 10.0;
  field[1] = 20.0;

  EXPECT_DOUBLE_EQ(field[0], 10.0);

  EXPECT_DOUBLE_EQ(field[1], 20.0);
}

TEST_F(FieldTest, CellFieldFillWorks) {
  CellField<double> field(mesh_);

  field.fill(7.0);

  for (const auto value : field) {

    EXPECT_DOUBLE_EQ(value, 7.0);
  }
}

TEST_F(FieldTest, CellFieldAtRejectsInvalidCell) {
  CellField<double> field(mesh_);

  EXPECT_THROW(field.at(static_cast<mesh::CellId>(mesh_.numCells())),
               std::out_of_range);
}

TEST_F(FieldTest, CellFieldKeepsMeshAssociation) {
  CellField<double> field(mesh_);

  EXPECT_EQ(&field.mesh(), &mesh_);
}

// ------------------------------------------------------------
// FaceField
// ------------------------------------------------------------

TEST_F(FieldTest, FaceFieldSizeMatchesMesh) {
  FaceField<double> field(mesh_);

  EXPECT_EQ(field.size(), mesh_.numFaces());

  EXPECT_EQ(field.size(), 7u);
}

TEST_F(FieldTest, FaceFieldSupportsInitialValue) {
  FaceField<double> field(mesh_, 2.0);

  for (mesh::FaceId f = 0; f < mesh_.numFaces(); ++f) {

    EXPECT_DOUBLE_EQ(field[f], 2.0);
  }
}

TEST_F(FieldTest, FaceFieldSupportsDenseIdAccess) {
  FaceField<double> field(mesh_);

  for (mesh::FaceId f = 0; f < mesh_.numFaces(); ++f) {

    field[f] = static_cast<double>(f);
  }

  for (mesh::FaceId f = 0; f < mesh_.numFaces(); ++f) {

    EXPECT_DOUBLE_EQ(field[f], static_cast<double>(f));
  }
}

TEST_F(FieldTest, FaceFieldAtRejectsInvalidFace) {
  FaceField<double> field(mesh_);

  EXPECT_THROW(field.at(static_cast<mesh::FaceId>(mesh_.numFaces())),
               std::out_of_range);
}

TEST_F(FieldTest, FaceFieldKeepsMeshAssociation) {
  FaceField<double> field(mesh_);

  EXPECT_EQ(&field.mesh(), &mesh_);
}

TEST_F(FieldTest, FieldSetGroupsTypedFieldsAndPropagatesMetadata) {
  using namespace mp_units;
  using namespace mp_units::si::unit_symbols;

  const auto metadata =
      makeFieldMetadata("inverse length", isq::repetency[one / m]);
  CellFieldSet<double, TestFieldId> fields(mesh_, 2, 1.5, metadata);

  EXPECT_EQ(fields.size(), 2u);
  EXPECT_FALSE(fields.empty());
  EXPECT_EQ(&fields.mesh(), &mesh_);
  EXPECT_EQ(fields.span().size(), 2u);
  const auto expected = pemu::unit::bridgeReference(isq::repetency[one / m]);
  EXPECT_EQ(*fields[TestFieldId{0}].metadata().physical_quantity, expected);
  EXPECT_EQ(*fields[TestFieldId{1}].metadata().physical_quantity, expected);

  fields[TestFieldId{0}][mesh::CellId{0}] = 9.0;
  EXPECT_DOUBLE_EQ(fields[TestFieldId{0}][mesh::CellId{0}], 9.0);
  EXPECT_DOUBLE_EQ(fields[TestFieldId{1}][mesh::CellId{0}], 1.5);

  fields.fill(4.0);
  for (const auto& field : fields) {
    for (const auto value : field) {
      EXPECT_DOUBLE_EQ(value, 4.0);
    }
  }

  FaceFieldSet<double, TestFieldId> face_fields(mesh_, 1, 2.0, metadata);
  EXPECT_EQ(face_fields[TestFieldId{0}].size(), mesh_.numFaces());
  EXPECT_DOUBLE_EQ(face_fields[TestFieldId{0}][mesh::FaceId{0}], 2.0);
}

// ------------------------------------------------------------
// Unit metadata and mp-units boundary conversion
// ------------------------------------------------------------

TEST_F(FieldTest, PhysicalQuantityIsMetadataAndRawStorageRemainsDouble) {
  using namespace mp_units;
  using namespace mp_units::si::unit_symbols;

  static_assert(
      FormsQuantityReference<decltype(isq::electric_potential), decltype(V)>);
  static_assert(!FormsQuantityReference<decltype(isq::electric_potential),
                                        decltype(cm / s)>);

  CellField<double> field(
      mesh_, 0.0, makeFieldMetadata("potential", isq::electric_potential[V]));

  static_assert(std::same_as<CellField<double>::value_type, double>);
  static_assert(std::same_as<decltype(field.data()), double*>);

  ASSERT_TRUE(field.metadata().hasPhysicalQuantity());
  EXPECT_EQ(field.metadata().name, "potential");
  EXPECT_EQ(field.metadata().physical_quantity->kind(),
            pemu::unit::QuantityKind::electric_potential);
  EXPECT_EQ(field.metadata().physical_quantity->unit(), units::precise::V);
  static_assert(
      std::is_trivially_copyable_v<pemu::unit::PhysicalQuantityMetadata>);
}

TEST_F(FieldTest, MpUnitsBridgeCoversCanonicalPlasmaFieldMetadata) {
  using namespace mp_units;
  using namespace mp_units::si::unit_symbols;

  constexpr auto density = pemu::unit::bridgeReference(
      pemu::unit::plasma_quantity::particle_number_density[one / cubic(cm)]);
  constexpr auto charge =
      pemu::unit::bridgeReference(isq::electric_charge_density[C / cubic(cm)]);
  constexpr auto electric_field = pemu::unit::bridgeReference(
      pemu::unit::plasma_quantity::normal_electric_field_strength[V / cm]);
  constexpr auto energy_density = pemu::unit::bridgeReference(
      pemu::unit::plasma_quantity::electron_energy_density[eV / cubic(cm)]);

  static_assert(density.kind() ==
                pemu::unit::QuantityKind::particle_number_density);
  EXPECT_EQ(density.unit(), units::precise::one / units::precise::cm.pow(3));
  static_assert(charge.kind() ==
                pemu::unit::QuantityKind::electric_charge_density);
  EXPECT_EQ(charge.unit(), units::precise::C / units::precise::cm.pow(3));
  static_assert(electric_field.kind() ==
                pemu::unit::QuantityKind::normal_electric_field_strength);
  EXPECT_EQ(electric_field.unit(), units::precise::V / units::precise::cm);
  static_assert(energy_density.kind() ==
                pemu::unit::QuantityKind::electron_energy_density);
  EXPECT_EQ(energy_density.unit(),
            units::precise::energy::eV / units::precise::cm.pow(3));

  EXPECT_EQ(pemu::to_string(density.kind()), "particle_number_density");
  EXPECT_EQ(pemu::to_string(charge.kind()), "electric_charge_density");
  EXPECT_EQ(pemu::to_string(electric_field.kind()),
            "normal_electric_field_strength");
  EXPECT_EQ(pemu::to_string(energy_density.kind()), "electron_energy_density");
}

TEST_F(FieldTest, QuantityBoundaryConvertsToFieldStorageUnit) {
  using namespace mp_units;
  using namespace mp_units::si::unit_symbols;

  constexpr auto centimetre_length = isq::length[cm];
  CellField<double> field(mesh_, 0.0,
                          makeFieldMetadata("x", centimetre_length));

  fillQuantity(field, 3.0 * m, centimetre_length);
  EXPECT_DOUBLE_EQ(field[0], 300.0);
  EXPECT_DOUBLE_EQ(field[1], 300.0);

  setQuantity(field, mesh::CellId{0}, 2.5 * m, centimetre_length);
  setQuantity(field, mesh::CellId{1}, 12.0 * mm, centimetre_length);

  EXPECT_DOUBLE_EQ(field[0], 250.0);
  EXPECT_DOUBLE_EQ(field[1], 1.2);
  const auto value = quantityAt(field, mesh::CellId{0}, centimetre_length);
  EXPECT_DOUBLE_EQ(value.numerical_value_ref_in(cm), 250.0);
}

TEST_F(FieldTest, QuantityBoundaryRejectsReferenceDifferentFromMetadata) {
  using namespace mp_units;
  using namespace mp_units::si::unit_symbols;

  FaceField<double> field(
      mesh_, 0.0,
      makeFieldMetadata(
          "normal electric field",
          pemu::unit::plasma_quantity::normal_electric_field_strength[V / cm]));

  EXPECT_THROW((void)quantityAt(field, mesh::FaceId{0}, isq::speed[m / s]),
               std::invalid_argument);

  CellField<double> potential(
      mesh_, 0.0, makeFieldMetadata("potential", isq::electric_potential[V]));
  EXPECT_THROW((void)quantityAt(potential, mesh::CellId{0},
                                isq::electric_potential_difference[V]),
               std::invalid_argument);
}

}  // namespace pemu::field::test
