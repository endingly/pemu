#pragma once

#include <pemu/boundary/boundary_condition_set.hpp>
#include <pemu/field/cell_field.hpp>
#include <pemu/linalg/types.hpp>
#include <pemu/mesh/i_mesh.hpp>

namespace pemu::discretization {

class PoissonFvm {
 public:
  PoissonFvm(const mesh::IMesh& mesh, const field::CellField<double>& source,
             double epsilon,
             const boundary::BoundaryConditionSet& boundary_conditions);

  [[nodiscard]]
  linalg::LinearSystem assemble() const;

  [[nodiscard]]
  const mesh::IMesh& mesh() const noexcept {
    return *mesh_;
  }

 private:
  [[nodiscard]]
  double internalFaceCoefficient(mesh::FaceId face) const;

  [[nodiscard]]
  double boundaryFaceCoefficient(mesh::FaceId face) const;

 private:
  const mesh::IMesh* mesh_;
  const field::CellField<double>* source_;
  const boundary::BoundaryConditionSet* boundary_conditions_;

  double epsilon_;
};

}  // namespace pemu::discretization