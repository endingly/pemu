#pragma once

#include <pemu/boundary/boundary_condition_set.hpp>
#include <pemu/discretization/operators/bernoulli.hpp>
#include <pemu/discretization/operators/divergence.hpp>
#include <pemu/discretization/operators/scharfetter_gummel_flux.hpp>
#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>
#include <pemu/mesh/geometry.hpp>
#include <pemu/mesh/i_mesh.hpp>

#include <algorithm>
#include <stdexcept>

namespace pemu::equation {

// The original fixed-timestep Scharfetter--Gummel continuity stepper.
// dt is immutable state of this object and every call to step() uses it.
class FixedStepExplicitSpeciesContinuityStepper {
 public:
  FixedStepExplicitSpeciesContinuityStepper(
      const mesh::IMesh& mesh,
      const field::FaceField<double>& normal_drift_velocity, double diffusivity,
      double dt, const boundary::BoundaryConditionSet& bc)
      : mesh_(&mesh),
        normal_drift_velocity_(&normal_drift_velocity),
        bc_(&bc),
        diffusivity_(diffusivity),
        dt_(dt),
        flux_(mesh, 0.0),
        divergence_(mesh, 0.0) {
    if (&normal_drift_velocity.mesh() != &mesh) {
      throw std::invalid_argument("drift velocity belongs to another mesh");
    }
    if (diffusivity <= 0.0) {
      throw std::invalid_argument("diffusivity must be positive");
    }
    if (dt <= 0.0) {
      throw std::invalid_argument("time step must be positive");
    }
  }

  [[nodiscard]] double maxTransportCfl() const {
    field::CellField<double> diagonal(*mesh_, 0.0);

    for (mesh::FaceId face = 0; face < mesh_->numFaces(); ++face) {
      const auto owner = mesh_->owner(face);
      const auto normal = mesh_->faceNormal(face);
      const double vn = (*normal_drift_velocity_)[face];

      if (!mesh_->isBoundary(face)) {
        const auto neighbor = mesh_->neighbor(face);
        const double distance = mesh::dot(
            mesh_->cellCenter(neighbor) - mesh_->cellCenter(owner), normal);
        if (distance <= 0.0) {
          throw std::runtime_error("invalid internal face distance");
        }
        const double pe = vn * distance / diffusivity_;
        const double scale = diffusivity_ * mesh_->faceArea(face) / distance;
        diagonal[owner] += scale * discretization::operators::bernoulli(-pe);
        diagonal[neighbor] += scale * discretization::operators::bernoulli(pe);
        continue;
      }

      const double distance =
          mesh::dot(mesh_->faceCenter(face) - mesh_->cellCenter(owner), normal);
      if (distance <= 0.0) {
        throw std::runtime_error("invalid boundary distance");
      }
      const double pe = vn * distance / diffusivity_;
      const double scale = diffusivity_ * mesh_->faceArea(face) / distance;
      diagonal[owner] += scale * discretization::operators::bernoulli(-pe);
    }

    double maximum = 0.0;
    for (mesh::CellId cell = 0; cell < mesh_->numCells(); ++cell) {
      const double local = dt_ * diagonal[cell] / mesh_->cellVolume(cell);
      maximum = std::max(maximum, local);
    }
    return maximum;
  }

  [[nodiscard]] const mesh::IMesh& mesh() const noexcept { return *mesh_; }
  [[nodiscard]] double timeStep() const noexcept { return dt_; }
  [[nodiscard]] double diffusivity() const noexcept { return diffusivity_; }

  void step(field::CellField<double>& density,
            const field::CellField<double>& source) {
    validateFields(density, source);
    discretization::operators::scharfetterGummelFlux(
        density, *normal_drift_velocity_, diffusivity_, *bc_, flux_);
    discretization::operators::divergence(flux_, divergence_);

    for (mesh::CellId cell = 0; cell < mesh_->numCells(); ++cell) {
      density[cell] += dt_ * (source[cell] - divergence_[cell]);
    }
  }

 private:
  void validateFields(const field::CellField<double>& density,
                      const field::CellField<double>& source) const {
    if (&density.mesh() != mesh_) {
      throw std::invalid_argument("density belongs to another mesh");
    }
    if (&source.mesh() != mesh_) {
      throw std::invalid_argument("source belongs to another mesh");
    }
  }

  const mesh::IMesh* mesh_;
  const field::FaceField<double>* normal_drift_velocity_;
  const boundary::BoundaryConditionSet* bc_;
  double diffusivity_;
  double dt_;
  field::FaceField<double> flux_;
  field::CellField<double> divergence_;
};

}  // namespace pemu::equation
