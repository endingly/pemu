#pragma once

#include <pemu/boundary/boundary_condition_set.hpp>
#include <pemu/discretization/operators/divergence.hpp>
#include <pemu/discretization/operators/upwind_advection_flux.hpp>
#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>

#include <algorithm>
#include <stdexcept>

namespace pemu::equation {

class ExplicitEulerAdvectionStepper {
 public:
  ExplicitEulerAdvectionStepper(const mesh::IMesh& mesh,
                                const field::FaceField<double>& normal_velocity,
                                double dt,
                                const boundary::BoundaryConditionSet& bc)
      : mesh_(&mesh),
        normal_velocity_(&normal_velocity),
        bc_(&bc),
        dt_(dt),
        flux_(mesh, 0.0),
        divergence_(mesh, 0.0) {
    if (&normal_velocity.mesh() != &mesh) {
      throw std::invalid_argument("velocity field belongs to another mesh");
    }

    if (dt <= 0.0) {
      throw std::invalid_argument("time step must be positive");
    }
  }

  [[nodiscard]]
  double maxCfl() const {
    field::CellField<double> outgoing(*mesh_, 0.0);

    for (mesh::FaceId face = 0; face < mesh_->numFaces(); ++face) {

      const double vn = (*normal_velocity_)[face];

      const double area = mesh_->faceArea(face);

      const auto owner = mesh_->owner(face);

      //
      // Relative to owner:
      //
      // vn > 0 means outgoing.
      //
      if (vn > 0.0) {
        outgoing[owner] += vn * area;
      }

      if (!mesh_->isBoundary(face)) {

        const auto neighbor = mesh_->neighbor(face);

        //
        // Relative to neighbor the outward
        // normal is -n_f.
        //
        // Therefore vn < 0 means outgoing
        // from neighbor.
        //
        if (vn < 0.0) {
          outgoing[neighbor] += (-vn) * area;
        }
      }
    }

    double maximum = 0.0;

    for (mesh::CellId cell = 0; cell < mesh_->numCells(); ++cell) {

      const double local = dt_ * outgoing[cell] / mesh_->cellVolume(cell);

      maximum = std::max(maximum, local);
    }

    return maximum;
  }

  void step(field::CellField<double>& state) {
    if (&state.mesh() != mesh_) {
      throw std::invalid_argument("state belongs to another mesh");
    }

    if (maxCfl() > 1.0) {
      throw std::runtime_error("explicit advection CFL condition violated");
    }

    //
    // Gamma^n
    //
    discretization::operators::upwindAdvectionFlux(state, *normal_velocity_,
                                                   *bc_, flux_);

    //
    // div(Gamma^n)
    //
    discretization::operators::divergence(flux_, divergence_);

    //
    // u^{n+1}
    //
    for (mesh::CellId cell = 0; cell < mesh_->numCells(); ++cell) {

      state[cell] -= dt_ * divergence_[cell];
    }
  }

 private:
  const mesh::IMesh* mesh_;

  const field::FaceField<double>* normal_velocity_;

  const boundary::BoundaryConditionSet* bc_;

  double dt_;

  field::FaceField<double> flux_;

  field::CellField<double> divergence_;
};

}  // namespace pemu::equation