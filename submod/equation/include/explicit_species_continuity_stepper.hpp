#pragma once

#include <pemu/boundary/boundary_condition_set.hpp>

#include <pemu/discretization/operators/divergence.hpp>
#include <pemu/discretization/operators/scharfetter_gummel_flux.hpp>

#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>

#include <pemu/mesh/geometry.hpp>
#include <pemu/mesh/i_mesh.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace pemu::equation {

class ExplicitSpeciesContinuityStepper {
 public:
  ExplicitSpeciesContinuityStepper(
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

      throw std::invalid_argument(
          "drift velocity belongs "
          "to another mesh");
    }

    if (diffusivity <= 0.0) {
      throw std::invalid_argument("diffusivity must be positive");
    }

    if (dt <= 0.0) {
      throw std::invalid_argument("time step must be positive");
    }
  }

  [[nodiscard]]
  double maxTransportCfl() const;

  [[nodiscard]]
  const mesh::IMesh& mesh() const noexcept {
    return *mesh_;
  }

  [[nodiscard]]
  double timeStep() const noexcept {
    return dt_;
  }

  [[nodiscard]]
  double diffusivity() const noexcept {
    return diffusivity_;
  }

  // ========================================================
  // One explicit timestep:
  //
  //     dn/dt + div(Gamma) = S
  //
  // therefore:
  //
  //     n^{k+1}
  //
  //       = n^k
  //         - dt div(Gamma^k)
  //         + dt S^k
  //
  // Gamma is evaluated using Scharfetter-Gummel.
  // ========================================================

  void step(field::CellField<double>& density,
            const field::CellField<double>& source) {
    validateFields(density, source);

    // ----------------------------------------------------
    // Face drift-diffusion flux:
    //
    //     Gamma
    //
    //       = v n
    //         - D grad(n)
    // ----------------------------------------------------

    discretization::operators::scharfetterGummelFlux(
        density, *normal_drift_velocity_, diffusivity_, *bc_, flux_);

    // ----------------------------------------------------
    // div(Gamma)
    // ----------------------------------------------------

    discretization::operators::divergence(flux_, divergence_);

    // ----------------------------------------------------
    // Explicit Euler update.
    // ----------------------------------------------------

    for (mesh::CellId cell = 0; cell < mesh_->numCells(); ++cell) {

      density[cell] += dt_ * (source[cell] - divergence_[cell]);
    }
  }

 private:
  void validateFields(const field::CellField<double>& density,
                      const field::CellField<double>& source) const {
    if (&density.mesh() != mesh_) {
      throw std::invalid_argument(
          "density belongs to "
          "another mesh");
    }

    if (&source.mesh() != mesh_) {
      throw std::invalid_argument(
          "source belongs to "
          "another mesh");
    }
  }

 private:
  const mesh::IMesh* mesh_;

  const field::FaceField<double>* normal_drift_velocity_;

  const boundary::BoundaryConditionSet* bc_;

  double diffusivity_;
  double dt_;

  field::FaceField<double> flux_;

  field::CellField<double> divergence_;
};

}  // namespace pemu::equation