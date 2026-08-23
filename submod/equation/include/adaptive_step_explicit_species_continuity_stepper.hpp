#pragma once

#include <pemu/boundary/boundary_condition.hpp>
#include <pemu/boundary/boundary_condition_set.hpp>

#include <pemu/discretization/operators/bernoulli.hpp>
#include <pemu/discretization/operators/divergence.hpp>
#include <pemu/discretization/operators/linear_boundary_flux.hpp>
#include <pemu/discretization/operators/scharfetter_gummel_flux.hpp>

#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>
#include <pemu/field/utils.hpp>

#include <pemu/mesh/geometry.hpp>
#include <pemu/mesh/i_mesh.hpp>

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <variant>

namespace pemu::equation {

class AdaptiveStepExplicitSpeciesContinuityStepper {
 public:
  AdaptiveStepExplicitSpeciesContinuityStepper(
      const mesh::IMesh& mesh,
      const field::FaceField<double>& normal_drift_velocity, double diffusivity,
      const boundary::BoundaryConditionSet& bc)
      : mesh_(&mesh),
        normal_drift_velocity_(&normal_drift_velocity),
        bc_(&bc),
        diffusivity_(diffusivity),
        flux_(mesh, 0.0),
        divergence_(mesh, 0.0),
        local_increment_(mesh, 0.0),
        wall_loss_rate_(mesh, 0.0) {
    if (&normal_drift_velocity.mesh() != &mesh) {
      throw std::invalid_argument("drift velocity belongs to another mesh");
    }

    if (diffusivity <= 0.0) {
      throw std::invalid_argument("diffusivity must be positive");
    }
  }

  /** @brief Creates an SG stepper with linear wall flux on selected faces. */
  AdaptiveStepExplicitSpeciesContinuityStepper(
      const mesh::IMesh& mesh,
      const field::FaceField<double>& normal_drift_velocity, double diffusivity,
      const boundary::BoundaryConditionSet& bc,
      discretization::operators::LinearBoundaryFluxView wall_flux)
      : AdaptiveStepExplicitSpeciesContinuityStepper(
            mesh, normal_drift_velocity, diffusivity, bc) {
    if (&wall_flux.mesh() != &mesh) {
      throw std::invalid_argument("wall flux belongs to another mesh");
    }
    wall_flux_.emplace(wall_flux);
  }

  [[nodiscard]]
  const mesh::IMesh& mesh() const noexcept {
    return *mesh_;
  }

  [[nodiscard]]
  double diffusivity() const noexcept {
    return diffusivity_;
  }

  /** @brief Computes the SG particle flux used by the continuity update. */
  void computeNormalFlux(const field::CellField<double>& density,
                         field::FaceField<double>& normal_flux) const {
    validateField(density);
    if (&normal_flux.mesh() != mesh_) {
      throw std::invalid_argument("normal flux belongs to another mesh");
    }
    if (wall_flux_.has_value()) {
      discretization::operators::scharfetterGummelFlux(
          density, *normal_drift_velocity_, diffusivity_, *bc_, *wall_flux_,
          normal_flux);
    } else {
      discretization::operators::scharfetterGummelFlux(
          density, *normal_drift_velocity_, diffusivity_, *bc_, normal_flux);
    }
  }

  // ========================================================
  // Compute:
  //
  //     lambda_P
  //
  // where the explicit transport operator contains:
  //
  //     -lambda_P n_P
  //
  // This is independent of dt.
  //
  // Therefore:
  //
  //     CFL_P = dt * lambda_P
  //
  // and the transport-only explicit limit is:
  //
  //     dt <= 1 / lambda_P
  // ========================================================

  void computeTransportLossRate(field::CellField<double>& loss_rate) const {
    if (&loss_rate.mesh() != mesh_) {
      throw std::invalid_argument("loss-rate field belongs to another mesh");
    }

    loss_rate.fill(0.0);

    for (mesh::FaceId face = 0; face < mesh_->numFaces(); ++face) {

      const auto owner = mesh_->owner(face);

      const auto normal = mesh_->faceNormal(face);

      const double vn = (*normal_drift_velocity_)[face];

      // =================================================
      // Internal face
      // =================================================

      if (!mesh_->isBoundary(face)) {

        const auto neighbor = mesh_->neighbor(face);

        const auto delta =
            mesh_->cellCenter(neighbor) - mesh_->cellCenter(owner);

        const double distance = mesh::dot(delta, normal);

        if (distance <= 0.0) {
          throw std::runtime_error(
              "invalid owner-neighbor "
              "normal distance");
        }

        const double pe = vn * distance / diffusivity_;

        const double scale = diffusivity_ * mesh_->faceArea(face) / distance;

        // Owner self-loss coefficient:
        //
        //     D A / d * B(-Pe)

        loss_rate[owner] += scale * discretization::operators::bernoulli(-pe) /
                            mesh_->cellVolume(owner);

        // Neighbor's outward normal is -n_f.
        //
        // Its self-loss coefficient is:
        //
        //     D A / d * B(+Pe)

        loss_rate[neighbor] += scale *
                               discretization::operators::bernoulli(pe) /
                               mesh_->cellVolume(neighbor);

        continue;
      }

      // =================================================
      // Boundary face
      // =================================================

      if (wall_flux_.has_value() && wall_flux_->active(face)) {
        continue;
      }

      const auto boundary_id = mesh_->boundaryId(face);

      if (boundary_id == mesh::invalid_boundary) {

        throw std::runtime_error("boundary face has no BoundaryId");
      }

      if (!bc_->contains(boundary_id)) {
        throw std::runtime_error("species boundary condition missing");
      }

      // Current SG implementation supports
      // Dirichlet boundaries only.

      const auto* dirichlet =
          std::get_if<boundary::Dirichlet>(&bc_->at(boundary_id));

      if (dirichlet == nullptr) {
        throw std::runtime_error(
            "explicit SG transport currently "
            "requires Dirichlet boundaries");
      }

      const auto delta = mesh_->faceCenter(face) - mesh_->cellCenter(owner);

      const double distance = mesh::dot(delta, normal);

      if (distance <= 0.0) {
        throw std::runtime_error(
            "invalid owner-boundary "
            "normal distance");
      }

      const double pe = vn * distance / diffusivity_;

      const double scale = diffusivity_ * mesh_->faceArea(face) / distance;

      loss_rate[owner] += scale * discretization::operators::bernoulli(-pe) /
                          mesh_->cellVolume(owner);
    }
    if (wall_flux_.has_value()) {
      discretization::operators::computeLinearBoundaryLossRate(*wall_flux_,
                                                               wall_loss_rate_);
      // Validate every merged value before publishing any wall contribution.
      for (mesh::CellId cell = 0; cell < mesh_->numCells(); ++cell) {
        if (!std::isfinite(loss_rate[cell] + wall_loss_rate_[cell])) {
          throw std::overflow_error("transport loss rate overflowed");
        }
      }
      for (mesh::CellId cell = 0; cell < mesh_->numCells(); ++cell) {
        loss_rate[cell] += wall_loss_rate_[cell];
      }
    }
  }

  // ========================================================
  // Compute the complete Explicit-Euler increment WITHOUT
  // modifying density:
  //
  //     delta n
  //
  //       = dt * (S - div Gamma)
  //
  // This separation is important for multi-species atomicity:
  // every species increment can be computed successfully
  // before any species field is modified.
  // ========================================================

  void computeIncrement(const field::CellField<double>& density,
                        const field::CellField<double>& source, double dt,
                        field::CellField<double>& increment) {
    validateField(density);
    validateField(source);
    validateField(increment);

    if (dt <= 0.0) {
      throw std::invalid_argument("time step must be positive");
    }

    // ----------------------------------------------------
    // Gamma = drift-diffusion SG flux
    // ----------------------------------------------------

    computeNormalFlux(density, flux_);

    // ----------------------------------------------------
    // div(Gamma)
    // ----------------------------------------------------

    discretization::operators::divergence(flux_, divergence_);

    // ----------------------------------------------------
    // delta n = dt * (S - div Gamma)
    // ----------------------------------------------------

    for (mesh::CellId cell = 0; cell < mesh_->numCells(); ++cell) {

      increment[cell] = dt * (source[cell] - divergence_[cell]);
    }
  }

  // ========================================================
  // Convenience API for a single species.
  //
  // Multi-species solver should normally use
  // computeIncrement() first and commit all species together.
  // ========================================================

  void step(field::CellField<double>& density,
            const field::CellField<double>& source, double dt) {
    computeIncrement(density, source, dt, local_increment_);

    for (mesh::CellId cell = 0; cell < mesh_->numCells(); ++cell) {

      density[cell] += local_increment_[cell];
    }
  }

 private:
  void validateField(const field::CellField<double>& field) const {
    if (&field.mesh() != mesh_) {
      throw std::invalid_argument("field belongs to another mesh");
    }
  }

 private:
  const mesh::IMesh* mesh_;

  const field::FaceField<double>* normal_drift_velocity_;

  const boundary::BoundaryConditionSet* bc_;

  std::optional<discretization::operators::LinearBoundaryFluxView> wall_flux_;

  double diffusivity_;

  field::FaceField<double> flux_;

  field::CellField<double> divergence_;

  field::CellField<double> local_increment_;

  mutable field::CellField<double> wall_loss_rate_;
};

}  // namespace pemu::equation
