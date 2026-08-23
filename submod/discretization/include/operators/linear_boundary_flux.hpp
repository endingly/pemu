#pragma once

#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>
#include <pemu/field/utils.hpp>

#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace pemu::discretization::operators {

/**
 * @brief Non-owning face law Gamma.n = v_loss*u_owner - Gamma_in.
 *
 * The byte mask selects boundary faces overridden by this law. Referenced
 * fields may be refreshed between evaluations without rebuilding the view.
 */
class LinearBoundaryFluxView {
 public:
  /** @brief Binds a stable mask, loss velocity, and inward-flux field. */
  LinearBoundaryFluxView(const field::FaceField<std::uint8_t>& active,
                         const field::FaceField<double>& loss_velocity,
                         const field::FaceField<double>& inward_flux)
      : active_(&active),
        loss_velocity_(&loss_velocity),
        inward_flux_(&inward_flux) {
    field::ensureSameMesh(active, loss_velocity);
    field::ensureSameMesh(active, inward_flux);
  }

  /** @brief Returns the mesh shared by all referenced face fields. */
  [[nodiscard]] const mesh::IMesh& mesh() const noexcept {
    return active_->mesh();
  }

  /** @brief Reports whether one face uses the linear boundary law. */
  [[nodiscard]] bool active(mesh::FaceId face) const noexcept {
    return (*active_)[face] != std::uint8_t{0};
  }

  /** @brief Returns the non-negative outward loss velocity on one face. */
  [[nodiscard]] double lossVelocity(mesh::FaceId face) const noexcept {
    return (*loss_velocity_)[face];
  }

  /** @brief Returns the non-negative prescribed inward flux on one face. */
  [[nodiscard]] double inwardFlux(mesh::FaceId face) const noexcept {
    return (*inward_flux_)[face];
  }

  /** @brief Validates all active faces and current coefficients atomically. */
  template <field::CellFieldLike StateField>
  void validate(const StateField& state) const {
    field::ensureSameMesh(*active_, state);
    for (mesh::FaceId face = 0; face < mesh().numFaces(); ++face) {
      if (!active(face)) {
        continue;
      }
      if (!mesh().isBoundary(face)) {
        throw std::invalid_argument(
            "linear boundary flux mask contains an internal face");
      }
      const double owner_state = state[mesh().owner(face)];
      const double loss_velocity = lossVelocity(face);
      const double inward_flux = inwardFlux(face);
      if (!std::isfinite(owner_state) || owner_state < 0.0 ||
          !std::isfinite(loss_velocity) || loss_velocity < 0.0 ||
          !std::isfinite(inward_flux) || inward_flux < 0.0) {
        throw std::invalid_argument(
            "linear boundary flux requires finite non-negative inputs");
      }
      const double candidate = loss_velocity * owner_state - inward_flux;
      if (!std::isfinite(candidate)) {
        throw std::overflow_error("linear boundary flux overflowed");
      }
    }
  }

  /** @brief Evaluates one active face after validate() has succeeded. */
  template <field::CellFieldLike StateField>
  [[nodiscard]] double normalFlux(const StateField& state,
                                  mesh::FaceId face) const noexcept {
    return lossVelocity(face) * state[mesh().owner(face)] - inwardFlux(face);
  }

 private:
  const field::FaceField<std::uint8_t>* active_;
  const field::FaceField<double>* loss_velocity_;
  const field::FaceField<double>* inward_flux_;
};

/**
 * @brief Overwrites a scratch field with sum(v_loss*A/V) for each cell.
 *
 * The destination is working storage and may be partially changed when an
 * exception is reported. Callers that publish a loss rate must merge this
 * scratch field only after the complete boundary sum succeeds.
 */
inline void computeLinearBoundaryLossRate(
    const LinearBoundaryFluxView& boundary_flux,
    field::CellField<double>& loss_rate) {
  if (&boundary_flux.mesh() != &loss_rate.mesh()) {
    throw std::invalid_argument(
        "linear boundary loss rate belongs to another mesh");
  }
  loss_rate.fill(0.0);
  for (mesh::FaceId face = 0; face < boundary_flux.mesh().numFaces(); ++face) {
    if (!boundary_flux.active(face)) {
      continue;
    }
    const auto& mesh = boundary_flux.mesh();
    if (!mesh.isBoundary(face) ||
        !std::isfinite(boundary_flux.lossVelocity(face)) ||
        boundary_flux.lossVelocity(face) < 0.0 ||
        !std::isfinite(mesh.faceArea(face)) || mesh.faceArea(face) <= 0.0) {
      throw std::invalid_argument("invalid linear boundary loss coefficient");
    }
    const auto owner = mesh.owner(face);
    if (!std::isfinite(mesh.cellVolume(owner)) ||
        mesh.cellVolume(owner) <= 0.0) {
      throw std::invalid_argument("invalid wall-owner cell volume");
    }
    const double contribution = boundary_flux.lossVelocity(face) *
                                mesh.faceArea(face) / mesh.cellVolume(owner);
    const double candidate = loss_rate[owner] + contribution;
    if (!std::isfinite(contribution) || !std::isfinite(candidate)) {
      throw std::overflow_error("linear boundary loss rate overflowed");
    }
    loss_rate[owner] = candidate;
  }
}

}  // namespace pemu::discretization::operators
