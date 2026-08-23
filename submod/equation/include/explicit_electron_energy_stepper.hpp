#pragma once

#include <pemu/boundary/boundary_condition_set.hpp>
#include <pemu/discretization/operators/linear_boundary_flux.hpp>
#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>
#include <pemu/mesh/i_mesh.hpp>

#include <memory>

namespace pemu::equation {

/**
 * @brief Advances electron energy density with explicit SG transport.
 *
 * The transported state is w_e = n_e * mean_energy. The caller supplies the
 * complete non-transport source, including field work and collisional terms.
 * This class deliberately has no Simulation dependency.
 */
class ExplicitElectronEnergyStepper {
 public:
  /** @brief Maxwellian closure coefficient for energy mobility/diffusivity. */
  static constexpr double maxwellian_transport_factor = 5.0 / 3.0;

  /**
   * @brief Creates an electron energy transport operator.
   * @param electron_normal_drift_velocity Current electron particle drift
   * velocity on faces; subsequent field updates are observed by the stepper.
   * @param electron_diffusivity Electron particle diffusivity.
   * @param energy_boundary_conditions Dirichlet energy density or prescribed
   * outward diffusive Neumann flux on each ordinary boundary.
   * @param energy_transport_factor Ratio mu_energy/mu_e and D_energy/D_e.
   */
  ExplicitElectronEnergyStepper(
      const mesh::IMesh& mesh,
      const field::FaceField<double>& electron_normal_drift_velocity,
      double electron_diffusivity,
      boundary::BoundaryConditionSet energy_boundary_conditions,
      double energy_transport_factor = maxwellian_transport_factor);

  /** @brief Creates energy transport with linear wall flux on selected faces. */
  ExplicitElectronEnergyStepper(
      const mesh::IMesh& mesh,
      const field::FaceField<double>& electron_normal_drift_velocity,
      double electron_diffusivity,
      boundary::BoundaryConditionSet energy_boundary_conditions,
      discretization::operators::LinearBoundaryFluxView wall_flux,
      double energy_transport_factor = maxwellian_transport_factor);

  /** @brief Releases the private transport workspace. */
  ~ExplicitElectronEnergyStepper();

  ExplicitElectronEnergyStepper(const ExplicitElectronEnergyStepper&) = delete;
  ExplicitElectronEnergyStepper& operator=(
      const ExplicitElectronEnergyStepper&) = delete;

  /** @brief Transfers ownership of the private transport workspace. */
  ExplicitElectronEnergyStepper(ExplicitElectronEnergyStepper&&) noexcept;

  /** @brief Replaces this stepper with another private transport workspace. */
  ExplicitElectronEnergyStepper& operator=(
      ExplicitElectronEnergyStepper&&) noexcept;

  /** @brief Returns the mesh on which the energy equation is defined. */
  [[nodiscard]] const mesh::IMesh& mesh() const noexcept;

  /** @brief Returns the closure multiplier applied to particle transport. */
  [[nodiscard]] double energyTransportFactor() const noexcept;

  /** @brief Returns D_energy = factor * D_e. */
  [[nodiscard]] double energyDiffusivity() const noexcept;

  /**
   * @brief Computes the diagonal loss rate of the explicit SG operator.
   *
   * The resulting inverse-time field is independent of the selected dt.
   */
  void computeTransportLossRate(field::CellField<double>& loss_rate);

  /** @brief Returns the largest positivity-preserving transport-only dt. */
  [[nodiscard]] double maxStableTransportTimeStep();

  /**
   * @brief Returns a conservative state-dependent positivity limit.
   *
   * Positive sources and incoming boundary fluxes are deliberately not
   * credited. Adaptive integration must combine this bound with the separate
   * transport-only limit, including when the current cell energy is zero.
   */
  [[nodiscard]] double maxPositiveTimeStep(
      const field::CellField<double>& energy_density,
      const field::CellField<double>& source);

  /**
   * @brief Computes dt * (source - div(Gamma_energy)) without changing state.
   */
  void computeIncrement(const field::CellField<double>& energy_density,
                        const field::CellField<double>& source, double dt,
                        field::CellField<double>& increment);

  /**
   * @brief Computes an increment only when dt satisfies all explicit limits.
   *
   * In addition to transport CFL and conservative source-depletion checks,
   * this method verifies the complete candidate state. It is intended for a
   * coupled workflow that must validate every equation before committing any
   * of them.
   */
  void computeStableIncrement(const field::CellField<double>& energy_density,
                              const field::CellField<double>& source, double dt,
                              field::CellField<double>& increment);

  /**
   * @brief Advances one stable explicit step and commits it atomically.
   *
   * The method rejects transport CFL violations and any update that would
   * produce non-finite or negative electron energy density.
   */
  void step(field::CellField<double>& energy_density,
            const field::CellField<double>& source, double dt);

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace pemu::equation
