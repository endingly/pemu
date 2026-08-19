#pragma once

#include <pemu/boundary/boundary_condition_set.hpp>

#include <pemu/discretization/operators/drift_velocity.hpp>
#include <pemu/discretization/operators/electric_field.hpp>
#include <pemu/discretization/poisson_fvm.hpp>

#include <pemu/equation/explicit_species_continuity_stepper.hpp>
#include <pemu/equation/poisson_solver.hpp>

#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>

#include <pemu/linalg/i_solver.hpp>

#include <pemu/mesh/i_mesh.hpp>

#include <pemu/physics/charge_density.hpp>
#include <pemu/physics/charged_species_transport.hpp>

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <utility>

namespace pemu::equation {

/**
 * @brief Explicitly coupled electrostatic drift-diffusion stepper.
 *
 * The timestep is
 *
 *   n_e^k, n_i^k
 *        |
 *        v
 *      rho^k
 *        |
 *        v
 *     Poisson
 *        |
 *        v
 *      phi^k
 *        |
 *        v
 *       E^k
 *      /   \
 *     v     v
 *   v_e^k  v_i^k
 *     |      |
 *     v      v
 *    SG     SG
 *     |      |
 *     v      v
 * n_e^{k+1} n_i^{k+1}
 *
 * Both species are advanced using the same frozen electric
 * field E^k.
 *
 * Note:
 *
 * After step(), density fields contain time level k+1, while
 * chargeDensity(), potential(), electricFieldNormal(), and the
 * drift-velocity fields still correspond to time level k.
 *
 * Call updateElectrostatics() again if synchronized k+1
 * electrostatic fields are required.
 */
class ElectrostaticDriftDiffusionStepper {
 public:
  ElectrostaticDriftDiffusionStepper(
      const mesh::IMesh& mesh, double permittivity, double dt,
      physics::ChargedSpeciesTransport electron,
      physics::ChargedSpeciesTransport ion,
      const boundary::BoundaryConditionSet& potential_bc,
      const boundary::BoundaryConditionSet& electron_bc,
      const boundary::BoundaryConditionSet& ion_bc,
      std::unique_ptr<linalg::ISolver> poisson_backend)
      : mesh_(&mesh),
        permittivity_(permittivity),

        electron_(electron),
        ion_(ion),

        potential_bc_(&potential_bc),
        electron_bc_(&electron_bc),
        ion_bc_(&ion_bc),

        charge_density_(mesh, 0.0),
        potential_(mesh, 0.0),

        electric_field_normal_(mesh, 0.0),

        electron_drift_velocity_normal_(mesh, 0.0),

        ion_drift_velocity_normal_(mesh, 0.0),

        poisson_solver_(discretization::PoissonFvm(mesh, charge_density_,
                                                   permittivity, potential_bc),
                        validateBackend(std::move(poisson_backend))),

        electron_stepper_(mesh, electron_drift_velocity_normal_,
                          electron.diffusivity, dt, electron_bc),

        ion_stepper_(mesh, ion_drift_velocity_normal_, ion.diffusivity, dt,
                     ion_bc) {
    validateConfiguration(dt);
  }

  ElectrostaticDriftDiffusionStepper(
      const ElectrostaticDriftDiffusionStepper&) = delete;

  ElectrostaticDriftDiffusionStepper& operator=(
      const ElectrostaticDriftDiffusionStepper&) = delete;

  ElectrostaticDriftDiffusionStepper(ElectrostaticDriftDiffusionStepper&&) =
      delete;

  ElectrostaticDriftDiffusionStepper& operator=(
      ElectrostaticDriftDiffusionStepper&&) = delete;

  // ========================================================
  // Electrostatic update
  // ========================================================

  /**
     * @brief Update electrostatic state from the current
     * electron and ion number densities.
     *
     * Computes:
     *
     *   rho = q_e n_e + q_i n_i
     *
     * then solves:
     *
     *   -div(epsilon grad(phi)) = rho
     *
     * followed by:
     *
     *   E_n = -grad(phi) dot n
     *
     * and:
     *
     *   v_{n,s} = sign(q_s) mu_s E_n
     */
  [[nodiscard]]
  linalg::SolverResult updateElectrostatics(
      const field::CellField<double>& electron_density,
      const field::CellField<double>& ion_density) {
    validateDensityFields(electron_density, ion_density);

    // ----------------------------------------------------
    // rho = q_e n_e + q_i n_i
    // ----------------------------------------------------

    charge_density_.fill(0.0);

    physics::addSpeciesChargeDensity(electron_density, electron_.charge,
                                     charge_density_);

    physics::addSpeciesChargeDensity(ion_density, ion_.charge, charge_density_);

    // ----------------------------------------------------
    // -div(epsilon grad(phi)) = rho
    // ----------------------------------------------------

    const auto poisson_result = poisson_solver_.solve(potential_);

    if (!poisson_result.success()) {
      return poisson_result;
    }

    // ----------------------------------------------------
    // phi -> E_n
    //
    //     E_n = -grad(phi) dot n
    // ----------------------------------------------------

    discretization::operators::electricFieldNormal(
        potential_, permittivity_, *potential_bc_, electric_field_normal_);

    // ----------------------------------------------------
    // Electron:
    //
    // q_e < 0
    //
    //     v_e = -mu_e E
    // ----------------------------------------------------

    discretization::operators::driftVelocityNormal(
        electric_field_normal_, electron_.mobility, electron_.polarity(),
        electron_drift_velocity_normal_);

    // ----------------------------------------------------
    // Positive ion:
    //
    // q_i > 0
    //
    //     v_i = +mu_i E
    // ----------------------------------------------------

    discretization::operators::driftVelocityNormal(
        electric_field_normal_, ion_.mobility, ion_.polarity(),
        ion_drift_velocity_normal_);

    return poisson_result;
  }

  // ========================================================
  // Coupled timestep
  // ========================================================

  /**
     * @brief Advance electron and ion densities by one explicit
     * electrostatic drift-diffusion timestep.
     *
     * Sequence:
     *
     *   1. rho^k = q_e n_e^k + q_i n_i^k
     *   2. solve phi^k
     *   3. compute E^k
     *   4. compute v_e^k and v_i^k
     *   5. verify both explicit SG CFL constraints
     *   6. advance electron density
     *   7. advance ion density
     *
     * Both species see the same frozen electrostatic state.
     */
  [[nodiscard]]
  linalg::SolverResult step(field::CellField<double>& electron_density,
                            const field::CellField<double>& electron_source,

                            field::CellField<double>& ion_density,
                            const field::CellField<double>& ion_source) {
    validateDensityFields(electron_density, ion_density);

    validateSourceFields(electron_source, ion_source);

    // ----------------------------------------------------
    // n^k
    //   ->
    // rho^k
    //   ->
    // phi^k
    //   ->
    // E^k
    //   ->
    // v_s^k
    // ----------------------------------------------------

    const auto poisson_result =
        updateElectrostatics(electron_density, ion_density);

    if (!poisson_result.success()) {
      return poisson_result;
    }

    // ----------------------------------------------------
    // IMPORTANT:
    //
    // Check BOTH species before modifying either density.
    //
    // This gives step-level atomicity with respect to CFL
    // rejection:
    //
    // either both species can advance,
    // or neither density is modified.
    // ----------------------------------------------------

    const double electron_cfl = electron_stepper_.maxTransportCfl();

    const double ion_cfl = ion_stepper_.maxTransportCfl();

    if (electron_cfl > 1.0) {
      throw std::runtime_error(
          "electron explicit SG transport "
          "CFL condition violated");
    }

    if (ion_cfl > 1.0) {
      throw std::runtime_error(
          "ion explicit SG transport "
          "CFL condition violated");
    }

    // ----------------------------------------------------
    // Both species now use the same frozen E^k.
    // ----------------------------------------------------

    electron_stepper_.step(electron_density, electron_source);

    ion_stepper_.step(ion_density, ion_source);

    return poisson_result;
  }

  // ========================================================
  // Diagnostics / state access
  // ========================================================

  [[nodiscard]]
  const mesh::IMesh& mesh() const noexcept {
    return *mesh_;
  }

  [[nodiscard]]
  double permittivity() const noexcept {
    return permittivity_;
  }

  [[nodiscard]]
  const physics::ChargedSpeciesTransport& electronProperties() const noexcept {
    return electron_;
  }

  [[nodiscard]]
  const physics::ChargedSpeciesTransport& ionProperties() const noexcept {
    return ion_;
  }

  /**
     * @brief Charge density corresponding to the most recent
     * electrostatic update.
     */
  [[nodiscard]]
  const field::CellField<double>& chargeDensity() const noexcept {
    return charge_density_;
  }

  /**
     * @brief Potential corresponding to the most recent
     * electrostatic update.
     */
  [[nodiscard]]
  const field::CellField<double>& potential() const noexcept {
    return potential_;
  }

  /**
     * @brief Face-normal electric field:
     *
     *     E_n = E dot n
     */
  [[nodiscard]]
  const field::FaceField<double>& electricFieldNormal() const noexcept {
    return electric_field_normal_;
  }

  /**
     * @brief Electron face-normal drift velocity.
     */
  [[nodiscard]]
  const field::FaceField<double>& electronDriftVelocityNormal() const noexcept {
    return electron_drift_velocity_normal_;
  }

  /**
     * @brief Ion face-normal drift velocity.
     */
  [[nodiscard]]
  const field::FaceField<double>& ionDriftVelocityNormal() const noexcept {
    return ion_drift_velocity_normal_;
  }

  /**
     * @brief Current generalized explicit SG CFL for electrons.
     *
     * This requires updateElectrostatics() to have been called
     * if the electric field has changed.
     */
  [[nodiscard]]
  double electronTransportCfl() const {
    return electron_stepper_.maxTransportCfl();
  }

  /**
     * @brief Current generalized explicit SG CFL for ions.
     */
  [[nodiscard]]
  double ionTransportCfl() const {
    return ion_stepper_.maxTransportCfl();
  }

  // ========================================================
  // Explicit synchronization helper
  // ========================================================

  /**
     * @brief Recompute rho, phi, E and drift velocities using
     * the supplied density state.
     *
     * Useful after step() when the caller wants all diagnostic
     * fields synchronized with n^{k+1}.
     */
  [[nodiscard]]
  linalg::SolverResult synchronizeElectrostatics(
      const field::CellField<double>& electron_density,
      const field::CellField<double>& ion_density) {
    return updateElectrostatics(electron_density, ion_density);
  }

 private:
  // ========================================================
  // Backend validation
  // ========================================================

  [[nodiscard]]
  static std::unique_ptr<linalg::ISolver> validateBackend(
      std::unique_ptr<linalg::ISolver> backend) {
    if (!backend) {
      throw std::invalid_argument(
          "Poisson linear solver backend "
          "must not be null");
    }

    return backend;
  }

  // ========================================================
  // Configuration validation
  // ========================================================

  void validateConfiguration(double dt) const {
    if (permittivity_ <= 0.0) {
      throw std::invalid_argument("permittivity must be positive");
    }

    if (dt <= 0.0) {
      throw std::invalid_argument("time step must be positive");
    }

    electron_.validate();
    ion_.validate();

    if (electron_.charge >= 0.0) {
      throw std::invalid_argument(
          "electron species must have "
          "negative charge");
    }

    if (ion_.charge <= 0.0) {
      throw std::invalid_argument(
          "ion species must have "
          "positive charge");
    }
  }

  // ========================================================
  // Field validation
  // ========================================================

  void validateDensityFields(
      const field::CellField<double>& electron_density,
      const field::CellField<double>& ion_density) const {
    if (&electron_density.mesh() != mesh_) {

      throw std::invalid_argument(
          "electron density belongs "
          "to another mesh");
    }

    if (&ion_density.mesh() != mesh_) {

      throw std::invalid_argument(
          "ion density belongs "
          "to another mesh");
    }
  }

  void validateSourceFields(const field::CellField<double>& electron_source,
                            const field::CellField<double>& ion_source) const {
    if (&electron_source.mesh() != mesh_) {

      throw std::invalid_argument(
          "electron source belongs "
          "to another mesh");
    }

    if (&ion_source.mesh() != mesh_) {

      throw std::invalid_argument(
          "ion source belongs "
          "to another mesh");
    }
  }

 private:
  // ========================================================
  // Configuration
  //
  // These must outlive all dependent internal objects.
  // ========================================================

  const mesh::IMesh* mesh_;

  double permittivity_;

  physics::ChargedSpeciesTransport electron_;

  physics::ChargedSpeciesTransport ion_;

  // --------------------------------------------------------
  // Boundary-condition objects are non-owning.
  //
  // Their lifetime must exceed this stepper's lifetime.
  //
  // IMPORTANT:
  //
  // Boundary CONDITION TYPES must not change after Poisson
  // factorization has been initialized, because switching
  // Dirichlet <-> Neumann changes the Poisson matrix.
  //
  // Changing only boundary values is compatible with the
  // current Poisson RHS reassembly design.
  // --------------------------------------------------------

  const boundary::BoundaryConditionSet* potential_bc_;

  const boundary::BoundaryConditionSet* electron_bc_;

  const boundary::BoundaryConditionSet* ion_bc_;

  // ========================================================
  // Electrostatic workspace
  // ========================================================

  field::CellField<double> charge_density_;

  field::CellField<double> potential_;

  field::FaceField<double> electric_field_normal_;

  // ========================================================
  // Species drift workspace
  //
  // These objects are intentionally stable in memory because
  // ExplicitSpeciesContinuityStepper stores references to
  // them.
  // ========================================================

  field::FaceField<double> electron_drift_velocity_normal_;

  field::FaceField<double> ion_drift_velocity_normal_;

  // ========================================================
  // Solvers
  //
  // Declaration order is important:
  //
  // charge_density_ must exist before PoissonFvm stores its
  // reference.
  //
  // drift velocity fields must exist before the continuity
  // steppers store their references.
  // ========================================================

  PoissonSolver poisson_solver_;

  ExplicitSpeciesContinuityStepper electron_stepper_;

  ExplicitSpeciesContinuityStepper ion_stepper_;
};

}  // namespace pemu::equation