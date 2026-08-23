#pragma once

#include <pemu/boundary/boundary_condition_set.hpp>
#include <pemu/equation/fixed_step_explicit_species_continuity_stepper.hpp>
#include <pemu/equation/poisson_solver.hpp>
#include <pemu/field/plasma_field_metadata.hpp>
#include <pemu/mesh/i_mesh.hpp>
#include <pemu/physics/species.hpp>
#include <pemu/physics/wall/flux_assembler.hpp>

#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>

namespace pemu::equation {

class FixedStepMultiSpeciesDriftDiffusionStepper {
 public:
  /**
   * @brief Creates a fixed-step multi-species drift-diffusion stepper.
   * @param pure_neumann_options Gauge and compatibility tolerances required
   * when every potential boundary is Neumann.
   */
  FixedStepMultiSpeciesDriftDiffusionStepper(
      const mesh::IMesh& mesh, const physics::SpeciesSet& species,
      double permittivity, double dt,
      boundary::BoundaryConditionSet potential_bc,
      std::vector<boundary::BoundaryConditionSet> species_bc,
      std::unique_ptr<linalg::ISolver> poisson_backend,
      field::PlasmaFieldMetadata field_metadata = {},
      std::optional<PureNeumannOptions> pure_neumann_options = std::nullopt,
      physics::wall::WallBoundarySet wall_boundaries = {});

  FixedStepMultiSpeciesDriftDiffusionStepper(
      const FixedStepMultiSpeciesDriftDiffusionStepper&) = delete;
  FixedStepMultiSpeciesDriftDiffusionStepper& operator=(
      const FixedStepMultiSpeciesDriftDiffusionStepper&) = delete;

  // Poisson and wall views retain addresses of this object's member fields.
  // Moving the aggregate would leave those non-owning references dangling.
  FixedStepMultiSpeciesDriftDiffusionStepper(
      FixedStepMultiSpeciesDriftDiffusionStepper&&) = delete;
  FixedStepMultiSpeciesDriftDiffusionStepper& operator=(
      FixedStepMultiSpeciesDriftDiffusionStepper&&) = delete;

  [[nodiscard]] linalg::SolverResult updateElectrostatics(
      const physics::SpeciesCellFields& density);

  [[nodiscard]] linalg::SolverResult step(
      physics::SpeciesCellFields& density,
      const physics::SpeciesCellFields& source);

  [[nodiscard]] const field::CellField<double>& chargeDensity() const noexcept {
    return charge_density_;
  }

  [[nodiscard]] const field::CellField<double>& potential() const noexcept {
    return potential_;
  }

  [[nodiscard]] const field::FaceField<double>& electricFieldNormal()
      const noexcept {
    return electric_field_normal_;
  }

  [[nodiscard]] const field::FaceField<double>& driftVelocityNormal(
      physics::SpeciesId id) const {
    auto _ = species_->at(id);
    return drift_velocity_[id];
  }

  /** @brief Computes one species' current SG particle flux without advancing. */
  void computeParticleFluxNormal(const physics::SpeciesCellFields& density,
                                 physics::SpeciesId id,
                                 field::FaceField<double>& normal_flux) const;

  [[nodiscard]] double timeStep() const noexcept { return dt_; }

  [[nodiscard]] const field::PlasmaFieldMetadata& fieldMetadata()
      const noexcept {
    return field_metadata_;
  }

  [[nodiscard]] const physics::SpeciesSet& species() const noexcept {
    return *species_;
  }

  [[nodiscard]] const mesh::IMesh& mesh() const noexcept { return *mesh_; }

  /** @brief Reports whether plasma-wall particle flux is configured. */
  [[nodiscard]] bool hasWallFluxAssembler() const noexcept {
    return wall_flux_assembler_.has_value();
  }

  /** @brief Returns current wall fields for electron-energy coupling. */
  [[nodiscard]] const physics::wall::WallFluxAssembler& wallFluxAssembler()
      const {
    if (!wall_flux_assembler_.has_value()) {
      throw std::logic_error("plasma-wall flux assembly is not configured");
    }
    return *wall_flux_assembler_;
  }

  [[nodiscard]] double transportCfl(physics::SpeciesId id) const {
    auto _ = species_->at(id);
    const auto index = static_cast<std::size_t>(id.value);
    if (!transport_steppers_[index]) {
      return 0.0;
    }
    return transport_steppers_[index]->maxTransportCfl();
  }

  void advanceTransport(physics::SpeciesCellFields& density,
                        const physics::SpeciesCellFields& source);

 private:
  static std::unique_ptr<linalg::ISolver> validateBackend(
      std::unique_ptr<linalg::ISolver> backend);

  /** @brief Builds the private per-species steppers and binds wall views. */
  void buildTransportSteppers();

  void validateFields(const physics::SpeciesCellFields& density) const;

  const mesh::IMesh* mesh_;
  const physics::SpeciesSet* species_;
  double permittivity_;
  double dt_;
  boundary::BoundaryConditionSet potential_bc_;
  std::vector<boundary::BoundaryConditionSet> species_bc_;
  field::PlasmaFieldMetadata field_metadata_;
  field::CellField<double> charge_density_;
  field::CellField<double> potential_;
  field::FaceField<double> electric_field_normal_;
  physics::SpeciesFaceFields drift_velocity_;
  std::optional<physics::wall::WallFluxAssembler> wall_flux_assembler_;
  physics::SpeciesCellFields increments_;
  equation::PoissonSolver poisson_solver_;
  std::vector<std::unique_ptr<FixedStepExplicitSpeciesContinuityStepper>>
      transport_steppers_;
  bool electrostatics_ready_{false};
};

}  // namespace pemu::equation
