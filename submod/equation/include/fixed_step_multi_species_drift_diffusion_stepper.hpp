#pragma once

#include <pemu/boundary/boundary_condition_set.hpp>
#include <pemu/equation/fixed_step_explicit_species_continuity_stepper.hpp>
#include <pemu/equation/poisson_solver.hpp>
#include <pemu/field/plasma_field_metadata.hpp>
#include <pemu/mesh/i_mesh.hpp>
#include <pemu/physics/species.hpp>

#include <memory>
#include <vector>

namespace pemu::equation {

class FixedStepMultiSpeciesDriftDiffusionStepper {
 public:
  FixedStepMultiSpeciesDriftDiffusionStepper(
      const mesh::IMesh& mesh, const physics::SpeciesSet& species,
      double permittivity, double dt,
      boundary::BoundaryConditionSet potential_bc,
      std::vector<boundary::BoundaryConditionSet> species_bc,
      std::unique_ptr<linalg::ISolver> poisson_backend,
      field::PlasmaFieldMetadata field_metadata = {});

  void buildTransportSteppers();

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

  [[nodiscard]] double timeStep() const noexcept { return dt_; }

  [[nodiscard]] const field::PlasmaFieldMetadata& fieldMetadata()
      const noexcept {
    return field_metadata_;
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
  equation::PoissonSolver poisson_solver_;
  std::vector<std::unique_ptr<FixedStepExplicitSpeciesContinuityStepper>>
      transport_steppers_;
  bool electrostatics_ready_{false};
};

}  // namespace pemu::equation
