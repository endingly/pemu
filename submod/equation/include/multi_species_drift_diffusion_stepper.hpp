#pragma once

// #include <algorithm>
// #include <memory>
#include <pemu/boundary/boundary_condition_set.hpp>
// #include <pemu/discretization/operators/drift_velocity.hpp>
// #include <pemu/discretization/operators/electric_field.hpp>
// #include <pemu/discretization/poisson_fvm.hpp>
#include <pemu/equation/explicit_species_continuity_stepper.hpp>
#include <pemu/equation/poisson_solver.hpp>
// #include <pemu/field/cell_field.hpp>
// #include <pemu/field/face_field.hpp>
// #include <pemu/linalg/i_solver.hpp>
#include <pemu/mesh/i_mesh.hpp>
// #include <pemu/physics/charge_density.hpp>
// #include <pemu/physics/charged_species_transport.hpp>
// #include <pemu/physics/reaction.hpp>
#include <pemu/physics/species.hpp>
// #include <stdexcept>
// #include <utility>

namespace pemu::equation {

class MultiSpeciesDriftDiffusionStepper {
 public:
  MultiSpeciesDriftDiffusionStepper(
      const mesh::IMesh& mesh, const physics::SpeciesSet& species,
      double permittivity, double dt,
      boundary::BoundaryConditionSet potential_bc,
      std::vector<boundary::BoundaryConditionSet> species_bc,
      std::unique_ptr<linalg::ISolver> poisson_backend);

  void buildTransportSteppers();

  linalg::SolverResult updateElectrostatics(
      const physics::SpeciesCellFields& density);

  linalg::SolverResult step(physics::SpeciesCellFields& density,
                            const physics::SpeciesCellFields& source);

  [[nodiscard]]
  const field::CellField<double>& chargeDensity() const noexcept {
    return charge_density_;
  }

  [[nodiscard]]
  const field::CellField<double>& potential() const noexcept {
    return potential_;
  }

  [[nodiscard]]
  const field::FaceField<double>& electricFieldNormal() const noexcept {
    return electric_field_normal_;
  }

  [[nodiscard]]
  const field::FaceField<double>& driftVelocityNormal(
      physics::SpeciesId id) const {
    auto _ = species_->at(id);
    return drift_velocity_[id];
  }

  [[nodiscard]]
  double transportCfl(physics::SpeciesId id) const {
    auto _ = species_->at(id);
    const auto index = static_cast<std::size_t>(id.value);
    if (!transport_steppers_[index]) {
      return 0.0;
    }
    return transport_steppers_[index]->maxTransportCfl();
  }

 private:
  static std::unique_ptr<linalg::ISolver> validateBackend(
      std::unique_ptr<linalg::ISolver> backend);

  void validateFields(const physics::SpeciesCellFields& density) const;

 private:
  const mesh::IMesh* mesh_;
  const physics::SpeciesSet* species_;
  double permittivity_;
  double dt_;
  boundary::BoundaryConditionSet potential_bc_;
  std::vector<boundary::BoundaryConditionSet> species_bc_;
  field::CellField<double> charge_density_;
  field::CellField<double> potential_;
  field::FaceField<double> electric_field_normal_;
  physics::SpeciesFaceFields drift_velocity_;
  equation::PoissonSolver poisson_solver_;
  std::vector<std::unique_ptr<ExplicitSpeciesContinuityStepper>>
      transport_steppers_;
};

};  // namespace pemu::equation