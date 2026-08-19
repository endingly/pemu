#include <pemu/discretization/operators/drift_velocity.hpp>
#include <pemu/discretization/operators/electric_field.hpp>
#include <pemu/equation/multi_species_drift_diffusion_stepper.hpp>
#include <pemu/physics/charge_density.hpp>

namespace pemu::equation {

std::unique_ptr<linalg::ISolver>
MultiSpeciesDriftDiffusionStepper::validateBackend(
    std::unique_ptr<linalg::ISolver> backend) {
  if (!backend) {
    throw std::invalid_argument(
        "Poisson linear solver backend "
        "must not be null");
  }
  return backend;
}

void MultiSpeciesDriftDiffusionStepper::buildTransportSteppers() {
  transport_steppers_.resize(species_->size());

  for (std::size_t i = 0; i < species_->size(); ++i) {

    const physics::SpeciesId id{static_cast<std::uint32_t>(i)};

    const auto& properties = species_->at(id);

    properties.validate();

    if (!properties.isTransported()) {
      continue;
    }

    //
    // Current electric drift model:
    //
    //     v = sign(q) mu E
    //
    // therefore a transported species using
    // electric drift must be charged.
    //
    if (!properties.isCharged()) {

      //
      // We could later support neutral diffusion
      // separately, but current SG drift-diffusion
      // path is intended for charged species.
      //
      throw std::invalid_argument(
          "current drift-diffusion solver "
          "supports charged species only");
    }

    transport_steppers_[i] = std::make_unique<ExplicitSpeciesContinuityStepper>(
        *mesh_, drift_velocity_[id], properties.diffusivity, dt_,
        species_bc_[i]);
  }
}

void MultiSpeciesDriftDiffusionStepper::validateFields(
    const physics::SpeciesCellFields& fields) const {
  if (&fields.mesh() != mesh_) {
    throw std::invalid_argument(
        "species fields belong "
        "to another mesh");
  }

  if (fields.size() != species_->size()) {

    throw std::invalid_argument(
        "species field count does "
        "not match SpeciesSet");
  }
}

linalg::SolverResult MultiSpeciesDriftDiffusionStepper::updateElectrostatics(
    const physics::SpeciesCellFields& density) {
  validateFields(density);

  // ========================================================
  // rho = sum_s q_s n_s
  // ========================================================

  charge_density_.fill(0.0);

  for (std::size_t i = 0; i < species_->size(); ++i) {

    const physics::SpeciesId id{static_cast<std::uint32_t>(i)};

    const auto& properties = species_->at(id);

    if (properties.charge == 0.0) {
      continue;
    }

    physics::addSpeciesChargeDensity(density[id], properties.charge,
                                     charge_density_);
  }

  // ========================================================
  // rho -> phi
  // ========================================================

  const auto result = poisson_solver_.solve(potential_);

  if (!result.success()) {
    return result;
  }

  // ========================================================
  // phi -> E_n
  // ========================================================

  discretization::operators::electricFieldNormal(
      potential_, permittivity_, potential_bc_, electric_field_normal_);

  // ========================================================
  // E_n -> v_n,s
  // ========================================================

  for (std::size_t i = 0; i < species_->size(); ++i) {

    const physics::SpeciesId id{static_cast<std::uint32_t>(i)};

    const auto& properties = species_->at(id);

    if (!properties.isTransported()) {

      drift_velocity_[id].fill(0.0);

      continue;
    }

    discretization::operators::driftVelocityNormal(
        electric_field_normal_, properties.mobility, properties.polarity(),
        drift_velocity_[id]);
  }

  return result;
}

MultiSpeciesDriftDiffusionStepper::MultiSpeciesDriftDiffusionStepper(
    const mesh::IMesh& mesh, const physics::SpeciesSet& species,
    double permittivity, double dt, boundary::BoundaryConditionSet potential_bc,
    std::vector<boundary::BoundaryConditionSet> species_bc,
    std::unique_ptr<linalg::ISolver> poisson_backend)
    : mesh_(&mesh),
      species_(&species),
      permittivity_(permittivity),
      dt_(dt),
      potential_bc_(std::move(potential_bc)),
      species_bc_(std::move(species_bc)),
      charge_density_(mesh, 0.0),
      potential_(mesh, 0.0),
      electric_field_normal_(mesh, 0.0),
      drift_velocity_(mesh, species.size(), 0.0),
      poisson_solver_(discretization::PoissonFvm(mesh, charge_density_,
                                                 permittivity, potential_bc_),
                      validateBackend(std::move(poisson_backend))) {
  if (permittivity <= 0.0) {
    throw std::invalid_argument("permittivity must be positive");
  }

  if (dt <= 0.0) {
    throw std::invalid_argument("time step must be positive");
  }

  if (species.size() == 0) {
    throw std::invalid_argument("species set must not be empty");
  }

  if (species_bc_.size() != species.size()) {

    throw std::invalid_argument(
        "boundary condition count "
        "does not match species count");
  }

  buildTransportSteppers();
}

linalg::SolverResult MultiSpeciesDriftDiffusionStepper::step(
    physics::SpeciesCellFields& density,
    const physics::SpeciesCellFields& source) {
  validateFields(density);
  validateFields(source);

  // ========================================================
  // n^k
  //
  // -> rho^k
  // -> phi^k
  // -> E^k
  // -> v_s^k
  // ========================================================

  const auto result = updateElectrostatics(density);

  if (!result.success()) {
    return result;
  }

  // ========================================================
  // PRECHECK ALL transported species.
  //
  // No density may be modified until all CFL checks pass.
  // ========================================================

  for (std::size_t i = 0; i < species_->size(); ++i) {

    if (!transport_steppers_[i]) {
      continue;
    }

    if (transport_steppers_[i]->maxTransportCfl() > 1.0) {

      throw std::runtime_error(
          "multi-species explicit "
          "transport CFL violated");
    }
  }

  // ========================================================
  // Advance all transported species using the same E^k.
  // ========================================================

  for (std::size_t i = 0; i < species_->size(); ++i) {

    if (!transport_steppers_[i]) {
      continue;
    }

    const physics::SpeciesId id{static_cast<std::uint32_t>(i)};

    transport_steppers_[i]->step(density[id], source[id]);
  }

  return result;
}

};  // namespace pemu::equation