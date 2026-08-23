#include <pemu/discretization/operators/drift_velocity.hpp>
#include <pemu/discretization/operators/electric_field.hpp>
#include <pemu/discretization/operators/linear_boundary_flux.hpp>
#include <pemu/equation/fixed_step_multi_species_drift_diffusion_stepper.hpp>
#include <pemu/physics/charge_density.hpp>

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>

namespace pemu::equation {

std::unique_ptr<linalg::ISolver>
FixedStepMultiSpeciesDriftDiffusionStepper::validateBackend(
    std::unique_ptr<linalg::ISolver> backend) {
  if (!backend) {
    throw std::invalid_argument(
        "Poisson linear solver backend must not be null");
  }
  return backend;
}

void FixedStepMultiSpeciesDriftDiffusionStepper::buildTransportSteppers() {
  transport_steppers_.resize(species_->size());

  for (std::size_t i = 0; i < species_->size(); ++i) {
    const physics::SpeciesId id{static_cast<std::uint32_t>(i)};
    const auto& properties = species_->at(id);
    properties.validate();

    if (!properties.isTransported()) {
      continue;
    }
    if (!properties.isCharged()) {
      throw std::invalid_argument(
          "current drift-diffusion solver supports charged species only");
    }

    if (wall_flux_assembler_.has_value()) {
      const discretization::operators::LinearBoundaryFluxView wall_flux(
          wall_flux_assembler_->activeFaces()[id],
          wall_flux_assembler_->particleLossVelocity()[id],
          wall_flux_assembler_->particleInwardFlux()[id]);
      transport_steppers_[i] =
          std::make_unique<FixedStepExplicitSpeciesContinuityStepper>(
              *mesh_, drift_velocity_[id], properties.diffusivity, dt_,
              species_bc_[i], wall_flux);
    } else {
      transport_steppers_[i] =
          std::make_unique<FixedStepExplicitSpeciesContinuityStepper>(
              *mesh_, drift_velocity_[id], properties.diffusivity, dt_,
              species_bc_[i]);
    }
  }
}

void FixedStepMultiSpeciesDriftDiffusionStepper::validateFields(
    const physics::SpeciesCellFields& fields) const {
  if (&fields.mesh() != mesh_) {
    throw std::invalid_argument("species fields belong to another mesh");
  }
  if (fields.size() != species_->size()) {
    throw std::invalid_argument(
        "species field count does not match SpeciesSet");
  }
}

linalg::SolverResult
FixedStepMultiSpeciesDriftDiffusionStepper::updateElectrostatics(
    const physics::SpeciesCellFields& density) {
  validateFields(density);
  electrostatics_ready_ = false;
  charge_density_.fill(0.0);

  for (std::size_t i = 0; i < species_->size(); ++i) {
    const physics::SpeciesId id{static_cast<std::uint32_t>(i)};
    const auto& properties = species_->at(id);
    if (properties.charge != 0.0) {
      physics::addSpeciesChargeDensity(density[id], properties.charge,
                                       charge_density_);
    }
  }

  const auto result = poisson_solver_.solve(potential_);
  if (!result.success()) {
    return result;
  }

  discretization::operators::electricFieldNormal(
      potential_, permittivity_, potential_bc_, electric_field_normal_);

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

  if (wall_flux_assembler_.has_value()) {
    wall_flux_assembler_->evaluate(density);
  }

  electrostatics_ready_ = true;
  return result;
}

void FixedStepMultiSpeciesDriftDiffusionStepper::computeParticleFluxNormal(
    const physics::SpeciesCellFields& density, physics::SpeciesId id,
    field::FaceField<double>& normal_flux) const {
  validateFields(density);
  if (!electrostatics_ready_) {
    throw std::logic_error(
        "electrostatics must be updated before computing particle flux");
  }
  const auto index = static_cast<std::size_t>(id.value);
  if (index >= transport_steppers_.size() || !transport_steppers_[index]) {
    throw std::invalid_argument("particle flux requires a transported species");
  }
  transport_steppers_[index]->computeNormalFlux(density[id], normal_flux);
}

void FixedStepMultiSpeciesDriftDiffusionStepper::advanceTransport(
    physics::SpeciesCellFields& density,
    const physics::SpeciesCellFields& source) {
  validateFields(density);
  validateFields(source);

  if (!electrostatics_ready_) {
    throw std::logic_error(
        "electrostatics must be updated before advancing transport");
  }

  // Preserve atomicity: validate every fixed-step CFL and candidate state
  // before changing any species field.
  for (const auto& stepper : transport_steppers_) {
    if (stepper && stepper->maxTransportCfl() > 1.0) {
      throw std::runtime_error("multi-species explicit transport CFL violated");
    }
  }

  increments_.fill(0.0);
  for (std::size_t i = 0; i < species_->size(); ++i) {
    if (!transport_steppers_[i]) {
      continue;
    }
    const physics::SpeciesId id{static_cast<std::uint32_t>(i)};
    transport_steppers_[i]->computeIncrement(density[id], source[id],
                                             increments_[id]);
  }

  for (std::size_t i = 0; i < species_->size(); ++i) {
    if (!transport_steppers_[i]) {
      continue;
    }
    const physics::SpeciesId id{static_cast<std::uint32_t>(i)};
    for (mesh::CellId cell = 0; cell < mesh_->numCells(); ++cell) {
      const double candidate = density[id][cell] + increments_[id][cell];
      const double scale = std::max(1.0, std::abs(density[id][cell]));
      if (!std::isfinite(candidate) || candidate < -1e-12 * scale) {
        throw std::runtime_error(
            "multi-species update would produce negative density");
      }
      if (candidate < 0.0) {
        increments_[id][cell] = -density[id][cell];
      }
    }
  }

  for (std::size_t i = 0; i < species_->size(); ++i) {
    if (!transport_steppers_[i]) {
      continue;
    }
    const physics::SpeciesId id{static_cast<std::uint32_t>(i)};
    for (mesh::CellId cell = 0; cell < mesh_->numCells(); ++cell) {
      density[id][cell] += increments_[id][cell];
    }
  }

  electrostatics_ready_ = false;
}

FixedStepMultiSpeciesDriftDiffusionStepper::
    FixedStepMultiSpeciesDriftDiffusionStepper(
        const mesh::IMesh& mesh, const physics::SpeciesSet& species,
        double permittivity, double dt,
        boundary::BoundaryConditionSet potential_bc,
        std::vector<boundary::BoundaryConditionSet> species_bc,
        std::unique_ptr<linalg::ISolver> poisson_backend,
        field::PlasmaFieldMetadata field_metadata,
        std::optional<PureNeumannOptions> pure_neumann_options,
        physics::wall::WallBoundarySet wall_boundaries)
    : mesh_(&mesh),
      species_(&species),
      permittivity_(permittivity),
      dt_(dt),
      potential_bc_(std::move(potential_bc)),
      species_bc_(std::move(species_bc)),
      field_metadata_(std::move(field_metadata)),
      charge_density_(mesh, 0.0, field_metadata_.charge_density),
      potential_(mesh, 0.0, field_metadata_.electric_potential),
      electric_field_normal_(mesh, 0.0, field_metadata_.electric_field),
      drift_velocity_(mesh, species.size(), 0.0,
                      field_metadata_.drift_velocity),
      increments_(mesh, species.size(), 0.0, field_metadata_.number_density),
      poisson_solver_(discretization::PoissonFvm(mesh, charge_density_,
                                                 permittivity, potential_bc_),
                      validateBackend(std::move(poisson_backend)),
                      std::move(pure_neumann_options)) {
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
        "boundary condition count does not match species count");
  }
  if (!wall_boundaries.empty()) {
    wall_flux_assembler_.emplace(mesh, species, std::move(wall_boundaries));
  }
  buildTransportSteppers();
}

linalg::SolverResult FixedStepMultiSpeciesDriftDiffusionStepper::step(
    physics::SpeciesCellFields& density,
    const physics::SpeciesCellFields& source) {
  const auto result = updateElectrostatics(density);
  if (!result.success()) {
    return result;
  }
  advanceTransport(density, source);
  return result;
}

}  // namespace pemu::equation
