#include <pemu/equation/explicit_electron_energy_stepper.hpp>

#include <pemu/equation/adaptive_step_explicit_species_continuity_stepper.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>

namespace pemu::equation {

struct ExplicitElectronEnergyStepper::Impl {
  const mesh::IMesh* mesh;
  const field::FaceField<double>* electron_normal_drift_velocity;
  double energy_transport_factor;
  double energy_diffusivity;
  boundary::BoundaryConditionSet energy_boundary_conditions;
  std::optional<discretization::operators::LinearBoundaryFluxView> wall_flux;
  field::FaceField<double> energy_normal_drift_velocity;
  std::unique_ptr<AdaptiveStepExplicitSpeciesContinuityStepper> transport;
  field::CellField<double> transport_loss_rate;
  field::CellField<double> increment;

  /** @brief Validates closure parameters and returns factor * D_e. */
  [[nodiscard]] static double validatedEnergyDiffusivity(
      double electron_diffusivity, double factor) {
    if (!std::isfinite(electron_diffusivity) || electron_diffusivity <= 0.0) {
      throw std::invalid_argument(
          "electron diffusivity must be finite and positive");
    }
    if (!std::isfinite(factor) || factor <= 0.0) {
      throw std::invalid_argument(
          "energy transport factor must be finite and positive");
    }
    const double result = electron_diffusivity * factor;
    if (!std::isfinite(result)) {
      throw std::invalid_argument("electron energy diffusivity overflowed");
    }
    return result;
  }

  /** @brief Builds workspaces and binds the SG operator to their addresses. */
  Impl(const mesh::IMesh& mesh_in,
       const field::FaceField<double>& electron_velocity,
       double electron_diffusivity,
       boundary::BoundaryConditionSet boundary_conditions, double factor,
       std::optional<discretization::operators::LinearBoundaryFluxView>
           wall_flux_in)
      : mesh(&mesh_in),
        electron_normal_drift_velocity(&electron_velocity),
        energy_transport_factor(factor),
        energy_diffusivity(
            validatedEnergyDiffusivity(electron_diffusivity, factor)),
        energy_boundary_conditions(std::move(boundary_conditions)),
        wall_flux(std::move(wall_flux_in)),
        energy_normal_drift_velocity(mesh_in, 0.0),
        transport_loss_rate(mesh_in, 0.0),
        increment(mesh_in, 0.0) {
    if (&electron_velocity.mesh() != &mesh_in) {
      throw std::invalid_argument(
          "electron drift velocity belongs to another mesh");
    }
    if (wall_flux.has_value() && &wall_flux->mesh() != &mesh_in) {
      throw std::invalid_argument(
          "electron energy wall flux belongs to another mesh");
    }
    validateBoundaryConditions();
    if (wall_flux.has_value()) {
      transport =
          std::make_unique<AdaptiveStepExplicitSpeciesContinuityStepper>(
              mesh_in, energy_normal_drift_velocity, energy_diffusivity,
              energy_boundary_conditions, *wall_flux);
    } else {
      transport =
          std::make_unique<AdaptiveStepExplicitSpeciesContinuityStepper>(
              mesh_in, energy_normal_drift_velocity, energy_diffusivity,
              energy_boundary_conditions);
    }
  }

  /** @brief Requires finite non-negative Dirichlet energy on every boundary. */
  void validateBoundaryConditions() const {
    for (mesh::FaceId face = 0; face < mesh->numFaces(); ++face) {
      if (!mesh->isBoundary(face)) {
        continue;
      }
      if (wall_flux.has_value() && wall_flux->active(face)) {
        continue;
      }
      const auto boundary_id = mesh->boundaryId(face);
      if (boundary_id == mesh::invalid_boundary ||
          !energy_boundary_conditions.contains(boundary_id)) {
        throw std::invalid_argument(
            "electron energy boundary condition is missing");
      }
      const auto* dirichlet = std::get_if<boundary::Dirichlet>(
          &energy_boundary_conditions.at(boundary_id));
      if (dirichlet == nullptr || !std::isfinite(dirichlet->value) ||
          dirichlet->value < 0.0) {
        throw std::invalid_argument(
            "electron energy transport requires finite non-negative "
            "Dirichlet boundaries");
      }
    }
  }

  /** @brief Validates that a cell field belongs to this equation mesh. */
  void validateMesh(const field::CellField<double>& value,
                    std::string_view name) const {
    if (&value.mesh() != mesh) {
      throw std::invalid_argument(std::string{name} +
                                  " belongs to another mesh");
    }
  }

  /** @brief Validates a finite, non-negative electron energy state. */
  void validateEnergyDensity(
      const field::CellField<double>& energy_density) const {
    validateMesh(energy_density, "electron energy density");
    for (const double value : energy_density) {
      if (!std::isfinite(value) || value < 0.0) {
        throw std::invalid_argument(
            "electron energy density must be finite and non-negative");
      }
    }
  }

  /** @brief Validates a finite signed energy source field. */
  void validateSource(const field::CellField<double>& source) const {
    validateMesh(source, "electron energy source");
    for (const double value : source) {
      if (!std::isfinite(value)) {
        throw std::invalid_argument("electron energy source must be finite");
      }
    }
  }

  /** @brief Refreshes v_energy = factor * v_e from the current face field. */
  void updateEnergyDriftVelocity() {
    for (mesh::FaceId face = 0; face < mesh->numFaces(); ++face) {
      const double electron_velocity = (*electron_normal_drift_velocity)[face];
      if (!std::isfinite(electron_velocity)) {
        throw std::invalid_argument("electron drift velocity must be finite");
      }
      const double energy_velocity =
          energy_transport_factor * electron_velocity;
      if (!std::isfinite(energy_velocity)) {
        throw std::overflow_error("electron energy drift velocity overflowed");
      }
      energy_normal_drift_velocity[face] = energy_velocity;
    }
  }

  /** @brief Refreshes and validates all SG diagonal loss coefficients. */
  void updateTransportLossRate() {
    updateEnergyDriftVelocity();
    transport->computeTransportLossRate(transport_loss_rate);
    for (const double loss_rate : transport_loss_rate) {
      if (!std::isfinite(loss_rate) || loss_rate < 0.0) {
        throw std::runtime_error(
            "electron energy transport loss rate is invalid");
      }
    }
  }

  /** @brief Reduces the current loss-rate field to its global CFL limit. */
  [[nodiscard]] double transportLimitFromLossRate() const {
    double limit = std::numeric_limits<double>::infinity();
    for (const double loss_rate : transport_loss_rate) {
      if (loss_rate > 0.0) {
        limit = std::min(limit, 1.0 / loss_rate);
      }
    }
    return limit;
  }

  /** @brief Computes the current transport-only stability limit. */
  [[nodiscard]] double maxStableTransportTimeStep() {
    updateTransportLossRate();
    return transportLimitFromLossRate();
  }

  /** @brief Reduces current loss rates and signed sources to a positive dt. */
  [[nodiscard]] double positivityLimitFromLossRate(
      const field::CellField<double>& energy_density,
      const field::CellField<double>& source) const {
    double limit = std::numeric_limits<double>::infinity();
    for (mesh::CellId cell = 0; cell < mesh->numCells(); ++cell) {
      const double depletion =
          transport_loss_rate[cell] * energy_density[cell] +
          std::max(-source[cell], 0.0);
      if (depletion > 0.0) {
        limit = std::min(limit, energy_density[cell] / depletion);
      }
    }
    return limit;
  }

  /** @brief Computes a conservative positivity limit from loss and source. */
  [[nodiscard]] double maxPositiveTimeStep(
      const field::CellField<double>& energy_density,
      const field::CellField<double>& source) {
    validateEnergyDensity(energy_density);
    validateSource(source);
    updateTransportLossRate();
    return positivityLimitFromLossRate(energy_density, source);
  }

  /** @brief Computes and validates an increment in the private scratch field. */
  void computeIncrement(const field::CellField<double>& energy_density,
                        const field::CellField<double>& source, double dt) {
    validateEnergyDensity(energy_density);
    validateSource(source);
    if (!std::isfinite(dt) || dt <= 0.0) {
      throw std::invalid_argument("electron energy time step must be positive");
    }
    updateEnergyDriftVelocity();
    transport->computeIncrement(energy_density, source, dt, increment);
    for (const double value : increment) {
      if (!std::isfinite(value)) {
        throw std::runtime_error("electron energy increment is not finite");
      }
    }
  }
};

/** @copydoc ExplicitElectronEnergyStepper::ExplicitElectronEnergyStepper */
ExplicitElectronEnergyStepper::ExplicitElectronEnergyStepper(
    const mesh::IMesh& mesh,
    const field::FaceField<double>& electron_normal_drift_velocity,
    double electron_diffusivity,
    boundary::BoundaryConditionSet energy_boundary_conditions,
    double energy_transport_factor)
    : impl_(std::make_unique<Impl>(mesh, electron_normal_drift_velocity,
                                   electron_diffusivity,
                                   std::move(energy_boundary_conditions),
                                   energy_transport_factor, std::nullopt)) {}

/** @copydoc ExplicitElectronEnergyStepper::ExplicitElectronEnergyStepper */
ExplicitElectronEnergyStepper::ExplicitElectronEnergyStepper(
    const mesh::IMesh& mesh,
    const field::FaceField<double>& electron_normal_drift_velocity,
    double electron_diffusivity,
    boundary::BoundaryConditionSet energy_boundary_conditions,
    discretization::operators::LinearBoundaryFluxView wall_flux,
    double energy_transport_factor)
    : impl_(std::make_unique<Impl>(mesh, electron_normal_drift_velocity,
                                   electron_diffusivity,
                                   std::move(energy_boundary_conditions),
                                   energy_transport_factor, wall_flux)) {}

/** @copydoc ExplicitElectronEnergyStepper::~ExplicitElectronEnergyStepper */
ExplicitElectronEnergyStepper::~ExplicitElectronEnergyStepper() = default;

/** @copydoc ExplicitElectronEnergyStepper::ExplicitElectronEnergyStepper */
ExplicitElectronEnergyStepper::ExplicitElectronEnergyStepper(
    ExplicitElectronEnergyStepper&&) noexcept = default;

/** @copydoc ExplicitElectronEnergyStepper::operator= */
ExplicitElectronEnergyStepper& ExplicitElectronEnergyStepper::operator=(
    ExplicitElectronEnergyStepper&&) noexcept = default;

/** @copydoc ExplicitElectronEnergyStepper::mesh */
const mesh::IMesh& ExplicitElectronEnergyStepper::mesh() const noexcept {
  return *impl_->mesh;
}

/** @copydoc ExplicitElectronEnergyStepper::energyTransportFactor */
double ExplicitElectronEnergyStepper::energyTransportFactor() const noexcept {
  return impl_->energy_transport_factor;
}

/** @copydoc ExplicitElectronEnergyStepper::energyDiffusivity */
double ExplicitElectronEnergyStepper::energyDiffusivity() const noexcept {
  return impl_->energy_diffusivity;
}

/** @copydoc ExplicitElectronEnergyStepper::computeTransportLossRate */
void ExplicitElectronEnergyStepper::computeTransportLossRate(
    field::CellField<double>& loss_rate) {
  impl_->validateMesh(loss_rate, "electron energy loss rate");
  impl_->updateTransportLossRate();
  for (mesh::CellId cell = 0; cell < impl_->mesh->numCells(); ++cell) {
    loss_rate[cell] = impl_->transport_loss_rate[cell];
  }
}

/** @copydoc ExplicitElectronEnergyStepper::maxStableTransportTimeStep */
double ExplicitElectronEnergyStepper::maxStableTransportTimeStep() {
  return impl_->maxStableTransportTimeStep();
}

/** @copydoc ExplicitElectronEnergyStepper::maxPositiveTimeStep */
double ExplicitElectronEnergyStepper::maxPositiveTimeStep(
    const field::CellField<double>& energy_density,
    const field::CellField<double>& source) {
  return impl_->maxPositiveTimeStep(energy_density, source);
}

/** @copydoc ExplicitElectronEnergyStepper::computeIncrement */
void ExplicitElectronEnergyStepper::computeIncrement(
    const field::CellField<double>& energy_density,
    const field::CellField<double>& source, double dt,
    field::CellField<double>& increment) {
  impl_->validateMesh(increment, "electron energy increment");
  impl_->computeIncrement(energy_density, source, dt);
  for (mesh::CellId cell = 0; cell < impl_->mesh->numCells(); ++cell) {
    increment[cell] = impl_->increment[cell];
  }
}

/** @copydoc ExplicitElectronEnergyStepper::computeStableIncrement */
void ExplicitElectronEnergyStepper::computeStableIncrement(
    const field::CellField<double>& energy_density,
    const field::CellField<double>& source, double dt,
    field::CellField<double>& increment) {
  impl_->validateMesh(increment, "electron energy increment");
  impl_->validateEnergyDensity(energy_density);
  impl_->validateSource(source);
  if (!std::isfinite(dt) || dt <= 0.0) {
    throw std::invalid_argument("electron energy time step must be positive");
  }

  impl_->updateTransportLossRate();
  const double transport_limit = impl_->transportLimitFromLossRate();
  if (std::isfinite(transport_limit) && dt > transport_limit * (1.0 + 1e-12)) {
    throw std::runtime_error("electron energy transport CFL violated");
  }
  const double positivity_limit =
      impl_->positivityLimitFromLossRate(energy_density, source);
  if (std::isfinite(positivity_limit) &&
      dt > positivity_limit * (1.0 + 1e-12)) {
    throw std::runtime_error("electron energy positivity limit violated");
  }

  impl_->computeIncrement(energy_density, source, dt);
  for (mesh::CellId cell = 0; cell < impl_->mesh->numCells(); ++cell) {
    const double candidate = energy_density[cell] + impl_->increment[cell];
    if (!std::isfinite(candidate) || candidate < 0.0) {
      throw std::runtime_error(
          "electron energy update would violate non-negativity");
    }
    increment[cell] = impl_->increment[cell];
  }
}

/** @copydoc ExplicitElectronEnergyStepper::step */
void ExplicitElectronEnergyStepper::step(
    field::CellField<double>& energy_density,
    const field::CellField<double>& source, double dt) {
  computeStableIncrement(energy_density, source, dt, impl_->increment);
  for (mesh::CellId cell = 0; cell < impl_->mesh->numCells(); ++cell) {
    energy_density[cell] += impl_->increment[cell];
  }
}

}  // namespace pemu::equation
