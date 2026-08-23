#pragma once

#include <pemu/boundary/boundary_condition_set.hpp>

#include <pemu/discretization/operators/drift_velocity.hpp>
#include <pemu/discretization/operators/electric_field.hpp>
#include <pemu/discretization/poisson_fvm.hpp>

#include <pemu/equation/adaptive_step_explicit_species_continuity_stepper.hpp>
#include <pemu/equation/poisson_solver.hpp>

#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>
#include <pemu/field/plasma_field_metadata.hpp>

#include <pemu/linalg/i_solver.hpp>

#include <pemu/physics/charge_density.hpp>

#include <pemu/physics/species.hpp>
#include <pemu/physics/wall/flux_assembler.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pemu::equation::detail {

class ExplicitMultiSpeciesDriftDiffusionOperator {
 public:
  /**
   * @brief Creates the explicit multi-species transport operator.
   * @param pure_neumann_options Gauge and compatibility tolerances required
   * when every potential boundary is Neumann.
   */
  ExplicitMultiSpeciesDriftDiffusionOperator(
      const mesh::IMesh& mesh, const physics::SpeciesSet& species,
      double permittivity,

      boundary::BoundaryConditionSet potential_bc,

      std::vector<boundary::BoundaryConditionSet> species_bc,

      std::unique_ptr<linalg::ISolver> poisson_backend,

      field::PlasmaFieldMetadata field_metadata = {},

      std::optional<PureNeumannOptions> pure_neumann_options = std::nullopt,

      physics::wall::WallBoundarySet wall_boundaries = {})

      : mesh_(&mesh),

        species_(&species),

        permittivity_(permittivity),

        potential_bc_(std::move(potential_bc)),

        species_bc_(std::move(species_bc)),

        field_metadata_(std::move(field_metadata)),

        charge_density_(mesh, 0.0, field_metadata_.charge_density),

        potential_(mesh, 0.0, field_metadata_.electric_potential),

        electric_field_normal_(mesh, 0.0, field_metadata_.electric_field),

        drift_velocity_(mesh, species.size(), 0.0,
                        field_metadata_.drift_velocity),

        transport_loss_rate_(mesh, species.size(), 0.0,
                             field_metadata_.inverse_time),

        prescribed_boundary_sink_rate_(mesh, species.size(), 0.0,
                                       field_metadata_.number_density_source),

        increments_(mesh, species.size(), 0.0, field_metadata_.number_density),

        poisson_solver_(discretization::PoissonFvm(mesh, charge_density_,
                                                   permittivity, potential_bc_),
                        validateBackend(std::move(poisson_backend)),
                        std::move(pure_neumann_options)) {
    if (permittivity_ <= 0.0) {
      throw std::invalid_argument("permittivity must be positive");
    }

    if (species_->size() == 0) {
      throw std::invalid_argument("species set must not be empty");
    }

    if (species_bc_.size() != species_->size()) {

      throw std::invalid_argument(
          "species BC count does not "
          "match SpeciesSet");
    }

    if (!wall_boundaries.empty()) {
      wall_flux_assembler_.emplace(*mesh_, *species_,
                                   std::move(wall_boundaries));
    }

    buildTransportSteppers();
  }

  ExplicitMultiSpeciesDriftDiffusionOperator(
      const ExplicitMultiSpeciesDriftDiffusionOperator&) = delete;

  ExplicitMultiSpeciesDriftDiffusionOperator& operator=(
      const ExplicitMultiSpeciesDriftDiffusionOperator&) = delete;

  //
  // Important:
  //
  // PoissonFvm and species steppers contain pointers to
  // members owned by THIS object.
  //
  // Therefore moving this object would invalidate those
  // internal object-address assumptions.
  //

  ExplicitMultiSpeciesDriftDiffusionOperator(
      ExplicitMultiSpeciesDriftDiffusionOperator&&) = delete;

  ExplicitMultiSpeciesDriftDiffusionOperator& operator=(
      ExplicitMultiSpeciesDriftDiffusionOperator&&) = delete;

  // ========================================================
  // n_s -> rho -> phi -> E_n -> v_n,s
  // ========================================================

  [[nodiscard]]
  linalg::SolverResult updateElectrostatics(
      const physics::SpeciesCellFields& density) {
    validateFields(density);

    electrostatics_ready_ = false;

    // ----------------------------------------------------
    // rho = sum_s q_s n_s
    // ----------------------------------------------------

    charge_density_.fill(0.0);

    forEachSpecies([&](physics::SpeciesId id) {
      const auto& properties = species_->at(id);

      if (properties.charge == 0.0) {

        return;
      }

      physics::addSpeciesChargeDensity(density[id], properties.charge,
                                       charge_density_);
    });

    // ----------------------------------------------------
    // rho -> phi
    // ----------------------------------------------------

    const auto result = poisson_solver_.solve(potential_);

    if (!result.success()) {
      return result;
    }

    // ----------------------------------------------------
    // phi -> E_n
    // ----------------------------------------------------

    discretization::operators::electricFieldNormal(
        potential_, permittivity_, potential_bc_, electric_field_normal_);

    // ----------------------------------------------------
    // E_n -> v_n,s
    // ----------------------------------------------------

    forEachSpecies([&](physics::SpeciesId id) {
      const auto& properties = species_->at(id);

      if (!properties.isTransported()) {

        drift_velocity_[id].fill(0.0);

        return;
      }

      discretization::operators::driftVelocityNormal(
          electric_field_normal_, properties.mobility, properties.polarity(),
          drift_velocity_[id]);
    });

    if (wall_flux_assembler_.has_value()) {
      wall_flux_assembler_->evaluate(density);
    }

    electrostatics_ready_ = true;

    return result;
  }

  // ========================================================
  // Compute transport lambda_s,P.
  // ========================================================

  void updateTransportLossRates() {
    requireElectrostaticsReady();

    transport_loss_rate_.fill(0.0);
    prescribed_boundary_sink_rate_.fill(0.0);

    forEachSpecies([&](physics::SpeciesId id) {
      auto* stepper = transportStepper(id);

      if (stepper == nullptr) {
        return;
      }

      stepper->computeTransportLossRate(transport_loss_rate_[id]);
      stepper->computePrescribedBoundarySinkRate(
          prescribed_boundary_sink_rate_[id]);
    });
  }

  // ========================================================
  // Transport-only explicit stability limit:
  //
  //     dt <= min 1/lambda_s,P
  // ========================================================

  [[nodiscard]]
  double maxStableTransportTimeStep() {
    requireElectrostaticsReady();

    updateTransportLossRates();

    double limit = std::numeric_limits<double>::infinity();

    forEachSpecies([&](physics::SpeciesId id) {
      if (transportStepper(id) == nullptr) {

        return;
      }

      for (mesh::CellId cell = 0; cell < mesh_->numCells(); ++cell) {

        const double lambda = transport_loss_rate_[id][cell];

        if (lambda > 0.0) {

          limit = std::min(limit, 1.0 / lambda);
        }
      }
    });

    return limit;
  }

  // ========================================================
  // Combined transport, source and prescribed-boundary positivity limit:
  //
  //     dt <=
  //
  //           n
  //     ----------------
  //     lambda*n + source_sink + boundary_sink
  //
  // where:
  //
  //     source_sink = max(-S, 0)
  //     boundary_sink = sum max(q_neumann, 0) A / V
  //
  // Positive source and incoming Neumann flux do not constrain positivity.
  // ========================================================

  [[nodiscard]]
  double positivityTimeStepLimit(const physics::SpeciesCellFields& density,
                                 const physics::SpeciesCellFields& source) {
    validateFields(density);

    validateFields(source);

    requireElectrostaticsReady();

    updateTransportLossRates();

    double limit = std::numeric_limits<double>::infinity();

    forEachSpecies([&](physics::SpeciesId id) {
      if (transportStepper(id) == nullptr) {

        return;
      }

      for (mesh::CellId cell = 0; cell < mesh_->numCells(); ++cell) {

        const double n = density[id][cell];

        if (n < 0.0) {

          throw std::runtime_error("negative species density");
        }

        const double lambda = transport_loss_rate_[id][cell];

        const double sink = std::max(-source[id][cell], 0.0);

        const double depletion_rate =
            lambda * n + sink + prescribed_boundary_sink_rate_[id][cell];

        if (depletion_rate == 0.0) {

          continue;
        }

        if (n == 0.0) {

          //
          // A sink is trying to remove
          // particles from an empty cell.
          //

          limit = 0.0;
          return;
        }

        limit = std::min(limit, n / depletion_rate);
      }
    });

    return limit;
  }

  // ========================================================
  // Advance using externally selected dt.
  //
  // This function does NOT select dt.
  //
  // Therefore fixed-step and adaptive-step wrappers can
  // share exactly the same physics implementation.
  // ========================================================

  void advanceTransport(physics::SpeciesCellFields& density,
                        const physics::SpeciesCellFields& source, double dt) {
    validateFields(density);

    validateFields(source);

    requireElectrostaticsReady();

    if (dt <= 0.0) {
      throw std::invalid_argument("time step must be positive");
    }

    // ----------------------------------------------------
    // Positivity pre-check.
    // ----------------------------------------------------

    const double limit = positivityTimeStepLimit(density, source);

    if (limit == 0.0) {
      throw std::runtime_error(
          "no positive timestep satisfies "
          "species positivity");
    }

    if (std::isfinite(limit) && dt > limit * (1.0 + 1e-12)) {

      throw std::runtime_error(
          "explicit species timestep "
          "violates positivity limit");
    }

    // ----------------------------------------------------
    // Phase 1:
    //
    // Compute EVERY increment.
    //
    // Do not modify any density yet.
    // ----------------------------------------------------

    increments_.fill(0.0);

    forEachSpecies([&](physics::SpeciesId id) {
      auto* stepper = transportStepper(id);

      if (stepper == nullptr) {
        return;
      }

      stepper->computeIncrement(density[id], source[id], dt, increments_[id]);
    });

    // ----------------------------------------------------
    // Phase 2:
    //
    // Validate EVERY candidate state before commit.
    // ----------------------------------------------------

    forEachSpecies([&](physics::SpeciesId id) {
      if (transportStepper(id) == nullptr) {

        return;
      }

      for (mesh::CellId cell = 0; cell < mesh_->numCells(); ++cell) {

        const double candidate = density[id][cell] + increments_[id][cell];

        const double scale = std::max(1.0, std::abs(density[id][cell]));

        if (candidate < -1e-12 * scale) {

          throw std::runtime_error(
              "species update would "
              "produce negative density");
        }
      }
    });

    // ----------------------------------------------------
    // Phase 3:
    //
    // Commit ALL transported species.
    // ----------------------------------------------------

    forEachSpecies([&](physics::SpeciesId id) {
      if (transportStepper(id) == nullptr) {

        return;
      }

      for (mesh::CellId cell = 0; cell < mesh_->numCells(); ++cell) {

        density[id][cell] += increments_[id][cell];
      }
    });

    //
    // density is now n^(k+1),
    // while electrostatics corresponded to n^k.
    //

    electrostatics_ready_ = false;
  }

  // ========================================================
  // Accessors
  // ========================================================

  [[nodiscard]]
  bool electrostaticsReady() const noexcept {
    return electrostatics_ready_;
  }

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

  /** @brief Computes one species' current SG particle flux without advancing. */
  void computeParticleFluxNormal(const physics::SpeciesCellFields& density,
                                 physics::SpeciesId id,
                                 field::FaceField<double>& normal_flux) const {
    validateFields(density);
    requireElectrostaticsReady();
    const auto* stepper = transportStepper(id);
    if (stepper == nullptr) {
      throw std::invalid_argument(
          "particle flux requires a transported species");
    }
    stepper->computeNormalFlux(density[id], normal_flux);
  }

  [[nodiscard]]
  const physics::SpeciesCellFields& transportLossRate() const noexcept {
    return transport_loss_rate_;
  }

  [[nodiscard]]
  const physics::SpeciesSet& species() const noexcept {
    return *species_;
  }

  [[nodiscard]]
  const mesh::IMesh& mesh() const noexcept {
    return *mesh_;
  }

  [[nodiscard]]
  const field::PlasmaFieldMetadata& fieldMetadata() const noexcept {
    return field_metadata_;
  }

  /** @brief Reports whether plasma-wall flux assembly is configured. */
  [[nodiscard]] bool hasWallFluxAssembler() const noexcept {
    return wall_flux_assembler_.has_value();
  }

  /** @brief Returns current wall flux fields for coupled energy transport. */
  [[nodiscard]] const physics::wall::WallFluxAssembler& wallFluxAssembler()
      const {
    if (!wall_flux_assembler_.has_value()) {
      throw std::logic_error("plasma-wall flux assembly is not configured");
    }
    return *wall_flux_assembler_;
  }

 private:
  static std::unique_ptr<linalg::ISolver> validateBackend(
      std::unique_ptr<linalg::ISolver> backend) {
    if (!backend) {
      throw std::invalid_argument(
          "Poisson solver backend "
          "must not be null");
    }

    return backend;
  }

  void buildTransportSteppers() {
    transport_steppers_.resize(species_->size());

    forEachSpecies([&](physics::SpeciesId id) {
      const auto& properties = species_->at(id);

      properties.validate();

      if (!properties.isTransported()) {
        return;
      }

      if (!properties.isCharged()) {

        throw std::invalid_argument(
            "current drift-diffusion "
            "core transports charged "
            "species only");
      }

      if (wall_flux_assembler_.has_value()) {
        const discretization::operators::LinearBoundaryFluxView wall_flux(
            wall_flux_assembler_->activeFaces()[id],
            wall_flux_assembler_->particleLossVelocity()[id],
            wall_flux_assembler_->particleInwardFlux()[id]);
        transport_steppers_[id.value] =
            std::make_unique<AdaptiveStepExplicitSpeciesContinuityStepper>(
                *mesh_, drift_velocity_[id], properties.diffusivity,
                species_bc_[id.value], wall_flux);
      } else {
        transport_steppers_[id.value] =
            std::make_unique<AdaptiveStepExplicitSpeciesContinuityStepper>(
                *mesh_, drift_velocity_[id], properties.diffusivity,
                species_bc_[id.value]);
      }
    });
  }

  [[nodiscard]]
  AdaptiveStepExplicitSpeciesContinuityStepper* transportStepper(
      physics::SpeciesId id) {
    return transport_steppers_.at(id.value).get();
  }

  [[nodiscard]]
  const AdaptiveStepExplicitSpeciesContinuityStepper* transportStepper(
      physics::SpeciesId id) const {
    return transport_steppers_.at(id.value).get();
  }

  void validateFields(const physics::SpeciesCellFields& fields) const {
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

  void requireElectrostaticsReady() const {
    if (!electrostatics_ready_) {

      throw std::logic_error(
          "electrostatics must be "
          "updated first");
    }
  }

  template <typename Function>
  void forEachSpecies(Function&& function) {
    for (std::size_t i = 0; i < species_->size(); ++i) {

      function(physics::SpeciesId{static_cast<std::uint32_t>(i)});
    }
  }

 private:
  const mesh::IMesh* mesh_;

  const physics::SpeciesSet* species_;

  double permittivity_;

  // --------------------------------------------------------
  // These objects own the BC storage.
  //
  // Internal solvers contain pointers/references to them.
  // Therefore this operator object must not be moved.
  // --------------------------------------------------------

  boundary::BoundaryConditionSet potential_bc_;

  std::vector<boundary::BoundaryConditionSet> species_bc_;

  field::PlasmaFieldMetadata field_metadata_;

  // --------------------------------------------------------
  // Electrostatic workspace
  // --------------------------------------------------------

  field::CellField<double> charge_density_;

  field::CellField<double> potential_;

  field::FaceField<double> electric_field_normal_;

  // --------------------------------------------------------
  // Per-species workspace
  // --------------------------------------------------------

  physics::SpeciesFaceFields drift_velocity_;

  std::optional<physics::wall::WallFluxAssembler> wall_flux_assembler_;

  physics::SpeciesCellFields transport_loss_rate_;

  physics::SpeciesCellFields prescribed_boundary_sink_rate_;

  physics::SpeciesCellFields increments_;

  // --------------------------------------------------------
  // Poisson
  // --------------------------------------------------------

  PoissonSolver poisson_solver_;

  // --------------------------------------------------------
  // One stepper for every transported species.
  //
  // nullptr => species is Immobile.
  // --------------------------------------------------------

  std::vector<std::unique_ptr<AdaptiveStepExplicitSpeciesContinuityStepper>>
      transport_steppers_;

  bool electrostatics_ready_{false};
};

}  // namespace pemu::equation::detail
