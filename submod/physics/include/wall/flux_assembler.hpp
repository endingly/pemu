#pragma once

#include <pemu/physics/wall/types.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace pemu::physics::wall {

/**
 * @brief Pre-resolves plasma-wall definitions into reusable face fields.
 *
 * For one wall face and species s, the assembled outward particle flux is
 *
 *   Gamma_s.n = v_loss,s*n_s - Gamma_in,s.
 *
 * Secondary emission contributes to Gamma_in. Energy uses the analogous law
 * v_energy,s*w_s - Q_in,s. Construction performs all allocation; evaluate()
 * only refreshes density-dependent emission fields.
 */
class WallFluxAssembler {
 public:
  /** @brief Validates and binds wall definitions to boundary faces. */
  WallFluxAssembler(const mesh::IMesh& mesh, const SpeciesSet& species,
                    const WallBoundarySet& boundaries)
      : mesh_(&mesh),
        species_(&species),
        active_(mesh, species.size(), std::uint8_t{0}),
        particle_loss_velocity_(mesh, species.size(), 0.0),
        energy_loss_velocity_(mesh, species.size(), 0.0),
        particle_inward_flux_(mesh, species.size(), 0.0),
        energy_inward_flux_(mesh, species.size(), 0.0) {
    resolveBoundaries(boundaries);
    candidate_particle_inward_flux_.resize(resolved_targets_.size());
    candidate_energy_inward_flux_.resize(resolved_targets_.size());
  }

  /**
   * @brief Refreshes secondary particle and energy influx from current density.
   *
   * The previous published influx remains unchanged if validation or
  * arithmetic fails.
   */
  void evaluate(const SpeciesCellFields& density) {
    // A failed refresh preserves diagnostic fields but invalidates the
    // complete-flux query, which must never mix new density with old emission.
    evaluated_ = false;
    validateDensityLayout(density);
    validateIncidentDensities(density);
    std::fill(candidate_particle_inward_flux_.begin(),
              candidate_particle_inward_flux_.end(), 0.0);
    std::fill(candidate_energy_inward_flux_.begin(),
              candidate_energy_inward_flux_.end(), 0.0);

    for (const auto& emission : resolved_emissions_) {
      const auto owner = mesh_->owner(emission.face);
      const double incident_flux = checkedProduct(
          particle_loss_velocity_[emission.incident_species][emission.face],
          density[emission.incident_species][owner],
          "incident wall particle flux overflowed");
      const double emitted_flux =
          checkedProduct(emission.yield, incident_flux,
                         "secondary particle emission overflowed");
      const double emitted_energy_flux =
          checkedProduct(emission.emitted_mean_energy, emitted_flux,
                         "secondary energy emission overflowed");

      checkedAccumulate(candidate_particle_inward_flux_[emission.target_index],
                        emitted_flux, "secondary particle influx overflowed");
      checkedAccumulate(candidate_energy_inward_flux_[emission.target_index],
                        emitted_energy_flux,
                        "secondary energy influx overflowed");
    }

    for (std::size_t index = 0; index < resolved_targets_.size(); ++index) {
      const auto& target = resolved_targets_[index];
      particle_inward_flux_[target.species][target.face] =
          candidate_particle_inward_flux_[index];
      energy_inward_flux_[target.species][target.face] =
          candidate_energy_inward_flux_[index];
    }
    evaluated_ = true;
  }

  /** @brief Returns which species faces use the wall-flux law. */
  [[nodiscard]] const SpeciesWallFaceMasks& activeFaces() const noexcept {
    return active_;
  }

  /** @brief Returns particle loss velocities multiplying owner densities. */
  [[nodiscard]] const SpeciesFaceFields& particleLossVelocity() const noexcept {
    return particle_loss_velocity_;
  }

  /** @brief Returns energy loss velocities multiplying owner energy density. */
  [[nodiscard]] const SpeciesFaceFields& energyLossVelocity() const noexcept {
    return energy_loss_velocity_;
  }

  /** @brief Returns assembled inward secondary particle fluxes. */
  [[nodiscard]] const SpeciesFaceFields& particleInwardFlux() const noexcept {
    return particle_inward_flux_;
  }

  /** @brief Returns assembled inward secondary energy fluxes. */
  [[nodiscard]] const SpeciesFaceFields& energyInwardFlux() const noexcept {
    return energy_inward_flux_;
  }

  /** @brief Evaluates one complete outward particle wall flux. */
  [[nodiscard]] double particleNormalFlux(const SpeciesCellFields& density,
                                          SpeciesId species,
                                          mesh::FaceId face) const {
    requireEvaluated();
    validateWallFace(species, face);
    if (&density.mesh() != mesh_ || density.size() != species_->size()) {
      throw std::invalid_argument("wall density layout differs");
    }
    const double owner_density = density[species][mesh_->owner(face)];
    if (!std::isfinite(owner_density) || owner_density < 0.0) {
      throw std::invalid_argument(
          "wall density must be finite and non-negative");
    }
    const double loss =
        checkedProduct(particle_loss_velocity_[species][face], owner_density,
                       "wall particle loss overflowed");
    return loss - particle_inward_flux_[species][face];
  }

  /** @brief Evaluates one complete outward energy wall flux. */
  [[nodiscard]] double energyNormalFlux(
      const field::CellField<double>& energy_density, SpeciesId species,
      mesh::FaceId face) const {
    requireEvaluated();
    validateWallFace(species, face);
    if (&energy_density.mesh() != mesh_) {
      throw std::invalid_argument(
          "wall energy density belongs to another mesh");
    }
    const double owner_energy = energy_density[mesh_->owner(face)];
    if (!std::isfinite(owner_energy) || owner_energy < 0.0) {
      throw std::invalid_argument(
          "wall energy density must be finite and non-negative");
    }
    const double loss =
        checkedProduct(energy_loss_velocity_[species][face], owner_energy,
                       "wall energy loss overflowed");
    return loss - energy_inward_flux_[species][face];
  }

 private:
  struct ResolvedTarget {
    mesh::FaceId face;
    SpeciesId species;
  };

  struct ResolvedEmission {
    mesh::FaceId face;
    SpeciesId incident_species;
    std::size_t target_index;
    double yield;
    double emitted_mean_energy;
  };

  /** @brief Multiplies finite non-negative wall quantities safely. */
  [[nodiscard]] static double checkedProduct(double lhs, double rhs,
                                             const char* message) {
    if (lhs == 0.0 || rhs == 0.0) {
      return 0.0;
    }
    if (lhs > std::numeric_limits<double>::max() / rhs) {
      throw std::overflow_error(message);
    }
    return lhs * rhs;
  }

  /** @brief Adds a non-negative contribution while detecting overflow. */
  static void checkedAccumulate(double& destination, double contribution,
                                const char* message) {
    if (destination > std::numeric_limits<double>::max() - contribution) {
      throw std::overflow_error(message);
    }
    destination += contribution;
  }

  /** @brief Rejects invalid species IDs and wall coefficients. */
  void validateBoundary(const WallBoundary& boundary) const {
    if (boundary.losses.empty()) {
      throw std::invalid_argument("wall requires at least one species loss");
    }
    std::vector<SpeciesId> loss_species;
    loss_species.reserve(boundary.losses.size());
    for (const auto& loss : boundary.losses) {
      static_cast<void>(species_->at(loss.species));
      if (!std::isfinite(loss.particle_loss_velocity) ||
          loss.particle_loss_velocity < 0.0 ||
          !std::isfinite(loss.energy_loss_velocity) ||
          loss.energy_loss_velocity < 0.0) {
        throw std::invalid_argument(
            "wall loss velocities must be finite and non-negative");
      }
      for (const auto existing : loss_species) {
        if (existing == loss.species) {
          throw std::invalid_argument("duplicate wall loss species");
        }
      }
      loss_species.push_back(loss.species);
    }
    for (const auto& emission : boundary.secondary_emissions) {
      static_cast<void>(species_->at(emission.incident_species));
      static_cast<void>(species_->at(emission.emitted_species));
      if (!std::isfinite(emission.yield) || emission.yield < 0.0 ||
          !std::isfinite(emission.emitted_mean_energy) ||
          emission.emitted_mean_energy < 0.0) {
        throw std::invalid_argument(
            "secondary yield and emitted energy must be finite and "
            "non-negative");
      }
      const auto has_loss = [&loss_species](SpeciesId id) {
        for (const auto loss_id : loss_species) {
          if (loss_id == id) {
            return true;
          }
        }
        return false;
      };
      if (!has_loss(emission.incident_species) ||
          !has_loss(emission.emitted_species)) {
        throw std::invalid_argument(
            "secondary-emission species require wall loss definitions");
      }
    }
  }

  /** @brief Finds or creates the compact accumulator for one emitted flux. */
  [[nodiscard]] std::size_t resolveTarget(mesh::FaceId face,
                                          SpeciesId species) {
    for (std::size_t index = 0; index < resolved_targets_.size(); ++index) {
      const auto& target = resolved_targets_[index];
      if (target.face == face && target.species == species) {
        return index;
      }
    }
    resolved_targets_.push_back({face, species});
    return resolved_targets_.size() - 1;
  }

  /** @brief Resolves boundary IDs into face coefficients and emission terms. */
  void resolveBoundaries(const WallBoundarySet& boundaries) {
    std::unordered_map<mesh::BoundaryId, bool> found;
    for (const auto& [boundary_id, boundary] : boundaries) {
      if (boundary_id == mesh::invalid_boundary) {
        throw std::invalid_argument("invalid wall boundary ID");
      }
      validateBoundary(boundary);
      found.emplace(boundary_id, false);
    }

    for (mesh::FaceId face = 0; face < mesh_->numFaces(); ++face) {
      if (!mesh_->isBoundary(face)) {
        continue;
      }
      const auto boundary_id = mesh_->boundaryId(face);
      if (!boundaries.contains(boundary_id)) {
        continue;
      }
      found.at(boundary_id) = true;
      const auto& boundary = boundaries.at(boundary_id);
      for (const auto& loss : boundary.losses) {
        active_[loss.species][face] = std::uint8_t{1};
        particle_loss_velocity_[loss.species][face] =
            loss.particle_loss_velocity;
        energy_loss_velocity_[loss.species][face] = loss.energy_loss_velocity;
      }
      for (const auto& emission : boundary.secondary_emissions) {
        const auto target_index = resolveTarget(face, emission.emitted_species);
        resolved_emissions_.push_back({face, emission.incident_species,
                                       target_index, emission.yield,
                                       emission.emitted_mean_energy});
      }
    }

    for (const auto& [boundary_id, was_found] : found) {
      static_cast<void>(boundary_id);
      if (!was_found) {
        throw std::invalid_argument(
            "configured wall boundary ID does not exist in mesh");
      }
    }
  }

  /** @brief Validates the species-density container without scanning cells. */
  void validateDensityLayout(const SpeciesCellFields& density) const {
    if (&density.mesh() != mesh_ || density.size() != species_->size()) {
      throw std::invalid_argument("wall density layout differs");
    }
  }

  /** @brief Validates only owner densities consumed by emission channels. */
  void validateIncidentDensities(const SpeciesCellFields& density) const {
    for (const auto& emission : resolved_emissions_) {
      const double value =
          density[emission.incident_species][mesh_->owner(emission.face)];
      if (!std::isfinite(value) || value < 0.0) {
        throw std::invalid_argument(
            "wall density must be finite and non-negative");
      }
    }
  }

  /** @brief Prevents reading secondary fluxes before their first refresh. */
  void requireEvaluated() const {
    if (!evaluated_) {
      throw std::logic_error(
          "wall flux assembler must be evaluated before complete flux query");
    }
  }

  /** @brief Requires a valid species ID and configured boundary face. */
  void validateWallFace(SpeciesId species, mesh::FaceId face) const {
    static_cast<void>(species_->at(species));
    if (face >= mesh_->numFaces() || !mesh_->isBoundary(face) ||
        active_[species][face] == std::uint8_t{0}) {
      throw std::invalid_argument("species face has no wall-flux definition");
    }
  }

  const mesh::IMesh* mesh_;
  const SpeciesSet* species_;
  SpeciesWallFaceMasks active_;
  SpeciesFaceFields particle_loss_velocity_;
  SpeciesFaceFields energy_loss_velocity_;
  SpeciesFaceFields particle_inward_flux_;
  SpeciesFaceFields energy_inward_flux_;
  std::vector<ResolvedTarget> resolved_targets_;
  std::vector<ResolvedEmission> resolved_emissions_;
  std::vector<double> candidate_particle_inward_flux_;
  std::vector<double> candidate_energy_inward_flux_;
  bool evaluated_{false};
};

}  // namespace pemu::physics::wall
