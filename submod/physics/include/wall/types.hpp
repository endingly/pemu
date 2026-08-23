#pragma once

#include <pemu/mesh/types.hpp>
#include <pemu/physics/species.hpp>

#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

namespace pemu::physics::wall {

/**
 * @brief Linear outward wall loss for one particle and energy state.
 *
 * The primary particle flux is v_loss*n and the primary energy flux is
 * v_energy*w. Reflection and sheath closures are represented by the effective
 * non-negative velocities supplied by the caller.
 */
struct SpeciesWallLoss {
  SpeciesId species;
  double particle_loss_velocity{};
  double energy_loss_velocity{};
};

/**
 * @brief Emits one target particle per incident wall loss with a fixed yield.
 *
 * emitted_mean_energy is energy carried by each emitted target particle in
 * the same energy unit used by the target energy-density equation.
 */
struct SecondaryEmission {
  SpeciesId incident_species;
  SpeciesId emitted_species;
  double yield{};
  double emitted_mean_energy{};
};

/** @brief Complete local interaction model assigned to one mesh boundary. */
struct WallBoundary {
  std::vector<SpeciesWallLoss> losses;
  std::vector<SecondaryEmission> secondary_emissions;
};

/** @brief Owns plasma-wall definitions indexed by mesh boundary ID. */
class WallBoundarySet {
 public:
  /** @brief Inserts or replaces the complete model for one boundary ID. */
  void set(mesh::BoundaryId id, WallBoundary boundary) {
    boundaries_.insert_or_assign(id, std::move(boundary));
  }

  /** @brief Reports whether an ID has a plasma-wall model. */
  [[nodiscard]] bool contains(mesh::BoundaryId id) const noexcept {
    return boundaries_.contains(id);
  }

  /** @brief Returns one configured wall or throws for an absent ID. */
  [[nodiscard]] const WallBoundary& at(mesh::BoundaryId id) const {
    return boundaries_.at(id);
  }

  /** @brief Returns the number of configured mesh boundary IDs. */
  [[nodiscard]] std::size_t size() const noexcept { return boundaries_.size(); }

  /** @brief Reports whether no plasma-wall boundary is configured. */
  [[nodiscard]] bool empty() const noexcept { return boundaries_.empty(); }

  /** @brief Removes all wall definitions. */
  void clear() noexcept { boundaries_.clear(); }

  /** @brief Provides read-only iteration for mesh/configuration validation. */
  [[nodiscard]] auto begin() const noexcept { return boundaries_.begin(); }

  /** @brief Returns the read-only end iterator. */
  [[nodiscard]] auto end() const noexcept { return boundaries_.end(); }

 private:
  std::unordered_map<mesh::BoundaryId, WallBoundary> boundaries_;
};

/** @brief Byte mask used to identify wall-overridden species face fluxes. */
using SpeciesWallFaceMasks = field::FaceFieldSet<std::uint8_t, SpeciesId>;

}  // namespace pemu::physics::wall
