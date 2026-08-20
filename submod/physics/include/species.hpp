#pragma once

#include <cstdint>
#include <limits>
#include <pemu/field/field_set.hpp>
#include <pemu/physics/charge_polarity.hpp>
#include <stdexcept>
#include <string>
#include <vector>

namespace pemu::physics {

enum class SpeciesTransportModel { Immobile, DriftDiffusion };

struct SpeciesId {
  std::uint32_t value{};
  friend constexpr bool operator==(SpeciesId, SpeciesId) noexcept = default;
};

inline constexpr SpeciesId invalid_species{
    std::numeric_limits<std::uint32_t>::max()};

using SpeciesCellFields = field::CellFieldSet<double, SpeciesId>;
using SpeciesFaceFields = field::FaceFieldSet<double, SpeciesId>;

struct SpeciesProperties {
  std::string name;

  double charge{};
  double mobility{};
  double diffusivity{};

  SpeciesTransportModel transport_model{SpeciesTransportModel::Immobile};

  [[nodiscard]]
  bool isCharged() const noexcept {
    return charge != 0.0;
  }

  [[nodiscard]]
  bool isTransported() const noexcept {
    return transport_model == SpeciesTransportModel::DriftDiffusion;
  }

  [[nodiscard]]
  ChargePolarity polarity() const {
    if (!isCharged()) {
      throw std::logic_error("neutral species has no charge polarity");
    }

    return charge > 0.0 ? ChargePolarity::Positive : ChargePolarity::Negative;
  }

  void validate() const {
    if (name.empty()) {
      throw std::invalid_argument("species name must not be empty");
    }

    if (mobility < 0.0) {
      throw std::invalid_argument("mobility must be non-negative");
    }

    if (diffusivity < 0.0) {
      throw std::invalid_argument("diffusivity must be non-negative");
    }

    if (isTransported() && diffusivity <= 0.0) {

      throw std::invalid_argument(
          "drift-diffusion species "
          "requires positive diffusivity");
    }

    if (isTransported() && !isCharged() && mobility != 0.0) {

      throw std::invalid_argument(
          "neutral drift-diffusion species "
          "cannot have electric mobility");
    }
  }
};

class SpeciesSet {
 public:
  [[nodiscard]]
  SpeciesId add(SpeciesProperties properties) {
    properties.validate();

    for (const auto& existing : species_) {

      if (existing.name == properties.name) {

        throw std::invalid_argument("duplicate species name");
      }
    }

    const SpeciesId id{static_cast<std::uint32_t>(species_.size())};

    species_.push_back(std::move(properties));

    return id;
  }

  [[nodiscard]]
  std::size_t size() const noexcept {
    return species_.size();
  }

  [[nodiscard]]
  const SpeciesProperties& at(SpeciesId id) const {
    const auto index = static_cast<std::size_t>(id.value);

    if (index >= species_.size()) {
      throw std::out_of_range("invalid SpeciesId");
    }

    return species_[index];
  }

  [[nodiscard]]
  SpeciesId find(std::string_view name) const {
    for (std::size_t i = 0; i < species_.size(); ++i) {

      if (species_[i].name == name) {
        return SpeciesId{static_cast<std::uint32_t>(i)};
      }
    }

    return invalid_species;
  }

 private:
  std::vector<SpeciesProperties> species_;
};

}  // namespace pemu::physics
