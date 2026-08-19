#pragma once

#include <pemu/physics/charge_polarity.hpp>

#include <stdexcept>

namespace pemu::physics {

struct ChargedSpeciesTransport {
  double charge{};
  double mobility{};
  double diffusivity{};

  void validate() const {
    if (charge == 0.0) {
      throw std::invalid_argument(
          "charged species must have "
          "non-zero charge");
    }
    if (mobility < 0.0) {
      throw std::invalid_argument("mobility must be non-negative");
    }
    if (diffusivity <= 0.0) {
      throw std::invalid_argument("diffusivity must be positive");
    }
  }

  [[nodiscard]]
  ChargePolarity polarity() const noexcept {
    return charge > 0.0 ? ChargePolarity::Positive : ChargePolarity::Negative;
  }
};

}  // namespace pemu::physics