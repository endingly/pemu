#pragma once

namespace pemu::physics {

enum class ChargePolarity { Positive, Negative };

[[nodiscard]]
constexpr double polaritySign(ChargePolarity polarity) noexcept {
  switch (polarity) {
    case ChargePolarity::Positive:
      return 1.0;
    case ChargePolarity::Negative:
      return -1.0;
  }
  return 0.0;
}

}  // namespace pemu::physics