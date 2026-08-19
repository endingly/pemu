#pragma once

#include <cmath>

namespace pemu::discretization::operators {

[[nodiscard]]
inline double bernoulli(const double x) noexcept {
  //
  // Taylor expansion around x = 0:
  //
  // B(x)
  //
  //   = 1
  //     - x/2
  //     + x^2/12
  //     - x^4/720
  //     + O(x^6)
  //
  if (std::abs(x) < 1e-6) {

    const double x2 = x * x;

    return 1.0 - 0.5 * x + x2 / 12.0 - x2 * x2 / 720.0;
  }

  //
  // For large positive x:
  //
  //     B(x)
  //
  //       = x exp(-x)
  //         ------------
  //         1 - exp(-x)
  //
  // avoids exp(+x) overflow.
  //
  if (x > 50.0) {

    const double e = std::exp(-x);

    return x * e / (1.0 - e);
  }

  //
  // expm1(x) computes:
  //
  //     exp(x) - 1
  //
  // accurately for moderately small x.
  //
  return x / std::expm1(x);
}

}  // namespace pemu::discretization::operators