#pragma once

#include <cstdint>
#include <string_view>

namespace pemu::trace {

// Diagnostic domains deliberately mirror pemu module names. Keep the spelling
// stable: sinks may use it for filtering and machine-readable output.
enum class DiagDomain : std::uint8_t {
  linalg,
  mesh,
  unit,
  field,
  trace,
  boundary,
  discretization,
  physics,
  equation,
  simulation,
};

[[nodiscard]] std::string_view diagDomainName(DiagDomain domain) noexcept;

}  // namespace pemu::trace
