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
  output,
  trace,
  boundary,
  discretization,
  physics,
  equation,
  simulation,
};

}  // namespace pemu::trace

namespace pemu {

[[nodiscard]] constexpr std::string_view to_string(
    trace::DiagDomain domain) noexcept {
  switch (domain) {
    case trace::DiagDomain::linalg:
      return "linalg";
    case trace::DiagDomain::mesh:
      return "mesh";
    case trace::DiagDomain::unit:
      return "unit";
    case trace::DiagDomain::field:
      return "field";
    case trace::DiagDomain::output:
      return "output";
    case trace::DiagDomain::trace:
      return "trace";
    case trace::DiagDomain::boundary:
      return "boundary";
    case trace::DiagDomain::discretization:
      return "discretization";
    case trace::DiagDomain::physics:
      return "physics";
    case trace::DiagDomain::equation:
      return "equation";
    case trace::DiagDomain::simulation:
      return "simulation";
  }
  return "Unknown";
}

}  // namespace pemu

namespace pemu::trace {

// Compatibility spelling for existing callers. New enum-to-text conversions
// should use pemu::to_string(enum_value).
[[nodiscard]] constexpr std::string_view diagDomainName(
    DiagDomain domain) noexcept {
  return pemu::to_string(domain);
}

}  // namespace pemu::trace
