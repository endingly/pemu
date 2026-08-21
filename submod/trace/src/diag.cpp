#include <magic_enum/magic_enum.hpp>
#include <pemu/trace/diag.hpp>

namespace pemu::trace {

std::string_view diagDomainName(DiagDomain domain) noexcept {
  return magic_enum::enum_name(domain);
}

}  // namespace pemu::trace
