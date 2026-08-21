#pragma once

#include <pemu/trace/trace.hpp>

#include <iosfwd>

namespace pemu::trace::renderer {

class OrdinaryEventRenderer {
 public:
  void render(std::ostream& stream, const TraceEvent& event,
              bool write_header) const;
};

}  // namespace pemu::trace::renderer
