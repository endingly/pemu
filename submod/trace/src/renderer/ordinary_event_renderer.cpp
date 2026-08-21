#include <pemu/trace/renderer/ordinary_event_renderer.hpp>

#include "event_table_format.hpp"

namespace pemu::trace::renderer {

void OrdinaryEventRenderer::render(std::ostream& stream,
                                   const TraceEvent& event,
                                   bool write_header) const {
  detail::renderEventTable(stream, event, write_header);
}

}  // namespace pemu::trace::renderer
