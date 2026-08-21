#include <pemu/trace/renderer/diagnostic_event_renderer.hpp>

#include "event_table_format.hpp"

namespace pemu::trace::renderer {

void DiagnosticEventRenderer::render(std::ostream& stream,
                                     const TraceEvent& event,
                                     bool write_header) const {
  detail::renderEventTable(stream, event, write_header);
}

}  // namespace pemu::trace::renderer
