#pragma once

#include <pemu/output/output_types.hpp>
#include <pemu/trace/trace.hpp>

#include <array>
#include <string>
#include <string_view>

namespace pemu::output {

template <trace::TraceSink Sink>
void traceOutputCompleted(Sink& sink, const OutputRecord& record) {
  const std::string path = record.path.string();
  const std::array attributes{
      trace::TraceAttribute{"path", std::string_view{path}},
      trace::TraceAttribute{"step", record.stamp.step},
      trace::TraceAttribute{"time", record.stamp.time},
  };
  const trace::TraceEvent event{
      .domain = trace::DiagDomain::output,
      .name = "completed",
      .severity = trace::Severity::Info,
      .attributes = attributes,
  };

  sink(event);
}

}  // namespace pemu::output
