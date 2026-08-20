#include <gtest/gtest.h>

#include <pemu/trace/ostream_trace_sink.hpp>
#include <pemu/trace/trace.hpp>

#include <array>
#include <sstream>

namespace pemu::trace::test {

TEST(TraceSinkTest, NullSinkAcceptsStructuredEvents) {
  const std::array attributes{
      TraceAttribute{"step", std::uint64_t{3}},
      TraceAttribute{"time", 0.25},
  };
  const TraceEvent event{
      .category = "simulation",
      .name = "step.completed",
      .severity = Severity::Trace,
      .attributes = attributes,
  };

  NullTraceSink sink;
  EXPECT_NO_THROW(sink(event));
}

TEST(TraceSinkTest, OstreamSinkFormatsOneStructuredEventPerLine) {
  const std::array attributes{
      TraceAttribute{"step", std::uint64_t{3}},
      TraceAttribute{"time", 0.25},
      TraceAttribute{"finished", true},
  };
  const TraceEvent event{
      .category = "simulation",
      .name = "step.completed",
      .severity = Severity::Info,
      .attributes = attributes,
  };
  std::ostringstream output;
  OstreamTraceSink sink(output);

  sink(event);

  EXPECT_EQ(output.str(),
            "[info] simulation.step.completed step=3 time=0.25 finished=1\n");
  EXPECT_FALSE(sink.failed());
}

}  // namespace pemu::trace::test
