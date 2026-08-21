#include <gtest/gtest.h>

#include <pemu/trace/diag.hpp>
#include <pemu/trace/ostream_trace_sink.hpp>
#include <pemu/trace/trace.hpp>

#include <array>
#include <sstream>
#include <streambuf>
#include <utility>

namespace pemu::trace::test {

namespace {

class FlushCountingBuffer final : public std::stringbuf {
 public:
  int flush_count{};

 protected:
  int sync() override {
    ++flush_count;
    return std::stringbuf::sync();
  }
};

}  // namespace

TEST(TraceSinkTest, NullSinkAcceptsStructuredEvents) {
  const std::array attributes{
      TraceAttribute{"step", std::uint64_t{3}},
      TraceAttribute{"time", 0.25},
  };
  const TraceEvent event{
      .domain = DiagDomain::simulation,
      .category = "fixed_step",
      .name = "step.completed",
      .severity = Severity::Trace,
      .attributes = attributes,
  };

  NullTraceSink sink;
  EXPECT_NO_THROW(sink(event));
}

TEST(TraceSinkTest, OstreamSinkFormatsTabularStructuredEvents) {
  const std::array attributes{
      TraceAttribute{"step", std::uint64_t{3}},
      TraceAttribute{"time", 0.25},
      TraceAttribute{"finished", true},
      TraceAttribute{"offset", std::int64_t{-2}},
      TraceAttribute{"state", std::string_view{"running"}},
      TraceAttribute{"solver_status", std::string_view{"PatternAnalysisFailed"}},
  };
  const TraceEvent event{
      .domain = DiagDomain::simulation,
      .category = "fixed_step",
      .name = "step.completed",
      .severity = Severity::Info,
      .attributes = attributes,
  };
  std::ostringstream output;
  OstreamTraceSink sink(output);

  sink(event);

  EXPECT_NE(output.str().find("LEVEL    | EVENT"), std::string::npos);
  EXPECT_NE(output.str().find("Info     | simulation.fixed_step.step.completed"),
            std::string::npos);
  EXPECT_NE(output.str().find("|    3 |            0.25"), std::string::npos);
  EXPECT_NE(output.str().find("finished=true; offset=-2; state=running"),
            std::string::npos);
  EXPECT_NE(output.str().find("|    PatternAnalysisFailed |"),
            std::string::npos);
  EXPECT_FALSE(sink.failed());
}

TEST(TraceSinkTest, OstreamSinkFlushesEveryTwoSimulationStepsAndOnRunEnd) {
  FlushCountingBuffer buffer;
  std::ostream output(&buffer);
  OstreamTraceSink sink(output);
  const TraceEvent step_completed{
      .domain = DiagDomain::simulation,
      .category = "adaptive_step",
      .name = "step.completed",
  };
  const TraceEvent run_completed{
      .domain = DiagDomain::simulation,
      .category = "adaptive_step",
      .name = "run.completed",
  };

  sink(step_completed);
  EXPECT_EQ(buffer.flush_count, 0);
  sink(step_completed);
  EXPECT_EQ(buffer.flush_count, 1);
  sink(step_completed);
  EXPECT_EQ(buffer.flush_count, 1);
  sink(run_completed);
  EXPECT_EQ(buffer.flush_count, 2);
  EXPECT_FALSE(sink.failed());
}

TEST(EnumStringTest, DomainsAndSeveritiesHaveStableNames) {
  constexpr std::array domains{
      std::pair{DiagDomain::linalg, "linalg"},
      std::pair{DiagDomain::mesh, "mesh"},
      std::pair{DiagDomain::unit, "unit"},
      std::pair{DiagDomain::field, "field"},
      std::pair{DiagDomain::trace, "trace"},
      std::pair{DiagDomain::boundary, "boundary"},
      std::pair{DiagDomain::discretization, "discretization"},
      std::pair{DiagDomain::physics, "physics"},
      std::pair{DiagDomain::equation, "equation"},
      std::pair{DiagDomain::simulation, "simulation"},
  };

  for (const auto& [domain, expected_name] : domains) {
    EXPECT_EQ(pemu::to_string(domain), expected_name);
    EXPECT_EQ(diagDomainName(domain), expected_name);
  }

  constexpr std::array severities{
      std::pair{Severity::Trace, "Trace"},
      std::pair{Severity::Debug, "Debug"},
      std::pair{Severity::Info, "Info"},
      std::pair{Severity::Warning, "Warning"},
      std::pair{Severity::Error, "Error"},
      std::pair{Severity::Critical, "Critical"},
  };
  for (const auto& [severity, expected_name] : severities) {
    EXPECT_EQ(pemu::to_string(severity), expected_name);
  }
}

TEST(TraceSinkTest, NullSinkAcceptsStructuredDiagnostics) {
  const TraceEvent diagnostic{
      .kind = EventKind::Diagnostic,
      .domain = DiagDomain::equation,
      .category = "poisson",
      .name = "incompatible_rhs",
      .message = "pure-Neumann source violates compatibility",
      .severity = Severity::Error,
  };

  NullTraceSink sink;
  EXPECT_NO_THROW(sink(diagnostic));
}

TEST(TraceSinkTest, OstreamSinkFormatsDiagnosticEvent) {
  const std::array attributes{
      TraceAttribute{"residual", 1.25e-6},
      TraceAttribute{"iterations", std::uint64_t{7}},
  };
  const TraceEvent diagnostic{
      .kind = EventKind::Diagnostic,
      .domain = DiagDomain::equation,
      .category = "poisson",
      .name = "not_converged",
      .message = "linear solve did not converge",
      .severity = Severity::Warning,
      .attributes = attributes,
  };
  std::ostringstream output;
  OstreamTraceSink sink(output);

  sink(diagnostic);

  EXPECT_NE(output.str().find("Warning  | equation.poisson.not_converged"),
            std::string::npos);
  EXPECT_NE(output.str().find("linear solve did not converge; residual=1.25e-06; "
                               "iterations=7"),
            std::string::npos);
  EXPECT_FALSE(sink.failed());
}

}  // namespace pemu::trace::test
