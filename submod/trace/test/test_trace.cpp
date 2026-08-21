#include <gtest/gtest.h>

#include <pemu/trace/diag.hpp>
#include <pemu/trace/ostream_trace_sink.hpp>
#include <pemu/trace/trace.hpp>

#include <array>
#include <sstream>
#include <utility>

namespace pemu::trace::test {

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

TEST(TraceSinkTest, OstreamSinkFormatsOneStructuredEventPerLine) {
  const std::array attributes{
      TraceAttribute{"step", std::uint64_t{3}},
      TraceAttribute{"time", 0.25},
      TraceAttribute{"finished", true},
      TraceAttribute{"offset", std::int64_t{-2}},
      TraceAttribute{"state", std::string_view{"running"}},
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

  EXPECT_EQ(output.str(),
            "[Info] simulation.fixed_step.step.completed step=3 time=0.25 "
            "finished=1 offset=-2 state=running\n");
  EXPECT_FALSE(sink.failed());
}

TEST(DiagDomainTest, NamesMatchModuleNames) {
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
    EXPECT_EQ(diagDomainName(domain), expected_name);
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

  EXPECT_EQ(output.str(),
            "[Warning] equation.poisson.not_converged: linear solve did not "
            "converge residual=1.25e-06 iterations=7\n");
  EXPECT_FALSE(sink.failed());
}

}  // namespace pemu::trace::test
