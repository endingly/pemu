#include <gtest/gtest.h>

#include <pemu/trace/any_trace_sink.hpp>
#include <pemu/trace/diag.hpp>
#include <pemu/trace/ostream_trace_sink.hpp>
#include <pemu/trace/split_trace_sink.hpp>
#include <pemu/trace/statistics.hpp>
#include <pemu/trace/trace.hpp>

#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <streambuf>
#include <utility>

namespace pemu::trace::test {

namespace {

struct CountingSink {
  std::size_t* count{};
  std::size_t* flush_count{};

  void operator()(const TraceEvent&) const noexcept { ++*count; }
  void flush() const noexcept { ++*flush_count; }
};

struct EventOnlySink {
  void operator()(const TraceEvent&) const noexcept {}
};

struct ThrowingSink {
  void operator()(const TraceEvent&) const {}
  void flush() const noexcept {}
};

static_assert(TraceSink<CountingSink>);
static_assert(!TraceSink<EventOnlySink>);
static_assert(!TraceSink<ThrowingSink>);
static_assert(std::constructible_from<AnyTraceSink, CountingSink>);
static_assert(!std::constructible_from<AnyTraceSink, ThrowingSink>);

class FlushCountingBuffer final : public std::stringbuf {
 public:
  int flush_count{};

 protected:
  int sync() override {
    ++flush_count;
    return std::stringbuf::sync();
  }
};

[[nodiscard]] std::size_t countOccurrences(std::string_view text,
                                           std::string_view needle) {
  std::size_t count = 0;
  std::size_t position = 0;
  while ((position = text.find(needle, position)) != std::string_view::npos) {
    ++count;
    position += needle.size();
  }
  return count;
}

}  // namespace

TEST(TraceSinkTest, AnyTraceSinkForwardsNoexceptSink) {
  std::size_t count{};
  std::size_t flush_count{};
  AnyTraceSink sink{CountingSink{.count = &count, .flush_count = &flush_count}};
  sink({.domain = DiagDomain::simulation, .name = "test"});
  sink.flush();
  AnyTraceSink shared_sink = sink;
  shared_sink({.domain = DiagDomain::simulation, .name = "shared"});
  shared_sink.flush();
  EXPECT_EQ(count, 2u);
  EXPECT_EQ(flush_count, 2u);

  AnyTraceSink empty;
  EXPECT_NO_THROW(empty({.domain = DiagDomain::simulation, .name = "ignored"}));
  EXPECT_NO_THROW(empty.flush());
}

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
  EXPECT_NO_THROW(sink.flush());
}

TEST(PhysicalVolumeSemanticsTest, Extrudes2DMeasuresAndKeeps3DMeasures) {
  using namespace units;
  const PhysicalVolumeSemantics planar{
      .mesh_dimension = 2,
      .planar_depth = 2.5,
      .mesh_length_unit = precise::cm,
  };
  EXPECT_TRUE(planar.valid());
  EXPECT_DOUBLE_EQ(planar.cellVolume(3.0), 7.5);
  EXPECT_DOUBLE_EQ(planar.faceArea(4.0), 10.0);
  EXPECT_STREQ(planar.name(), "planar_extrusion");
  EXPECT_EQ(planar.cellWeightUnit(), precise::cm.pow(3));
  EXPECT_EQ(planar.faceWeightUnit(), precise::cm.pow(2));

  const PhysicalVolumeSemantics native_3d{
      .mesh_dimension = 3,
      .planar_depth = 99.0,
      .mesh_length_unit = precise::m,
  };
  EXPECT_TRUE(native_3d.valid());
  EXPECT_DOUBLE_EQ(native_3d.cellVolume(3.0), 3.0);
  EXPECT_DOUBLE_EQ(native_3d.faceArea(4.0), 4.0);
  EXPECT_STREQ(native_3d.name(), "native_3d");
}

TEST(StatisticsOptionsTest, SamplesOnlyAtConfiguredStepInterval) {
  using namespace mp_units;
  using namespace mp_units::si::unit_symbols;
  constexpr StatisticsOptions options{true, 3, 2.0 * cm, isq::length[cm]};
  static_assert(options.valid());
  EXPECT_TRUE(options.shouldSample(0));
  EXPECT_FALSE(options.shouldSample(1));
  EXPECT_FALSE(options.shouldSample(2));
  EXPECT_TRUE(options.shouldSample(3));

  constexpr StatisticsOptions disabled{};
  EXPECT_FALSE(disabled.shouldSample(0));
  constexpr StatisticsOptions invalid{true, 0, 1.0 * cm, isq::length[cm]};
  EXPECT_FALSE(invalid.valid());
  EXPECT_FALSE(invalid.shouldSample(0));
}

TEST(ScalarFieldStatisticsTest, ComputesWeightedExtremaIntegralMeanAndRms) {
  ScalarFieldStatisticsAccumulator accumulator;
  accumulator.add(1.0, 2.0);
  accumulator.add(3.0, 6.0);

  const auto statistics = accumulator.finish();
  EXPECT_TRUE(statistics.valid());
  EXPECT_EQ(statistics.sample_count, 2u);
  EXPECT_EQ(statistics.finite_count, 2u);
  EXPECT_EQ(statistics.non_finite_count, 0u);
  EXPECT_EQ(statistics.negative_count, 0u);
  EXPECT_DOUBLE_EQ(statistics.minimum, 1.0);
  EXPECT_DOUBLE_EQ(statistics.maximum, 3.0);
  EXPECT_DOUBLE_EQ(statistics.max_abs, 3.0);
  EXPECT_DOUBLE_EQ(statistics.weight_sum, 8.0);
  EXPECT_DOUBLE_EQ(statistics.integral, 20.0);
  EXPECT_DOUBLE_EQ(statistics.l1_integral, 20.0);
  EXPECT_DOUBLE_EQ(statistics.weighted_mean, 2.5);
  EXPECT_DOUBLE_EQ(statistics.weighted_rms, std::sqrt(7.0));
  EXPECT_EQ(statistics.value_unit, units::precise::one);
  EXPECT_EQ(statistics.weight_unit, units::precise::one);
  EXPECT_EQ(statistics.integralUnit(), units::precise::one);
}

TEST(ScalarFieldStatisticsTest, DerivesIntegralUnitOutsideSampleLoop) {
  ScalarFieldStatisticsAccumulator accumulator;
  accumulator.add(2.0, 4.0);

  const auto statistics = accumulator.finish(
      units::precise::C / units::precise::cm.pow(3), units::precise::cm.pow(3));
  EXPECT_EQ(statistics.value_unit,
            units::precise::C / units::precise::cm.pow(3));
  EXPECT_EQ(statistics.weight_unit, units::precise::cm.pow(3));
  EXPECT_EQ(statistics.integralUnit(), units::precise::C);
  EXPECT_DOUBLE_EQ(statistics.integral, 8.0);
}

TEST(ScalarFieldStatisticsTest, ReportsNonFiniteValuesAndInvalidWeights) {
  ScalarFieldStatisticsAccumulator accumulator;
  accumulator.add(-2.0, 1.0);
  accumulator.add(std::numeric_limits<double>::infinity(), 1.0);
  accumulator.add(4.0, 0.0);

  const auto statistics = accumulator.finish();
  EXPECT_FALSE(statistics.valid());
  EXPECT_EQ(statistics.sample_count, 3u);
  EXPECT_EQ(statistics.finite_count, 2u);
  EXPECT_EQ(statistics.non_finite_count, 1u);
  EXPECT_EQ(statistics.negative_count, 1u);
  EXPECT_EQ(statistics.invalid_weight_count, 1u);
  EXPECT_DOUBLE_EQ(statistics.minimum, -2.0);
  EXPECT_DOUBLE_EQ(statistics.maximum, 4.0);
}

TEST(TraceSinkTest, OstreamSinkFormatsTabularStructuredEvents) {
  const std::array attributes{
      TraceAttribute{"step", std::uint64_t{3}},
      TraceAttribute{"time", 0.25},
      TraceAttribute{"finished", true},
      TraceAttribute{"offset", std::int64_t{-2}},
      TraceAttribute{"state", std::string_view{"running"}},
      TraceAttribute{"solver_status",
                     std::string_view{"PatternAnalysisFailed"}},
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
  EXPECT_NE(
      output.str().find("Info     | simulation.fixed_step.step.completed"),
      std::string::npos);
  EXPECT_NE(output.str().find("|    3 |            0.25"), std::string::npos);
  EXPECT_NE(output.str().find("finished=true; offset=-2; state=running"),
            std::string::npos);
  EXPECT_NE(output.str().find("|    PatternAnalysisFailed |"),
            std::string::npos);
  EXPECT_NE(output.str().find("END TIME"), std::string::npos);
  EXPECT_NE(output.str().find("REL RESIDUAL"), std::string::npos);
  EXPECT_FALSE(sink.failed());
}

TEST(TraceSinkTest, StatisticsRendererGroupsRowsAndFormatsDedicatedColumns) {
  const std::array species_attributes{
      TraceAttribute{"step", std::uint64_t{3}},
      TraceAttribute{"time", 0.25, units::precise::s},
      TraceAttribute{"species", std::string_view{"e"}},
      TraceAttribute{"samples", std::uint64_t{4096}},
      TraceAttribute{"non_finite", std::uint64_t{0}},
      TraceAttribute{"negative", std::uint64_t{0}},
      TraceAttribute{"minimum", 1.0,
                     units::precise::one / units::precise::cm.pow(3)},
      TraceAttribute{"maximum", 3.0,
                     units::precise::one / units::precise::cm.pow(3)},
      TraceAttribute{"volume_mean", 2.0,
                     units::precise::one / units::precise::cm.pow(3)},
      TraceAttribute{"rms", 2.5,
                     units::precise::one / units::precise::cm.pow(3)},
      TraceAttribute{"total_number", 8.0, units::precise::one},
      TraceAttribute{"physical_volume", 1.0, units::precise::cm.pow(3)},
      TraceAttribute{"volume_semantics", std::string_view{"planar_extrusion"}},
  };
  const std::array field_attributes{
      TraceAttribute{"step", std::uint64_t{3}},
      TraceAttribute{"time", 0.25, units::precise::s},
      TraceAttribute{"samples", std::uint64_t{4096}},
      TraceAttribute{"non_finite", std::uint64_t{0}},
      TraceAttribute{"minimum", 1.25, units::precise::V},
      TraceAttribute{"maximum", 400.0, units::precise::V},
      TraceAttribute{"volume_mean", 200.0, units::precise::V},
      TraceAttribute{"rms", 230.0, units::precise::V},
      TraceAttribute{"volume_integral", 5.0,
                     units::precise::V * units::precise::cm.pow(3)},
  };
  const std::array next_step_attributes{
      TraceAttribute{"step", std::uint64_t{4}},
      TraceAttribute{"minimum", 2.0, units::precise::V},
      TraceAttribute{"maximum", 401.0, units::precise::V},
      TraceAttribute{"volume_mean", 201.0, units::precise::V},
      TraceAttribute{"rms", 231.0, units::precise::V},
      TraceAttribute{"volume_integral", 6.0,
                     units::precise::V * units::precise::cm.pow(3)},
  };
  const TraceEvent species_event{
      .output_channel = OutputChannel::Statistics,
      .domain = DiagDomain::physics,
      .category = "species",
      .name = "statistics",
      .attributes = species_attributes,
  };
  const TraceEvent field_event{
      .output_channel = OutputChannel::Statistics,
      .domain = DiagDomain::field,
      .category = "potential",
      .name = "statistics",
      .attributes = field_attributes,
  };
  const TraceEvent next_step_event{
      .output_channel = OutputChannel::Statistics,
      .domain = DiagDomain::field,
      .category = "potential",
      .name = "statistics",
      .attributes = next_step_attributes,
  };
  std::ostringstream output;
  OstreamTraceSink sink(output);

  sink(species_event);
  sink(field_event);
  sink(next_step_event);
  sink.flush();

  EXPECT_NE(output.str().find("STATISTICS | STEP=3 | TIME=0.25 [s]"),
            std::string::npos);
  EXPECT_NE(output.str().find("STATISTICS | STEP=4"), std::string::npos);
  EXPECT_EQ(countOccurrences(output.str(), "STATISTICS | STEP="), 2u);
  EXPECT_NE(output.str().find("FIELD"), std::string::npos);
  EXPECT_NE(output.str().find("UNIT"), std::string::npos);
  EXPECT_NE(output.str().find("MIN"), std::string::npos);
  EXPECT_NE(output.str().find("MAX"), std::string::npos);
  EXPECT_NE(output.str().find("MEAN"), std::string::npos);
  EXPECT_NE(output.str().find("RMS"), std::string::npos);
  EXPECT_NE(output.str().find("INTEGRAL"), std::string::npos);
  EXPECT_NE(output.str().find("species[e]"), std::string::npos);
  EXPECT_NE(output.str().find("potential"), std::string::npos);
  EXPECT_NE(output.str().find("1/mL"), std::string::npos);
  EXPECT_NE(output.str().find("8 [1]"), std::string::npos);
  EXPECT_EQ(countOccurrences(output.str(), "CELL_SAMPLES=4096"), 1u);
  EXPECT_EQ(countOccurrences(output.str(), "VOLUME_SEMANTICS="), 1u);
  EXPECT_EQ(output.str().find("DETAILS"), std::string::npos);
  EXPECT_EQ(output.str().find("non_finite=0"), std::string::npos);
  EXPECT_EQ(output.str().find("negative=0"), std::string::npos);
}

TEST(TraceSinkTest, StatisticsRendererOnlyShowsNonZeroCountsAsDiagnostics) {
  const std::array attributes{
      TraceAttribute{"step", std::uint64_t{4}},
      TraceAttribute{"species", std::string_view{"e"}},
      TraceAttribute{"non_finite", std::uint64_t{0}},
      TraceAttribute{"negative", std::uint64_t{2}},
      TraceAttribute{"minimum", -1.0,
                     units::precise::one / units::precise::cm.pow(3)},
      TraceAttribute{"maximum", 3.0,
                     units::precise::one / units::precise::cm.pow(3)},
      TraceAttribute{"volume_mean", 1.0,
                     units::precise::one / units::precise::cm.pow(3)},
      TraceAttribute{"rms", 2.0,
                     units::precise::one / units::precise::cm.pow(3)},
      TraceAttribute{"total_number", 4.0, units::precise::one},
  };
  const TraceEvent event{
      .kind = EventKind::Diagnostic,
      .output_channel = OutputChannel::Statistics,
      .domain = DiagDomain::physics,
      .category = "species",
      .name = "statistics",
      .severity = Severity::Warning,
      .attributes = attributes,
  };
  std::ostringstream output;
  OstreamTraceSink sink(output);

  sink(event);
  sink.flush();

  EXPECT_NE(output.str().find("STATISTICS DIAGNOSTICS"), std::string::npos);
  EXPECT_NE(output.str().find("Warning"), std::string::npos);
  EXPECT_NE(output.str().find("Diagnostic"), std::string::npos);
  EXPECT_NE(output.str().find("negative=2"), std::string::npos);
  EXPECT_EQ(output.str().find("non_finite=0"), std::string::npos);
}

TEST(TraceSinkTest, OstreamSinkFlushesOnlyWhenContractRequestsIt) {
  FlushCountingBuffer buffer;
  std::ostream output(&buffer);
  OstreamTraceSink sink(output);
  const TraceEvent step_completed{
      .domain = DiagDomain::simulation,
      .category = "adaptive_step",
      .name = "step.completed",
  };
  sink(step_completed);
  EXPECT_EQ(buffer.flush_count, 0);
  sink(step_completed);
  EXPECT_EQ(buffer.flush_count, 0);
  sink.flush();
  EXPECT_EQ(buffer.flush_count, 1);
  EXPECT_FALSE(sink.failed());
}

TEST(TraceSinkTest, SplitSinkRoutesChannelsAndForwardsFlush) {
  std::ostringstream diagnostic_output;
  FlushCountingBuffer statistics_buffer;
  std::ostream statistics_output(&statistics_buffer);
  SplitTraceSink sink{OstreamTraceSink{diagnostic_output},
                      OstreamTraceSink{statistics_output}};

  const TraceEvent abnormal_statistics{
      .kind = EventKind::Diagnostic,
      .output_channel = OutputChannel::Statistics,
      .domain = DiagDomain::physics,
      .category = "species",
      .name = "statistics",
      .severity = Severity::Warning,
  };
  const TraceEvent step_completed{
      .domain = DiagDomain::simulation,
      .category = "adaptive_step",
      .name = "step.completed",
  };

  sink(abnormal_statistics);
  sink(step_completed);
  EXPECT_EQ(statistics_buffer.flush_count, 0);
  sink(step_completed);
  EXPECT_EQ(statistics_buffer.flush_count, 0);
  sink.flush();

  EXPECT_NE(statistics_buffer.str().find("STATISTICS"), std::string::npos);
  EXPECT_NE(statistics_buffer.str().find("species"), std::string::npos);
  EXPECT_EQ(statistics_buffer.str().find("simulation.adaptive_step"),
            std::string::npos);
  EXPECT_NE(diagnostic_output.str().find("simulation.adaptive_step"),
            std::string::npos);
  EXPECT_EQ(diagnostic_output.str().find("physics.species.statistics"),
            std::string::npos);
  EXPECT_EQ(statistics_buffer.flush_count, 1);
}

TEST(EnumStringTest, DomainsAndSeveritiesHaveStableNames) {
  constexpr std::array domains{
      std::pair{DiagDomain::linalg, "linalg"},
      std::pair{DiagDomain::mesh, "mesh"},
      std::pair{DiagDomain::unit, "unit"},
      std::pair{DiagDomain::field, "field"},
      std::pair{DiagDomain::output, "output"},
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

  EXPECT_EQ(pemu::to_string(EventKind::Trace), "Trace");
  EXPECT_EQ(pemu::to_string(EventKind::Diagnostic), "Diagnostic");
  EXPECT_EQ(pemu::to_string(OutputChannel::Diagnostic), "Diagnostic");
  EXPECT_EQ(pemu::to_string(OutputChannel::Statistics), "Statistics");
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
  const TraceEvent ordinary_before{
      .domain = DiagDomain::simulation,
      .category = "adaptive_step",
      .name = "step.started",
  };
  const TraceEvent ordinary_after{
      .domain = DiagDomain::simulation,
      .category = "adaptive_step",
      .name = "step.failed",
      .severity = Severity::Error,
  };

  sink(ordinary_before);
  sink(diagnostic);
  sink(ordinary_after);

  EXPECT_EQ(countOccurrences(output.str(), "LEVEL    | EVENT"), 1u);
  EXPECT_EQ(output.str().find("\nTRACE\n"), std::string::npos);
  EXPECT_EQ(output.str().find("\nDIAGNOSTIC\n"), std::string::npos);
  EXPECT_NE(output.str().find("Warning  | equation.poisson.not_converged"),
            std::string::npos);
  EXPECT_NE(
      output.str().find("linear solve did not converge; residual=1.25e-06; "
                        "iterations=7"),
      std::string::npos);
  EXPECT_NE(output.str().find("simulation.adaptive_step.step.failed"),
            std::string::npos);
  EXPECT_FALSE(sink.failed());
}

}  // namespace pemu::trace::test
