#pragma once

#include <pemu/linalg/i_solver.hpp>
#include <pemu/output/dump/trace.hpp>
#include <pemu/physics/electron_energy.hpp>
#include <pemu/simulation/checkpoint.hpp>
#include <pemu/simulation/detail/plasma_checkpoint.hpp>
#include <pemu/simulation/detail/plasma_field_output.hpp>
#include <pemu/simulation/detail/plasma_statistics.hpp>
#include <pemu/simulation/electron_energy.hpp>
#include <pemu/simulation/plasma_reaction_rate_evaluator.hpp>
#include <pemu/simulation/simulation_state.hpp>
#include <pemu/trace/any_trace_sink.hpp>

#include <llnl-units/units.hpp>

#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace pemu::simulation::detail {

inline constexpr std::size_t trace_flush_step_interval = 2;

template <std::size_t N>
void emitTrace(
    const trace::AnyTraceSink& sink, std::string_view category,
    std::string_view name, trace::Severity severity,
    const std::array<trace::TraceAttribute, N>& attributes) noexcept {
  sink({.domain = trace::DiagDomain::simulation,
        .category = category,
        .name = name,
        .severity = severity,
        .attributes = attributes});
}

inline void emitDiagnostic(const trace::AnyTraceSink& sink,
                           const linalg::SolverResult& result, std::size_t step,
                           double time) noexcept {
  if (!result.diagnostic.has_value()) {
    return;
  }
  const std::array attributes{
      trace::TraceAttribute{"step", static_cast<std::uint64_t>(step)},
      trace::TraceAttribute{"solver_status", to_string(result.status)},
      trace::TraceAttribute{"time", time},
      trace::TraceAttribute{"residual_norm", result.residual_norm},
      trace::TraceAttribute{"relative_residual", result.relative_residual},
  };
  auto diagnostic = *result.diagnostic;
  diagnostic.attributes = attributes;
  sink(diagnostic);
}

inline void emitSolverResult(const trace::AnyTraceSink& sink,
                             std::string_view category, std::string_view name,
                             trace::Severity severity,
                             const linalg::SolverResult& result,
                             std::size_t step) noexcept {
  const std::array attributes{
      trace::TraceAttribute{"step", static_cast<std::uint64_t>(step)},
      trace::TraceAttribute{"solver_status", to_string(result.status)},
      trace::TraceAttribute{"residual_norm", result.residual_norm},
      trace::TraceAttribute{"relative_residual", result.relative_residual},
  };
  emitTrace(sink, category, name, severity, attributes);
}

inline std::array<output::checkpoint::ScalarTarget, 2>
plasmaWorkflowScalarTargets(double& schema_version, double& workflow_kind) {
  return {
      output::checkpoint::ScalarTarget{
          .value = &schema_version,
          .key = std::string{plasma_checkpoint_schema_key},
          .metadata = {.name = std::string{plasma_checkpoint_schema_key}}},
      output::checkpoint::ScalarTarget{
          .value = &workflow_kind,
          .key = std::string{plasma_checkpoint_kind_key},
          .metadata = {.name = std::string{plasma_checkpoint_kind_key}}},
  };
}

inline void validatePlasmaWorkflowSchema(double schema_version,
                                         double workflow_kind,
                                         double expected_kind) {
  if (schema_version != plasma_checkpoint_schema_version ||
      workflow_kind != expected_kind) {
    throw std::invalid_argument("checkpoint workflow schema differs");
  }
}

/** @brief Owns the mandatory electron-energy operator and workflow workspaces. */
template <typename Stepper>
struct ElectronEnergySubsystem {
  field::CellField<double>* energy_density;
  physics::SpeciesId electron;
  double density_floor;
  equation::ExplicitElectronEnergyStepper stepper;
  double field_power_conversion_factor;
  ElectronEnergyAdditionalSourceEvaluator additional_source_evaluator;
  field::CellField<double> mean_energy;
  field::FaceField<double> particle_flux_normal;
  field::CellField<double> additional_source;
  field::CellField<double> source;
  field::CellField<double> increment;

  /** @brief Validates configuration and creates the energy transport operator. */
  [[nodiscard]] static equation::ExplicitElectronEnergyStepper makeStepper(
      const ElectronEnergyConfiguration& configuration,
      const physics::SpeciesCellFields& density, Stepper& transport) {
    if (&configuration.energy_density.mesh() != &density.mesh() ||
        &density.mesh() != &transport.mesh()) {
      throw std::invalid_argument(
          "electron energy and species density must share the simulation mesh");
    }
    if (configuration.electron.value >= transport.species().size()) {
      throw std::invalid_argument("electron energy species id is invalid");
    }
    const auto& properties = transport.species().at(configuration.electron);
    if (properties.charge >= 0.0 || !properties.isTransported()) {
      throw std::invalid_argument(
          "electron energy species must be negative and drift-diffusion "
          "transported");
    }
    if (!configuration.additional_source_evaluator) {
      throw std::invalid_argument(
          "electron energy additional source evaluator must not be empty");
    }
    if (!std::isfinite(configuration.density_floor) ||
        configuration.density_floor < 0.0) {
      throw std::invalid_argument(
          "electron density floor must be finite and non-negative");
    }
    return equation::ExplicitElectronEnergyStepper(
        transport.mesh(), transport.driftVelocityNormal(configuration.electron),
        properties.diffusivity, configuration.boundary_conditions,
        configuration.transport_factor);
  }

  /** @brief Validates one coherent runtime unit system and returns V-to-energy scaling. */
  [[nodiscard]] static double validateMetadata(
      const ElectronEnergyConfiguration& configuration,
      const physics::SpeciesCellFields& density, const Stepper& transport) {
    const auto& declared = transport.fieldMetadata();
    const auto has_quantity = [](const field::FieldMetadata& metadata) {
      return metadata.physical_quantity.has_value();
    };
    bool any_quantity = has_quantity(configuration.energy_density.metadata());
    for (const auto& species_density : density) {
      any_quantity = any_quantity || has_quantity(species_density.metadata());
    }
    any_quantity = any_quantity || has_quantity(declared.number_density) ||
                   has_quantity(declared.number_density_source) ||
                   has_quantity(declared.reaction_rate) ||
                   has_quantity(declared.electric_potential) ||
                   has_quantity(declared.electric_field) ||
                   has_quantity(declared.drift_velocity) ||
                   has_quantity(declared.inverse_time) ||
                   has_quantity(declared.electron_mean_energy) ||
                   has_quantity(declared.electron_energy_density) ||
                   has_quantity(declared.electron_energy_density_source);

    if (!any_quantity) {
      if (!configuration.allow_unitless_raw_values) {
        throw std::invalid_argument(
            "electron energy requires physical metadata or an explicit "
            "unitless raw-value opt-in");
      }
      return 1.0;
    }

    const auto require_kind = [](const field::FieldMetadata& metadata,
                                 unit::QuantityKind kind,
                                 std::string_view name) {
      if (!metadata.physical_quantity.has_value() ||
          metadata.physical_quantity->kind() != kind) {
        throw std::invalid_argument(std::string{name} +
                                    " metadata quantity kind differs");
      }
    };
    const auto require_match = [](const field::FieldMetadata& actual,
                                  const field::FieldMetadata& expected,
                                  std::string_view name) {
      if (!actual.physical_quantity.has_value() ||
          !expected.physical_quantity.has_value() ||
          actual.physical_quantity != expected.physical_quantity) {
        throw std::invalid_argument(std::string{name} +
                                    " metadata unit differs from transport");
      }
    };

    require_kind(declared.number_density,
                 unit::QuantityKind::particle_number_density, "number density");
    require_kind(declared.number_density_source,
                 unit::QuantityKind::particle_number_density_rate,
                 "number density source");
    require_kind(declared.reaction_rate,
                 unit::QuantityKind::reaction_rate_density, "reaction rate");
    require_kind(declared.electric_potential,
                 unit::QuantityKind::electric_potential, "electric potential");
    require_kind(declared.electric_field,
                 unit::QuantityKind::normal_electric_field_strength,
                 "electric field");
    require_kind(declared.drift_velocity,
                 unit::QuantityKind::normal_drift_velocity, "drift velocity");
    require_kind(declared.inverse_time, unit::QuantityKind::frequency,
                 "inverse time");
    require_kind(declared.electron_mean_energy,
                 unit::QuantityKind::electron_mean_energy,
                 "electron mean energy");
    require_kind(declared.electron_energy_density,
                 unit::QuantityKind::electron_energy_density,
                 "electron energy density");
    require_kind(declared.electron_energy_density_source,
                 unit::QuantityKind::electron_energy_density_rate,
                 "electron energy source");

    for (const auto& species_density : density) {
      require_match(species_density.metadata(), declared.number_density,
                    "species density");
    }
    require_match(configuration.energy_density.metadata(),
                  declared.electron_energy_density, "electron energy density");

    const auto number_density_unit =
        declared.number_density.physical_quantity->unit();
    const auto inverse_time_unit =
        declared.inverse_time.physical_quantity->unit();
    const auto potential_unit =
        declared.electric_potential.physical_quantity->unit();
    const auto electric_field_unit =
        declared.electric_field.physical_quantity->unit();
    const auto drift_velocity_unit =
        declared.drift_velocity.physical_quantity->unit();
    const auto mean_energy_unit =
        declared.electron_mean_energy.physical_quantity->unit();
    const auto energy_density_unit =
        declared.electron_energy_density.physical_quantity->unit();

    if (declared.number_density_source.physical_quantity->unit() !=
            number_density_unit * inverse_time_unit ||
        declared.reaction_rate.physical_quantity->unit() !=
            number_density_unit * inverse_time_unit ||
        energy_density_unit != number_density_unit * mean_energy_unit ||
        declared.electron_energy_density_source.physical_quantity->unit() !=
            energy_density_unit * inverse_time_unit) {
      throw std::invalid_argument(
          "plasma field metadata units are not algebraically coherent");
    }

    const auto mesh_length_unit = potential_unit / electric_field_unit;
    if (!mesh_length_unit.is_convertible(units::precise::m) ||
        drift_velocity_unit / mesh_length_unit != inverse_time_unit) {
      throw std::invalid_argument(
          "electric-field and transport metadata use inconsistent length/time "
          "units");
    }

    const double potential_to_volt =
        units::convert(1.0, potential_unit, units::precise::V);
    const double electronvolt_to_storage =
        units::convert(1.0, units::precise::energy::eV, mean_energy_unit);
    const double factor = potential_to_volt * electronvolt_to_storage;
    if (!std::isfinite(factor) || factor <= 0.0) {
      throw std::invalid_argument(
          "electron field-power unit conversion is invalid");
    }
    return factor;
  }

  /** @brief Creates workspaces and validates the initial energy state. */
  ElectronEnergySubsystem(ElectronEnergyConfiguration configuration,
                          const physics::SpeciesCellFields& density,
                          Stepper& transport)
      : energy_density(&configuration.energy_density),
        electron(configuration.electron),
        density_floor(configuration.density_floor),
        stepper(makeStepper(configuration, density, transport)),
        field_power_conversion_factor(
            validateMetadata(configuration, density, transport)),
        additional_source_evaluator(
            std::move(configuration.additional_source_evaluator)),
        mean_energy(transport.mesh(), 0.0,
                    transport.fieldMetadata().electron_mean_energy),
        particle_flux_normal(transport.mesh(), 0.0),
        additional_source(
            transport.mesh(), 0.0,
            transport.fieldMetadata().electron_energy_density_source),
        source(transport.mesh(), 0.0,
               transport.fieldMetadata().electron_energy_density_source),
        increment(transport.mesh(), 0.0,
                  transport.fieldMetadata().electron_energy_density) {
    refreshMeanEnergy(density);
  }

  /** @brief Recomputes mean energy from synchronized density and energy. */
  void refreshMeanEnergy(const physics::SpeciesCellFields& density) {
    physics::computeElectronMeanEnergy(*energy_density, density[electron],
                                       density_floor, mean_energy);
  }

  /** @brief Assembles mandatory field power and additional energy sources. */
  void evaluateSource(const physics::SpeciesCellFields& density,
                      Stepper& transport,
                      const field::CellField<double>& potential,
                      const field::FaceField<double>& electric_field_normal,
                      const physics::ReactionRateFields& reaction_rates,
                      const physics::SpeciesCellFields& species_source) {
    transport.computeParticleFluxNormal(density, electron,
                                        particle_flux_normal);
    physics::computeElectronFieldPowerDensity(particle_flux_normal,
                                              electric_field_normal, source);
    for (double& value : source) {
      value *= field_power_conversion_factor;
    }

    additional_source.fill(std::numeric_limits<double>::quiet_NaN());
    additional_source_evaluator(
        {.density = density,
         .potential = potential,
         .electric_field_normal = electric_field_normal,
         .electron_particle_flux_normal = particle_flux_normal,
         .reaction_rates = reaction_rates,
         .species_source = species_source,
         .mean_energy = mean_energy},
        additional_source);
    for (const double value : additional_source) {
      if (!std::isfinite(value)) {
        throw std::invalid_argument(
            "electron energy additional source evaluator must write finite "
            "values");
      }
    }
    for (mesh::CellId cell = 0; cell < source.mesh().numCells(); ++cell) {
      source[cell] += additional_source[cell];
    }
    for (const double value : source) {
      if (!std::isfinite(value)) {
        throw std::invalid_argument(
            "assembled electron energy source must be finite");
      }
    }
  }

  /** @brief Returns the current transport CFL limit of the energy equation. */
  [[nodiscard]] double maxStableTransportTimeStep() {
    return stepper.maxStableTransportTimeStep();
  }

  /** @brief Returns the current source-aware energy positivity limit. */
  [[nodiscard]] double maxPositiveTimeStep() {
    return stepper.maxPositiveTimeStep(*energy_density, source);
  }

  /** @brief Validates dt and prepares an energy increment without mutation. */
  void prepareIncrement(double dt) {
    stepper.computeStableIncrement(*energy_density, source, dt, increment);
  }

  /** @brief Commits the previously validated energy increment. */
  void commitIncrement() {
    for (mesh::CellId cell = 0; cell < energy_density->mesh().numCells();
         ++cell) {
      (*energy_density)[cell] += increment[cell];
    }
  }
};

template <typename Stepper, typename Clock>
struct PlasmaWorkflowCore {
  physics::SpeciesCellFields* density;
  const physics::ReactionNetwork* reaction_network;
  Stepper* transport_stepper;
  PlasmaReactionRateEvaluator rate_evaluator;
  Clock clock;
  trace::AnyTraceSink trace_sink;
  trace::StatisticsOptions statistics_options;
  FieldOutputOptions field_output_options;
  CheckpointOptions checkpoint_options;
  physics::ReactionRateFields reaction_rates;
  physics::SpeciesCellFields source;
  ElectronEnergySubsystem<Stepper> electron_energy;
  std::unique_ptr<output::dump::ISeries> field_output_series;
  bool field_output_has_snapshots{};
  SimulationState state{SimulationState::ready};
  std::string_view category;
  double workflow_kind;

  PlasmaWorkflowCore(physics::SpeciesCellFields& density_in,
                     const physics::ReactionNetwork& reactions,
                     Stepper& stepper, PlasmaReactionRateEvaluator evaluator,
                     ElectronEnergyConfiguration energy_configuration,
                     Clock clock_in, trace::AnyTraceSink sink,
                     trace::StatisticsOptions statistics,
                     FieldOutputOptions field_output,
                     CheckpointOptions checkpoint,
                     std::string_view trace_category, double checkpoint_kind)
      : density(&density_in),
        reaction_network(&reactions),
        transport_stepper(&stepper),
        rate_evaluator(std::move(evaluator)),
        clock(std::move(clock_in)),
        trace_sink(std::move(sink)),
        statistics_options(statistics),
        field_output_options(std::move(field_output)),
        checkpoint_options(std::move(checkpoint)),
        reaction_rates(density_in.mesh(), reactions.size(), 0.0,
                       stepper.fieldMetadata().reaction_rate),
        source(density_in.mesh(), density_in.size(), 0.0,
               stepper.fieldMetadata().number_density_source),
        electron_energy(std::move(energy_configuration), density_in, stepper),
        category(trace_category),
        workflow_kind(checkpoint_kind) {
    if (!rate_evaluator) {
      throw std::invalid_argument(
          "simulation reaction evaluator must not be empty");
    }
    if (density->size() == 0) {
      throw std::invalid_argument("simulation density must contain species");
    }
    if (density->size() != stepper.species().size()) {
      throw std::invalid_argument(
          "simulation density and transport species counts differ");
    }
    validatePlasmaStatisticsConfiguration(
        statistics_options, *density, stepper.chargeDensity(),
        stepper.potential(), stepper.electricFieldNormal(),
        *electron_energy.energy_density, electron_energy.mean_energy,
        electron_energy.source);
    validateFieldOutputOptions(field_output_options);
    validateCheckpointOptions(checkpoint_options);
  }

  void emitRun(std::string_view name, trace::Severity severity) const noexcept {
    const std::array attributes{
        trace::TraceAttribute{"step", static_cast<std::uint64_t>(clock.step())},
        trace::TraceAttribute{"time", clock.time()},
        trace::TraceAttribute{"end_time", clock.endTime()},
    };
    emitTrace(trace_sink, category, name, severity, attributes);
  }

  void ensureOutput() {
    if (!field_output_options.enabled() || field_output_series) {
      return;
    }
    field_output_series = field_output_options.writer->openSeries(
        {.mesh = &transport_stepper->mesh(),
         .path = plasmaFieldOutputPath(field_output_options),
         .overwrite = field_output_options.overwrite});
  }

  void finishOutput() {
    if (!field_output_series) {
      return;
    }
    if (field_output_has_snapshots) {
      const auto record = field_output_series->finish();
      output::dump::traceCompleted(trace_sink, record);
    }
    field_output_series.reset();
  }

  void writeOutput(bool terminal) {
    if (!field_output_options.enabled()) {
      return;
    }
    if (shouldWritePlasmaFieldSnapshot(field_output_options, clock.step(),
                                       terminal)) {
      (void)writePlasmaFieldSnapshot(
          *density, *transport_stepper, *electron_energy.energy_density,
          electron_energy.mean_energy, field_output_options,
          *field_output_series,
          {.step = static_cast<std::uint64_t>(clock.step()),
           .time = clock.time()});
      field_output_has_snapshots = true;
    }
    if (terminal && field_output_has_snapshots) {
      finishOutput();
    }
  }

  /** @brief Evaluates the energy source from the synchronized workflow state. */
  void evaluateElectronEnergySource() {
    electron_energy.evaluateSource(
        *density, *transport_stepper, transport_stepper->potential(),
        transport_stepper->electricFieldNormal(), reaction_rates, source);
    const std::array attributes{
        trace::TraceAttribute{"step", static_cast<std::uint64_t>(clock.step())},
        trace::TraceAttribute{"field_count", std::uint64_t{1}},
    };
    emitTrace(trace_sink, category, "electron_energy_source.completed",
              trace::Severity::Trace, attributes);
  }

  /** @brief Commits prepared energy and synchronizes its derived mean field. */
  void commitElectronEnergy(double dt) {
    electron_energy.commitIncrement();
    electron_energy.refreshMeanEnergy(*density);
    const std::array attributes{
        trace::TraceAttribute{"step",
                              static_cast<std::uint64_t>(clock.step() + 1)},
        trace::TraceAttribute{"dt", dt},
    };
    emitTrace(trace_sink, category, "electron_energy.completed",
              trace::Severity::Trace, attributes);
  }

  [[nodiscard]] output::OutputRecord save(
      const output::checkpoint::IWriter& writer,
      const std::filesystem::path& path, bool overwrite,
      std::span<const output::checkpoint::ScalarSource> extra_scalars = {})
      const {
    if (state == SimulationState::failed || state == SimulationState::stopped) {
      throw std::logic_error(
          "terminal simulation state cannot be checkpointed");
    }
    auto sources =
        plasmaCheckpointSources(*density, transport_stepper->species());
    sources.push_back(
        {.field = electron_energy.energy_density,
         .key = electronEnergyDensityCheckpointKey(
             electron_energy.electron,
             transport_stepper->species().at(electron_energy.electron).name)});
    const double schema_version = plasma_checkpoint_schema_version;
    std::vector<output::checkpoint::ScalarSource> scalars{
        {.value = &schema_version,
         .key = std::string{plasma_checkpoint_schema_key},
         .metadata = {.name = std::string{plasma_checkpoint_schema_key}}},
        {.value = &workflow_kind,
         .key = std::string{plasma_checkpoint_kind_key},
         .metadata = {.name = std::string{plasma_checkpoint_kind_key}}},
    };
    scalars.insert(scalars.end(), extra_scalars.begin(), extra_scalars.end());
    return writer.write(
        {.mesh = &transport_stepper->mesh(),
         .path = path,
         .stamp = {.step = static_cast<std::uint64_t>(clock.step()),
                   .time = clock.time()},
         .cell_fields = sources,
         .scalars = scalars,
         .overwrite = overwrite});
  }

  void complete() {
    if (state != SimulationState::completed) {
      state = SimulationState::completed;
      emitRun("run.completed", trace::Severity::Info);
      trace_sink.flush();
    }
  }

  void fail() noexcept {
    if (state != SimulationState::failed) {
      state = SimulationState::failed;
      emitRun("run.failed", trace::Severity::Error);
      trace_sink.flush();
    }
  }
};

template <typename Core>
void startWorkflow(Core& core) {
  if (core.state == SimulationState::running) {
    return;
  }
  if (core.state != SimulationState::ready &&
      core.state != SimulationState::paused) {
    throw std::logic_error("simulation cannot start from its current state");
  }
  const bool resumed = core.state == SimulationState::paused;
  if (!resumed && !core.clock.finished()) {
    try {
      core.ensureOutput();
    } catch (...) {
      core.fail();
      throw;
    }
  }
  core.state = SimulationState::running;
  core.emitRun(resumed ? "run.resumed" : "run.started", trace::Severity::Info);
  if (core.clock.finished()) {
    core.complete();
  }
}

template <typename Core>
void pauseWorkflow(Core& core) {
  if (core.state == SimulationState::paused) {
    return;
  }
  if (core.state != SimulationState::running) {
    throw std::logic_error("only a running simulation can be paused");
  }
  core.state = SimulationState::paused;
  core.emitRun("run.paused", trace::Severity::Info);
  core.trace_sink.flush();
}

template <typename Core>
void stopWorkflow(Core& core) {
  if (core.state == SimulationState::stopped) {
    return;
  }
  if (isTerminal(core.state)) {
    throw std::logic_error("terminal simulation cannot be stopped");
  }
  try {
    core.finishOutput();
  } catch (...) {
    core.fail();
    throw;
  }
  core.state = SimulationState::stopped;
  core.emitRun("run.stopped", trace::Severity::Info);
  core.trace_sink.flush();
}

template <typename Core, typename Perform, typename SaveCheckpoint>
[[nodiscard]] linalg::SolverResult advanceWorkflow(
    Core& core, Perform&& perform, SaveCheckpoint&& save_checkpoint) {
  if (core.state == SimulationState::ready) {
    startWorkflow(core);
  }
  if (core.state != SimulationState::running) {
    throw std::logic_error("simulation step requires the running state");
  }
  try {
    const auto result = std::invoke(std::forward<Perform>(perform));
    if (!result.success()) {
      core.fail();
      return result;
    }
    // Flush cadence is a workflow policy. Keeping it here prevents generic
    // sinks and channel routers from parsing simulation-specific event names.
    if (!core.clock.finished() &&
        core.clock.step() % trace_flush_step_interval == 0) {
      core.trace_sink.flush();
    }
    if (shouldWriteCheckpoint(core.checkpoint_options, core.clock.step(),
                              core.clock.finished())) {
      std::invoke(std::forward<SaveCheckpoint>(save_checkpoint));
    }
    if (core.clock.finished()) {
      core.complete();
    }
    return result;
  } catch (...) {
    core.fail();
    throw;
  }
}

template <typename Core, typename Advance>
void runWorkflow(Core& core, Advance&& advance,
                 std::string_view failure_message) {
  if (core.state == SimulationState::ready ||
      core.state == SimulationState::paused) {
    startWorkflow(core);
  } else if (core.state == SimulationState::completed) {
    return;
  } else if (core.state != SimulationState::running) {
    throw std::logic_error("simulation cannot run from its current state");
  }
  while (core.state == SimulationState::running) {
    if (!std::invoke(advance).success()) {
      throw std::runtime_error(std::string{failure_message});
    }
  }
}

}  // namespace pemu::simulation::detail
