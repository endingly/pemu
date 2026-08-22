#pragma once

#include <pemu/equation/fixed_step_multi_species_drift_diffusion_stepper.hpp>
#include <pemu/output/checkpoint/i_reader.hpp>
#include <pemu/physics/reaction.hpp>
#include <pemu/physics/species.hpp>
#include <pemu/simulation/checkpoint.hpp>
#include <pemu/simulation/field_output.hpp>
#include <pemu/simulation/fixed_step_clock.hpp>
#include <pemu/simulation/plasma_reaction_rate_evaluator.hpp>
#include <pemu/simulation/simulation_state.hpp>
#include <pemu/trace/any_trace_sink.hpp>
#include <pemu/trace/statistics.hpp>

#include <cstddef>
#include <filesystem>
#include <memory>

namespace pemu::simulation {

/** @brief Stateful workflow for fixed-step multi-species plasma simulation. */
class FixedStepPlasmaSimulation {
 public:
  using Stepper = equation::FixedStepMultiSpeciesDriftDiffusionStepper;

  FixedStepPlasmaSimulation(physics::SpeciesCellFields& density,
                            const physics::ReactionNetwork& reaction_network,
                            Stepper& transport_stepper,
                            PlasmaReactionRateEvaluator rate_evaluator,
                            FixedStepClock clock,
                            trace::AnyTraceSink trace_sink = {},
                            trace::StatisticsOptions statistics_options = {},
                            FieldOutputOptions field_output_options = {},
                            CheckpointOptions checkpoint_options = {});
  ~FixedStepPlasmaSimulation();

  FixedStepPlasmaSimulation(const FixedStepPlasmaSimulation&) = delete;
  FixedStepPlasmaSimulation& operator=(const FixedStepPlasmaSimulation&) =
      delete;
  FixedStepPlasmaSimulation(FixedStepPlasmaSimulation&&) noexcept;
  FixedStepPlasmaSimulation& operator=(FixedStepPlasmaSimulation&&) noexcept;

  void start();
  void pause();
  void stop();
  [[nodiscard]] linalg::SolverResult advance();
  void run();

  [[nodiscard]] SimulationState state() const noexcept;
  [[nodiscard]] double time() const noexcept;
  [[nodiscard]] double timeStep() const noexcept;
  [[nodiscard]] std::size_t step() const noexcept;
  [[nodiscard]] std::size_t totalSteps() const noexcept;
  [[nodiscard]] bool finished() const noexcept;

  [[nodiscard]] output::OutputRecord saveCheckpoint() const;
  [[nodiscard]] output::OutputRecord saveCheckpoint(
      const output::checkpoint::IWriter& writer,
      const std::filesystem::path& path, bool overwrite = false) const;
  [[nodiscard]] output::OutputRecord restoreCheckpoint(
      const output::checkpoint::IReader& reader,
      const std::filesystem::path& path);

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace pemu::simulation
