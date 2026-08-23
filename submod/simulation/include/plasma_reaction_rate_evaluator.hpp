#pragma once

#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>
#include <pemu/physics/reaction/types.hpp>
#include <pemu/physics/species.hpp>

#include <functional>

namespace pemu::simulation {

/** @brief Read-only synchronized plasma state used to evaluate reaction rates. */
struct PlasmaReactionRateContext {
  const physics::SpeciesCellFields& density;
  const field::CellField<double>& potential;
  const field::FaceField<double>& electric_field_normal;
  const field::CellField<double>& electron_mean_energy;
};

/**
 * @brief Overwrites reaction-rate fields from one synchronized plasma state.
 *
 * The context form keeps the evaluator contract stable as chemistry gains
 * state variables and avoids positional callbacks with multiple similar field
 * types. Implementations must write every reaction field they own; Simulation
 * clears the complete set before invoking the callback.
 */
using PlasmaReactionRateEvaluator = std::function<void(
    const PlasmaReactionRateContext&, physics::reaction::ReactionRateFields&)>;

}  // namespace pemu::simulation
