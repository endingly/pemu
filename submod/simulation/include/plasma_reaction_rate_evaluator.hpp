#pragma once

#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>
#include <pemu/physics/reaction.hpp>
#include <pemu/physics/species.hpp>

#include <functional>

namespace pemu::simulation {

using PlasmaReactionRateEvaluator = std::function<void(
    const physics::SpeciesCellFields&, const field::CellField<double>&,
    const field::FaceField<double>&, physics::ReactionRateFields&)>;

}  // namespace pemu::simulation
