#pragma once

#include <pemu/boundary/boundary_condition.hpp>
#include <pemu/boundary/boundary_condition_set.hpp>
#include <pemu/field/types.hpp>
#include <pemu/field/utils.hpp>
#include <pemu/mesh/i_mesh.hpp>

#include <concepts>
#include <stdexcept>
#include <type_traits>
#include <variant>

namespace pemu::discretization::operators {

template <field::CellFieldLike StateField, field::FaceFieldLike VelocityField,
          field::FaceFieldLike FluxField>
  requires std::same_as<typename StateField::value_type,
                        typename VelocityField::value_type> &&
           std::same_as<typename StateField::value_type,
                        typename FluxField::value_type> &&
           std::floating_point<typename StateField::value_type>
void upwindAdvectionFlux(const StateField& state,
                         const VelocityField& normal_velocity,
                         const boundary::BoundaryConditionSet& bc,
                         FluxField& flux) {
  using Scalar = typename StateField::value_type;

  field::ensureSameMesh(state, normal_velocity);

  field::ensureSameMesh(state, flux);

  const auto& mesh = state.mesh();

  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {

    const Scalar vn = normal_velocity[face];

    const auto owner = mesh.owner(face);

    // ====================================================
    // Internal face
    //
    // vn > 0:
    //
    //       P ------> N
    //
    // upwind = P
    //
    // vn < 0:
    //
    //       P <------ N
    //
    // upwind = N
    // ====================================================

    if (!mesh.isBoundary(face)) {

      const auto neighbor = mesh.neighbor(face);

      const Scalar upwind_value =
          vn >= Scalar{0} ? state[owner] : state[neighbor];

      flux[face] = vn * upwind_value;

      continue;
    }

    // ====================================================
    // Boundary face
    //
    // n points owner -> outside.
    //
    // vn >= 0:
    //     outflow
    //
    // vn < 0:
    //     inflow
    // ====================================================

    // ----------------------------------------------------
    // Outflow:
    //
    // Domain value itself is upwind.
    // No external boundary value is required.
    // ----------------------------------------------------

    if (vn >= Scalar{0}) {

      flux[face] = vn * state[owner];

      continue;
    }

    // ----------------------------------------------------
    // Inflow:
    //
    // Need externally prescribed state.
    // ----------------------------------------------------

    const auto boundary_id = mesh.boundaryId(face);

    if (boundary_id == mesh::invalid_boundary) {

      throw std::runtime_error("inflow face has no BoundaryId");
    }

    if (!bc.contains(boundary_id)) {

      throw std::runtime_error(
          "inflow boundary has no "
          "boundary condition");
    }

    const auto& condition = bc.at(boundary_id);

    const auto* dirichlet = std::get_if<boundary::Dirichlet>(&condition);

    if (dirichlet == nullptr) {

      throw std::runtime_error(
          "advection inflow boundary "
          "requires Dirichlet condition");
    }

    flux[face] = vn * static_cast<Scalar>(dirichlet->value);
  }
}

}  // namespace pemu::discretization::operators