#pragma once

#include <pemu/boundary/boundary_condition.hpp>
#include <pemu/boundary/boundary_condition_set.hpp>
#include <pemu/discretization/operators/bernoulli.hpp>
#include <pemu/field/types.hpp>
#include <pemu/field/utils.hpp>
#include <pemu/mesh/geometry.hpp>
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
void scharfetterGummelFlux(const StateField& state,
                           const VelocityField& normal_velocity,
                           const double diffusivity,
                           const boundary::BoundaryConditionSet& bc,
                           FluxField& flux) {
  using Scalar = typename StateField::value_type;

  field::ensureSameMesh(state, normal_velocity);

  field::ensureSameMesh(state, flux);

  if (diffusivity <= 0.0) {
    throw std::invalid_argument(
        "Scharfetter-Gummel requires "
        "positive diffusivity");
  }

  const auto& mesh = state.mesh();

  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {

    const Scalar vn = normal_velocity[face];

    const auto owner = mesh.owner(face);

    const auto normal = mesh.faceNormal(face);

    // ====================================================
    // Internal face
    // ====================================================

    if (!mesh.isBoundary(face)) {

      const auto neighbor = mesh.neighbor(face);

      const auto delta = mesh.cellCenter(neighbor) - mesh.cellCenter(owner);

      const double distance = mesh::dot(delta, normal);

      if (distance <= 0.0) {
        throw std::runtime_error(
            "invalid owner-neighbor "
            "normal distance");
      }

      const double pe = static_cast<double>(vn) * distance / diffusivity;

      const double bp = bernoulli(pe);

      const double bm = bernoulli(-pe);

      flux[face] = static_cast<Scalar>(
          diffusivity / distance * (bm * state[owner] - bp * state[neighbor]));

      continue;
    }

    // ====================================================
    // Boundary face
    //
    // First SG version deliberately supports only
    // Dirichlet boundary states.
    // ====================================================

    const auto boundary_id = mesh.boundaryId(face);

    if (boundary_id == mesh::invalid_boundary) {

      throw std::runtime_error("boundary face has no BoundaryId");
    }

    if (!bc.contains(boundary_id)) {

      throw std::runtime_error(
          "Scharfetter-Gummel boundary "
          "condition is missing");
    }

    const auto& condition = bc.at(boundary_id);

    const auto* dirichlet = std::get_if<boundary::Dirichlet>(&condition);

    if (dirichlet == nullptr) {

      throw std::runtime_error(
          "Scharfetter-Gummel currently "
          "requires Dirichlet boundary "
          "conditions");
    }

    // ----------------------------------------------------
    // owner -> boundary
    // ----------------------------------------------------

    const auto delta = mesh.faceCenter(face) - mesh.cellCenter(owner);

    const double distance = mesh::dot(delta, normal);

    if (distance <= 0.0) {
      throw std::runtime_error(
          "invalid owner-boundary "
          "normal distance");
    }

    const double pe = static_cast<double>(vn) * distance / diffusivity;

    const double bp = bernoulli(pe);

    const double bm = bernoulli(-pe);

    flux[face] = static_cast<Scalar>(
        diffusivity / distance * (bm * state[owner] - bp * dirichlet->value));
  }
}

}  // namespace pemu::discretization::operators