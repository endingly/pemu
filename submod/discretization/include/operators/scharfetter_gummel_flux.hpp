#pragma once

#include <pemu/boundary/boundary_condition.hpp>
#include <pemu/boundary/boundary_condition_set.hpp>
#include <pemu/discretization/operators/bernoulli.hpp>
#include <pemu/discretization/operators/linear_boundary_flux.hpp>
#include <pemu/field/types.hpp>
#include <pemu/field/utils.hpp>
#include <pemu/mesh/geometry.hpp>
#include <pemu/mesh/i_mesh.hpp>

#include <cmath>
#include <concepts>
#include <stdexcept>
#include <type_traits>
#include <variant>

namespace pemu::discretization::operators {
namespace detail {

/** @brief Implements SG transport with an optional linear boundary override. */
template <field::CellFieldLike StateField, field::FaceFieldLike VelocityField,
          field::FaceFieldLike FluxField>
  requires std::same_as<typename StateField::value_type,
                        typename VelocityField::value_type> &&
           std::same_as<typename StateField::value_type,
                        typename FluxField::value_type> &&
           std::floating_point<typename StateField::value_type>
void scharfetterGummelFluxImpl(
    const StateField& state, const VelocityField& normal_velocity,
    const double diffusivity, const boundary::BoundaryConditionSet& bc,
    const LinearBoundaryFluxView* linear_boundary_flux, FluxField& flux) {
  using Scalar = typename StateField::value_type;

  field::ensureSameMesh(state, normal_velocity);

  field::ensureSameMesh(state, flux);

  if (diffusivity <= 0.0) {
    throw std::invalid_argument(
        "Scharfetter-Gummel requires "
        "positive diffusivity");
  }

  const auto& mesh = state.mesh();

  if (linear_boundary_flux != nullptr) {
    linear_boundary_flux->validate(state);
  }

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
    // Faces without a linear override use scalar boundary conditions.
    // ====================================================

    if (linear_boundary_flux != nullptr && linear_boundary_flux->active(face)) {
      flux[face] =
          static_cast<Scalar>(linear_boundary_flux->normalFlux(state, face));
      continue;
    }

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

    if (const auto* neumann = std::get_if<boundary::Neumann>(&condition)) {
      // Neumann prescribes the outward diffusive part -D grad(n).n.
      // The owner-state drift contribution remains part of the total flux.
      if (!std::isfinite(static_cast<double>(vn)) ||
          !std::isfinite(static_cast<double>(state[owner])) ||
          !std::isfinite(neumann->value)) {
        throw std::invalid_argument(
            "Scharfetter-Gummel Neumann flux requires finite inputs");
      }
      const double candidate =
          static_cast<double>(vn) * static_cast<double>(state[owner]) +
          neumann->value;
      if (!std::isfinite(candidate)) {
        throw std::overflow_error("Scharfetter-Gummel Neumann flux overflowed");
      }
      flux[face] = static_cast<Scalar>(candidate);
      continue;
    }

    const auto* dirichlet = std::get_if<boundary::Dirichlet>(&condition);

    if (dirichlet == nullptr) {
      throw std::logic_error("unsupported Scharfetter-Gummel boundary type");
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

}  // namespace detail

/** @brief Computes SG flux with scalar boundary conditions only. */
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
  detail::scharfetterGummelFluxImpl(state, normal_velocity, diffusivity, bc,
                                    nullptr, flux);
}

/** @brief Computes SG interior flux and linear wall flux on selected faces. */
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
                           const LinearBoundaryFluxView& linear_boundary_flux,
                           FluxField& flux) {
  detail::scharfetterGummelFluxImpl(state, normal_velocity, diffusivity, bc,
                                    &linear_boundary_flux, flux);
}

}  // namespace pemu::discretization::operators
