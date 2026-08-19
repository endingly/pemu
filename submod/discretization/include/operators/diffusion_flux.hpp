#pragma once

#include <pemu/boundary/boundary_condition.hpp>
#include <pemu/boundary/boundary_condition_set.hpp>
#include <pemu/field/types.hpp>
#include <pemu/field/utils.hpp>
#include <pemu/mesh/geometry.hpp>
#include <pemu/mesh/i_mesh.hpp>

#include <concepts>
#include <stdexcept>
#include <type_traits>
#include <variant>

namespace pemu::discretization::operators {

template <field::CellFieldLike CellField, field::FaceFieldLike FluxField>
  requires std::same_as<typename CellField::value_type,
                        typename FluxField::value_type> &&
           std::floating_point<typename CellField::value_type>
void diffusionFlux(const CellField& u, double diffusivity,
                   const boundary::BoundaryConditionSet& bc, FluxField& flux) {
  field::ensureSameMesh(u, flux);

  if (diffusivity < 0.0) {
    throw std::invalid_argument("diffusivity must be non-negative");
  }

  const auto& mesh = u.mesh();

  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {

    const auto owner = mesh.owner(face);

    const auto normal = mesh.faceNormal(face);

    // ====================================================
    // Internal face
    //
    //      u_N - u_P
    //
    // grad(u)·n ≈ -----------
    //              d_PN
    //
    // Gamma·n =
    //
    //      -D grad(u)·n
    // ====================================================

    if (!mesh.isBoundary(face)) {

      const auto neighbor = mesh.neighbor(face);

      const auto xp = mesh.cellCenter(owner);

      const auto xn = mesh.cellCenter(neighbor);

      const double distance = mesh::dot(xn - xp, normal);

      if (distance <= 0.0) {
        throw std::runtime_error(
            "invalid owner-neighbor "
            "normal distance");
      }

      const double gradient_normal = (u[neighbor] - u[owner]) / distance;

      flux[face] = -diffusivity * gradient_normal;

      continue;
    }

    // ====================================================
    // Boundary face
    // ====================================================

    const auto boundary_id = mesh.boundaryId(face);

    if (boundary_id == mesh::invalid_boundary) {

      throw std::runtime_error("boundary face has no BoundaryId");
    }

    const auto& condition = bc.at(boundary_id);

    std::visit(
        [&](const auto& condition_value) {
          using BcType = std::remove_cvref_t<decltype(condition_value)>;

          // --------------------------------------------
          // Dirichlet
          //
          //              u_b - u_P
          // Gamma·n = -D -----------
          //                 d_Pb
          // --------------------------------------------

          if constexpr (std::same_as<BcType, boundary::Dirichlet>) {

            const auto xp = mesh.cellCenter(owner);

            const auto xf = mesh.faceCenter(face);

            const double distance = mesh::dot(xf - xp, normal);

            if (distance <= 0.0) {
              throw std::runtime_error(
                  "invalid cell-boundary "
                  "normal distance");
            }

            const double gradient_normal =
                (condition_value.value - u[owner]) / distance;

            flux[face] = -diffusivity * gradient_normal;
          }

          // --------------------------------------------
          // Neumann
          //
          // Boundary condition already specifies:
          //
          //     -D grad(u) · n
          //
          // i.e. exactly the quantity stored by
          // FaceField.
          // --------------------------------------------

          else if constexpr (std::same_as<BcType, boundary::Neumann>) {

            flux[face] = condition_value.value;
          }
        },
        condition);
  }
}

}  // namespace pemu::discretization::operators