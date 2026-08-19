#pragma once

#include <pemu/boundary/boundary_condition_set.hpp>
#include <pemu/discretization/operators/diffusion_flux.hpp>
#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>

#include <stdexcept>

namespace pemu::discretization::operators {

inline void electricFieldNormal(
    const field::CellField<double>& potential, double permittivity,
    const boundary::BoundaryConditionSet& potential_bc,
    field::FaceField<double>& normal_electric_field) {
  if (permittivity <= 0.0) {
    throw std::invalid_argument("permittivity must be positive");
  }

  field::ensureSameMesh(potential, normal_electric_field);

  //
  // diffusionFlux computes:
  //
  //     q_f = -epsilon grad(phi) · n
  //
  // which is:
  //
  //     q_f = epsilon E_n
  //
  discretization::operators::diffusionFlux(potential, permittivity,
                                           potential_bc, normal_electric_field);

  //
  // Convert electric displacement-like flux
  // to electric field:
  //
  //     E_n = q_f / epsilon
  //
  for (mesh::FaceId face = 0; face < potential.mesh().numFaces(); ++face) {

    normal_electric_field[face] /= permittivity;
  }
}

}  // namespace pemu::discretization::operators