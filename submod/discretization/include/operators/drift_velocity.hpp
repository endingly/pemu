#pragma once

#include <pemu/field/face_field.hpp>
#include <pemu/field/utils.hpp>

#include <pemu/physics/charge_polarity.hpp>

#include <stdexcept>

namespace pemu::discretization::operators {

inline void driftVelocityNormal(
    const field::FaceField<double>& normal_electric_field, double mobility,
    physics::ChargePolarity polarity,
    field::FaceField<double>& normal_drift_velocity) {
  field::ensureSameMesh(normal_electric_field, normal_drift_velocity);

  if (mobility < 0.0) {
    throw std::invalid_argument("mobility must be non-negative");
  }

  const double sign = physics::polaritySign(polarity);

  for (mesh::FaceId face = 0; face < normal_electric_field.mesh().numFaces();
       ++face) {

    normal_drift_velocity[face] = sign * mobility * normal_electric_field[face];
  }
}

}  // namespace pemu::discretization::operators