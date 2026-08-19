#include <pemu/equation/explicit_species_continuity_stepper.hpp>

namespace pemu::equation {

double ExplicitSpeciesContinuityStepper::maxTransportCfl() const {
  field::CellField<double> diagonal(*mesh_, 0.0);

  for (mesh::FaceId face = 0; face < mesh_->numFaces(); ++face) {

    const auto owner = mesh_->owner(face);

    const auto normal = mesh_->faceNormal(face);

    const double vn = (*normal_drift_velocity_)[face];

    // ====================================================
    // Internal face
    // ====================================================

    if (!mesh_->isBoundary(face)) {

      const auto neighbor = mesh_->neighbor(face);

      const double distance = mesh::dot(
          mesh_->cellCenter(neighbor) - mesh_->cellCenter(owner), normal);

      if (distance <= 0.0) {
        throw std::runtime_error("invalid internal face distance");
      }

      const double pe = vn * distance / diffusivity_;

      const double scale = diffusivity_ * mesh_->faceArea(face) / distance;

      //
      // Owner diagonal.
      //
      diagonal[owner] += scale * discretization::operators::bernoulli(-pe);

      //
      // Neighbor diagonal.
      //
      diagonal[neighbor] += scale * discretization::operators::bernoulli(pe);

      continue;
    }

    // ====================================================
    // Boundary Dirichlet face
    //
    // owner diagonal coefficient:
    //
    //     D A / d * B(-Pe)
    // ====================================================

    const double distance =
        mesh::dot(mesh_->faceCenter(face) - mesh_->cellCenter(owner), normal);

    if (distance <= 0.0) {
      throw std::runtime_error("invalid boundary distance");
    }

    const double pe = vn * distance / diffusivity_;

    const double scale = diffusivity_ * mesh_->faceArea(face) / distance;

    diagonal[owner] += scale * discretization::operators::bernoulli(-pe);
  }

  double maximum = 0.0;

  for (mesh::CellId cell = 0; cell < mesh_->numCells(); ++cell) {

    const double local = dt_ * diagonal[cell] / mesh_->cellVolume(cell);

    maximum = std::max(maximum, local);
  }

  return maximum;
}

};  // namespace pemu::equation