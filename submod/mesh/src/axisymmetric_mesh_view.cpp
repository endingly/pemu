#include <pemu/mesh/axisymmetric_mesh_view.hpp>

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>

namespace pemu::mesh {

namespace {

// Verifies that a point can be interpreted as a finite `(r, z)` coordinate.
void validateMeridionalPoint(const Vec3 point, const char* entity_name) {
  if (!std::isfinite(point.x) || !std::isfinite(point.y) ||
      !std::isfinite(point.z)) {
    throw std::invalid_argument(std::string{"axisymmetric "} + entity_name +
                                " coordinates must be finite");
  }
  if (point.x < 0.0) {
    throw std::invalid_argument(std::string{"axisymmetric "} + entity_name +
                                " radius must be non-negative");
  }
}

// Revolves a planar measure whose centroid lies at `radius` around the axis.
[[nodiscard]] double revolvedMeasure(const double planar_measure,
                                     const double radius,
                                     const char* entity_name,
                                     const bool may_be_zero) {
  if (!std::isfinite(planar_measure) || planar_measure <= 0.0) {
    throw std::invalid_argument(std::string{"axisymmetric base "} +
                                entity_name + " measure must be positive");
  }

  const double measure = 2.0 * std::numbers::pi * radius * planar_measure;
  if (!std::isfinite(measure) || (!may_be_zero && measure <= 0.0)) {
    throw std::invalid_argument(std::string{"axisymmetric "} + entity_name +
                                " measure must be finite and valid");
  }
  return measure;
}

}  // namespace

AxisymmetricMeshView::AxisymmetricMeshView(const IMesh& mesh) : mesh_(&mesh) {
  if (mesh.dimension() != 2) {
    throw std::invalid_argument(
        "axisymmetric mesh view requires a two-dimensional base mesh");
  }

  for (VertexId vertex_id = 0; vertex_id < mesh.numVertices(); ++vertex_id) {
    validateMeridionalPoint(mesh.vertex(vertex_id), "vertex");
  }

  cell_volumes_.reserve(mesh.numCells());
  for (CellId cell = 0; cell < mesh.numCells(); ++cell) {
    const Vec3 center = mesh.cellCenter(cell);
    validateMeridionalPoint(center, "cell-center");
    cell_volumes_.push_back(
        revolvedMeasure(mesh.cellVolume(cell), center.x, "cell", false));
  }

  face_areas_.reserve(mesh.numFaces());
  for (FaceId face = 0; face < mesh.numFaces(); ++face) {
    const Vec3 center = mesh.faceCenter(face);
    validateMeridionalPoint(center, "face-center");
    face_areas_.push_back(
        revolvedMeasure(mesh.faceArea(face), center.x, "face", true));
  }
}

int AxisymmetricMeshView::dimension() const noexcept {
  return mesh_->dimension();
}

std::size_t AxisymmetricMeshView::numCells() const noexcept {
  return mesh_->numCells();
}

std::size_t AxisymmetricMeshView::numFaces() const noexcept {
  return mesh_->numFaces();
}

std::size_t AxisymmetricMeshView::numVertices() const noexcept {
  return mesh_->numVertices();
}

CellId AxisymmetricMeshView::owner(const FaceId face) const {
  return mesh_->owner(face);
}

CellId AxisymmetricMeshView::neighbor(const FaceId face) const {
  return mesh_->neighbor(face);
}

bool AxisymmetricMeshView::isBoundary(const FaceId face) const {
  return mesh_->isBoundary(face);
}

std::span<const FaceId> AxisymmetricMeshView::cellFaces(
    const CellId cell) const {
  return mesh_->cellFaces(cell);
}

std::span<const VertexId> AxisymmetricMeshView::cellVertices(
    const CellId cell) const {
  return mesh_->cellVertices(cell);
}

Vec3 AxisymmetricMeshView::cellCenter(const CellId cell) const {
  return mesh_->cellCenter(cell);
}

Vec3 AxisymmetricMeshView::vertex(const VertexId vertex_id) const {
  return mesh_->vertex(vertex_id);
}

double AxisymmetricMeshView::cellVolume(const CellId cell) const {
  return cell_volumes_.at(cell);
}

Vec3 AxisymmetricMeshView::faceCenter(const FaceId face) const {
  return mesh_->faceCenter(face);
}

double AxisymmetricMeshView::faceArea(const FaceId face) const {
  return face_areas_.at(face);
}

Vec3 AxisymmetricMeshView::faceNormal(const FaceId face) const {
  return mesh_->faceNormal(face);
}

BoundaryId AxisymmetricMeshView::boundaryId(const FaceId face) const {
  return mesh_->boundaryId(face);
}

}  // namespace pemu::mesh
