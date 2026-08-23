#pragma once

#include <pemu/mesh/i_mesh.hpp>

#include <vector>

namespace pemu::mesh {

/**
 * @brief Interprets a planar `(r, z)` mesh as an axisymmetric control-volume
 * mesh.
 *
 * The wrapped mesh uses `x` as the non-negative radius and `y` as the axial
 * coordinate. Topology, coordinates, centers and meridional unit normals are
 * preserved. Cell volumes and face areas are replaced by their surfaces-of-
 * revolution measures around `r = 0`.
 *
 * This is a non-owning view. The wrapped mesh must outlive this object and must
 * not change its geometry after construction.
 */
class AxisymmetricMeshView final : public IMesh {
 public:
  /**
   * @brief Validates and precomputes axisymmetric measures for a planar mesh.
   * @param mesh Two-dimensional mesh whose `x` coordinates are radii.
   * @throws std::invalid_argument if the geometry cannot represent an
   * axisymmetric meridional section.
   */
  explicit AxisymmetricMeshView(const IMesh& mesh);

  /** @brief Returns the meridional mesh dimension, always two. */
  [[nodiscard]] int dimension() const noexcept override;

  /** @brief Returns the wrapped mesh cell count. */
  [[nodiscard]] std::size_t numCells() const noexcept override;

  /** @brief Returns the wrapped mesh face count. */
  [[nodiscard]] std::size_t numFaces() const noexcept override;

  /** @brief Returns the wrapped mesh vertex count. */
  [[nodiscard]] std::size_t numVertices() const noexcept override;

  /** @brief Returns the wrapped owner cell of a face. */
  [[nodiscard]] CellId owner(FaceId face) const override;

  /** @brief Returns the wrapped neighbor cell of a face. */
  [[nodiscard]] CellId neighbor(FaceId face) const override;

  /** @brief Reports whether a wrapped face lies on the domain boundary. */
  [[nodiscard]] bool isBoundary(FaceId face) const override;

  /** @brief Returns the wrapped cell-to-face connectivity. */
  [[nodiscard]] std::span<const FaceId> cellFaces(CellId cell) const override;

  /** @brief Returns the wrapped cell-to-vertex connectivity. */
  [[nodiscard]] std::span<const VertexId> cellVertices(
      CellId cell) const override;

  /** @brief Returns the wrapped meridional cell center `(r, z)`. */
  [[nodiscard]] Vec3 cellCenter(CellId cell) const override;

  /** @brief Returns a wrapped meridional vertex `(r, z)`. */
  [[nodiscard]] Vec3 vertex(VertexId vertex) const override;

  /** @brief Returns the precomputed volume obtained by revolving a cell. */
  [[nodiscard]] double cellVolume(CellId cell) const override;

  /** @brief Returns the wrapped meridional face center `(r, z)`. */
  [[nodiscard]] Vec3 faceCenter(FaceId face) const override;

  /** @brief Returns the precomputed area obtained by revolving a face. */
  [[nodiscard]] double faceArea(FaceId face) const override;

  /** @brief Returns the wrapped unit normal in the meridional plane. */
  [[nodiscard]] Vec3 faceNormal(FaceId face) const override;

  /** @brief Returns the wrapped boundary identifier of a face. */
  [[nodiscard]] BoundaryId boundaryId(FaceId face) const override;

 private:
  const IMesh* mesh_;
  std::vector<double> cell_volumes_;
  std::vector<double> face_areas_;
};

}  // namespace pemu::mesh
