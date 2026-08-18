#pragma once

#include "types.hpp"

#include <cstddef>
#include <span>

namespace pemu::mesh {

class IMesh {
 public:
  virtual ~IMesh() = default;

  // --------------------------------------------------------
  // Basic information
  // --------------------------------------------------------

  [[nodiscard]]
  virtual int dimension() const noexcept = 0;

  [[nodiscard]]
  virtual std::size_t numCells() const noexcept = 0;

  [[nodiscard]]
  virtual std::size_t numFaces() const noexcept = 0;

  [[nodiscard]]
  virtual std::size_t numVertices() const noexcept = 0;

  // --------------------------------------------------------
  // Topology
  // --------------------------------------------------------

  [[nodiscard]]
  virtual CellId owner(FaceId face) const = 0;

  [[nodiscard]]
  virtual CellId neighbor(FaceId face) const = 0;

  [[nodiscard]]
  virtual bool isBoundary(FaceId face) const = 0;

  [[nodiscard]]
  virtual std::span<const FaceId> cellFaces(CellId cell) const = 0;

  // --------------------------------------------------------
  // Geometry
  // --------------------------------------------------------

  [[nodiscard]]
  virtual Vec3 cellCenter(CellId cell) const = 0;

  [[nodiscard]]
  virtual double cellVolume(CellId cell) const = 0;

  [[nodiscard]]
  virtual Vec3 faceCenter(FaceId face) const = 0;

  [[nodiscard]]
  virtual double faceArea(FaceId face) const = 0;

  //
  // Convention:
  //
  // faceNormal(face) always points:
  //
  //     owner -> neighbor
  //
  // or outward if boundary face.
  //
  [[nodiscard]]
  virtual Vec3 faceNormal(FaceId face) const = 0;

  // --------------------------------------------------------
  // Boundary metadata
  // --------------------------------------------------------

  [[nodiscard]]
  virtual BoundaryId boundaryId(FaceId face) const = 0;
};

}  // namespace pemu::mesh