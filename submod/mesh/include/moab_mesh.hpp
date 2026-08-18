#pragma once

#include <pemu/mesh/i_mesh.hpp>

#include <moab/Core.hpp>

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace pemu::mesh {

class MoabMesh final : public IMesh {
 public:
  explicit MoabMesh(const std::string& filename);

  ~MoabMesh() override = default;

  MoabMesh(const MoabMesh&) = delete;

  MoabMesh& operator=(const MoabMesh&) = delete;

  MoabMesh(MoabMesh&&) noexcept = default;

  MoabMesh& operator=(MoabMesh&&) noexcept = default;

  // --------------------------------------------------------
  // Basic information
  // --------------------------------------------------------

  int dimension() const noexcept override;

  std::size_t numCells() const noexcept override;

  std::size_t numFaces() const noexcept override;

  std::size_t numVertices() const noexcept override;

  // --------------------------------------------------------
  // Topology
  // --------------------------------------------------------

  CellId owner(FaceId face) const override;

  CellId neighbor(FaceId face) const override;

  bool isBoundary(FaceId face) const override;

  std::span<const FaceId> cellFaces(CellId cell) const override;

  // --------------------------------------------------------
  // Geometry
  // --------------------------------------------------------

  Vec3 cellCenter(CellId cell) const override;

  double cellVolume(CellId cell) const override;

  Vec3 faceCenter(FaceId face) const override;

  double faceArea(FaceId face) const override;

  Vec3 faceNormal(FaceId face) const override;

  // --------------------------------------------------------
  // Boundary
  // --------------------------------------------------------

  BoundaryId boundaryId(FaceId face) const override;

 private:
  void load(const std::string& filename);

  void buildEntities();

  void buildHandleMaps();

  void buildTopology();

  void buildGeometry();

  void buildBoundaryMetadata();

 private:
  std::unique_ptr<moab::Core> core_;

  int dimension_{0};

  // --------------------------------------------------------
  // Original MOAB entity handles
  //
  // Only backend construction uses them.
  // Runtime hot paths use dense IDs below.
  // --------------------------------------------------------

  std::vector<moab::EntityHandle> cell_handles_;

  std::vector<moab::EntityHandle> face_handles_;

  std::vector<moab::EntityHandle> vertex_handles_;

  // --------------------------------------------------------
  // MOAB handle -> dense pemu ID
  // --------------------------------------------------------

  std::unordered_map<moab::EntityHandle, CellId> cell_id_;

  std::unordered_map<moab::EntityHandle, FaceId> face_id_;

  // --------------------------------------------------------
  // Topology
  // --------------------------------------------------------

  std::vector<CellId> face_owner_;

  std::vector<CellId> face_neighbor_;

  //
  // CSR-like cell -> face connectivity.
  //
  std::vector<std::uint32_t> cell_face_offsets_;

  std::vector<FaceId> cell_faces_;

  // --------------------------------------------------------
  // Geometry
  // --------------------------------------------------------

  std::vector<Vec3> cell_centers_;

  //
  // For 2D mesh this stores cell area.
  //
  std::vector<double> cell_volumes_;

  std::vector<Vec3> face_centers_;

  //
  // Unit normal:
  //
  // internal: owner -> neighbor
  // boundary: owner -> outside
  //
  std::vector<Vec3> face_normals_;

  //
  // For 2D mesh this stores edge length.
  //
  std::vector<double> face_areas_;

  // --------------------------------------------------------
  // Boundary
  // --------------------------------------------------------

  std::vector<BoundaryId> boundary_ids_;
};

}  // namespace pemu::mesh