#include <pemu/mesh/moab_mesh.hpp>
#include <pemu/mesh/geometry.hpp>

#include <MBTagConventions.hpp>
#include <moab/CN.hpp>
#include <moab/Interface.hpp>
#include <moab/Range.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace pemu::mesh {

namespace {

// ============================================================
// MOAB error handling
// ============================================================

void checkMoab(const moab::ErrorCode error, const char* message) {
  if (error != moab::MB_SUCCESS) {
    throw std::runtime_error(message);
  }
}

double norm(const Vec3& v) noexcept {
  return std::sqrt(mesh::dot(v, v));
}

Vec3 normalize(const Vec3& v) {
  const double length = norm(v);

  if (length <= std::numeric_limits<double>::epsilon()) {
    throw std::runtime_error("cannot normalize zero-length vector");
  }

  return v / length;
}

// ============================================================
// Polygon geometry
//
// Current pemu mesh backend assumes a planar XY 2D mesh.
//
// For a polygon with ordered vertices:
//
// A = 1/2 sum_i (x_i y_i+1 - x_i+1 y_i)
//
// Cx = 1/(6A) sum_i
//      (x_i + x_i+1)
//      (x_i y_i+1 - x_i+1 y_i)
//
// Cy analogous.
//
// This works for triangle, quad, and general simple polygon.
// ============================================================

struct PolygonGeometry {
  Vec3 centroid{};
  double area{};
};

PolygonGeometry computePolygonGeometry(const std::vector<Vec3>& vertices) {
  if (vertices.size() < 3) {
    throw std::runtime_error("cell has fewer than three vertices");
  }

  double twice_signed_area = 0.0;

  double cx_numerator = 0.0;
  double cy_numerator = 0.0;

  double z_sum = 0.0;

  const std::size_t n = vertices.size();

  for (std::size_t i = 0; i < n; ++i) {

    const Vec3& p = vertices[i];

    const Vec3& q = vertices[(i + 1) % n];

    const double cross = p.x * q.y - q.x * p.y;

    twice_signed_area += cross;

    cx_numerator += (p.x + q.x) * cross;

    cy_numerator += (p.y + q.y) * cross;

    z_sum += p.z;
  }

  const double signed_area = 0.5 * twice_signed_area;

  if (std::abs(signed_area) <= std::numeric_limits<double>::epsilon()) {
    throw std::runtime_error("degenerate zero-area cell");
  }

  const double cx = cx_numerator / (6.0 * signed_area);

  const double cy = cy_numerator / (6.0 * signed_area);

  //
  // For a planar XY mesh this is normally zero.
  // Keep the average z coordinate so Vec3 remains useful.
  //
  const double cz = z_sum / static_cast<double>(n);

  return {.centroid = {cx, cy, cz}, .area = std::abs(signed_area)};
}

std::uint32_t checkedUint32(const std::size_t value) {
  if (value > std::numeric_limits<std::uint32_t>::max()) {
    throw std::overflow_error("mesh connectivity exceeds uint32_t range");
  }

  return static_cast<std::uint32_t>(value);
}

}  // namespace

// ============================================================
// Construction
// ============================================================

MoabMesh::MoabMesh(const std::string& filename)
    : core_(std::make_unique<moab::Core>()) {
  load(filename);

  buildEntities();
  buildHandleMaps();

  buildTopology();
  buildVertexTopology();
  buildGeometry();

  buildBoundaryMetadata();
}

// ============================================================
// Load
// ============================================================

void MoabMesh::load(const std::string& filename) {
  if (filename.empty()) {
    throw std::invalid_argument("mesh filename must not be empty");
  }

  const auto error = core_->load_mesh(filename.c_str());

  checkMoab(error, "MOAB failed to load mesh");
}

// ============================================================
// Entity discovery
//
// Important:
//
// Some mesh files contain only cells + vertices and do not
// explicitly store all edges/faces.
//
// Therefore:
//
//     cells
//       ↓
// get_adjacencies(...,
//                 dimension - 1,
//                 create_if_missing = true)
//       ↓
//     faces
//
// MOAB creates missing lower-dimensional entities.
// ============================================================

void MoabMesh::buildEntities() {
  moab::Range cells_3d;
  moab::Range cells_2d;

  //
  // Determine topological cell dimension.
  //
  checkMoab(core_->get_entities_by_dimension(0, 3, cells_3d),
            "MOAB failed to query 3D cells");

  if (!cells_3d.empty()) {

    //
    // Current implementation intentionally stops here.
    //
    // Geometry implementation below is 2D FVM geometry.
    //
    throw std::runtime_error("MoabMesh currently supports only 2D meshes");

  } else {

    checkMoab(core_->get_entities_by_dimension(0, 2, cells_2d),
              "MOAB failed to query 2D cells");

    if (cells_2d.empty()) {
      throw std::runtime_error("mesh contains no supported 2D cells");
    }

    dimension_ = 2;
  }

  //
  // Copy cells.
  //
  cell_handles_.assign(cells_2d.begin(), cells_2d.end());

  //
  // Ask MOAB to create/retrieve all cell edges.
  //
  // For a 2D mesh:
  //
  //     cell dimension = 2
  //     face dimension = 1
  //
  moab::Range faces;

  const auto face_error = core_->get_adjacencies(cells_2d, dimension_ - 1, true,
                                                 faces, moab::Interface::UNION);

  checkMoab(face_error, "MOAB failed to construct/query cell faces");

  if (faces.empty()) {
    throw std::runtime_error("mesh contains no faces");
  }

  face_handles_.assign(faces.begin(), faces.end());

  //
  // Obtain vertices from cells.
  //
  moab::Range vertices;

  const auto vertex_error = core_->get_connectivity(cells_2d, vertices);

  checkMoab(vertex_error, "MOAB failed to query mesh vertices");

  if (vertices.empty()) {
    throw std::runtime_error("mesh contains no vertices");
  }

  vertex_handles_.assign(vertices.begin(), vertices.end());
}

// ============================================================
// Dense local ID maps
//
// MOAB EntityHandle is deliberately kept inside this backend.
//
// pemu runtime:
//
//     CellId = [0, numCells)
//     FaceId = [0, numFaces)
//
// This lets fields remain dense:
//
//     phi[cell]
//     rho[cell]
//     flux[face]
// ============================================================

void MoabMesh::buildHandleMaps() {
  cell_id_.clear();
  face_id_.clear();
  vertex_id_.clear();

  cell_id_.reserve(cell_handles_.size());

  face_id_.reserve(face_handles_.size());
  vertex_id_.reserve(vertex_handles_.size());

  for (std::size_t i = 0; i < cell_handles_.size(); ++i) {

    if (i > std::numeric_limits<CellId>::max()) {
      throw std::overflow_error("number of cells exceeds CellId range");
    }

    cell_id_.emplace(cell_handles_[i], static_cast<CellId>(i));
  }

  for (std::size_t i = 0; i < face_handles_.size(); ++i) {

    if (i > std::numeric_limits<FaceId>::max()) {
      throw std::overflow_error("number of faces exceeds FaceId range");
    }

    face_id_.emplace(face_handles_[i], static_cast<FaceId>(i));
  }

  for (std::size_t i = 0; i < vertex_handles_.size(); ++i) {
    if (i > std::numeric_limits<VertexId>::max()) {
      throw std::overflow_error("number of vertices exceeds VertexId range");
    }

    vertex_id_.emplace(vertex_handles_[i], static_cast<VertexId>(i));
  }
}

// ============================================================
// Topology
//
// Build:
//
//     face -> owner
//     face -> neighbor
//
// and:
//
//     cell -> faces
//
// cell -> faces is stored in CSR-like form:
//
//     cell_face_offsets_[cell]
//     cell_faces_[offset ... next_offset)
//
// A manifold interior face must touch exactly two cells.
// A boundary face must touch exactly one cell.
// ============================================================

void MoabMesh::buildTopology() {
  const std::size_t num_faces = face_handles_.size();

  const std::size_t num_cells = cell_handles_.size();

  face_owner_.assign(num_faces, invalid_cell);

  face_neighbor_.assign(num_faces, invalid_cell);

  std::vector<std::vector<FaceId>> temporary_cell_faces(num_cells);

  for (std::size_t f = 0; f < num_faces; ++f) {

    const auto face_handle = face_handles_[f];

    moab::Range adjacent_cells;

    const auto error =
        core_->get_adjacencies(&face_handle, 1, dimension_, false,
                               adjacent_cells, moab::Interface::UNION);

    checkMoab(error, "MOAB failed to query face-cell adjacency");

    //
    // FVM manifold assumption:
    //
    // boundary face -> 1 adjacent cell
    // internal face -> 2 adjacent cells
    //
    if (adjacent_cells.empty()) {
      throw std::runtime_error("orphan face has no adjacent cell");
    }

    if (adjacent_cells.size() > 2) {
      throw std::runtime_error(
          "non-manifold face has more than two adjacent cells");
    }

    const FaceId face_id = static_cast<FaceId>(f);

    auto it = adjacent_cells.begin();

    //
    // Owner selection is arbitrary at this point.
    //
    // Geometry construction later guarantees normal:
    //
    //     owner -> neighbor
    //
    const auto owner_handle = *it;

    const auto owner_map_it = cell_id_.find(owner_handle);

    if (owner_map_it == cell_id_.end()) {
      throw std::runtime_error("adjacent owner cell missing from dense map");
    }

    const CellId owner = owner_map_it->second;

    face_owner_[face_id] = owner;

    temporary_cell_faces[owner].push_back(face_id);

    ++it;

    if (it != adjacent_cells.end()) {

      const auto neighbor_handle = *it;

      const auto neighbor_map_it = cell_id_.find(neighbor_handle);

      if (neighbor_map_it == cell_id_.end()) {
        throw std::runtime_error(
            "adjacent neighbor cell missing from dense map");
      }

      const CellId neighbor = neighbor_map_it->second;

      face_neighbor_[face_id] = neighbor;

      temporary_cell_faces[neighbor].push_back(face_id);
    }
  }

  //
  // Flatten connectivity.
  //
  cell_face_offsets_.resize(num_cells + 1);

  std::size_t total_faces = 0;

  for (std::size_t cell = 0; cell < num_cells; ++cell) {

    cell_face_offsets_[cell] = checkedUint32(total_faces);

    total_faces += temporary_cell_faces[cell].size();
  }

  cell_face_offsets_[num_cells] = checkedUint32(total_faces);

  cell_faces_.clear();
  cell_faces_.reserve(total_faces);

  for (const auto& faces : temporary_cell_faces) {

    cell_faces_.insert(cell_faces_.end(), faces.begin(), faces.end());
  }
}

// ============================================================
// Entity vertex coordinates
//
// This helper intentionally remains local to this translation
// unit because MOAB connectivity must never escape the backend.
// ============================================================

static std::vector<Vec3> getEntityVertices(const moab::Core& core,
                                           const moab::EntityHandle entity) {
  const moab::EntityHandle* connectivity = nullptr;

  int num_vertices = 0;

  const auto connectivity_error =
      core.get_connectivity(entity, connectivity, num_vertices, true);

  checkMoab(connectivity_error, "MOAB failed to query entity connectivity");

  if (connectivity == nullptr || num_vertices <= 0) {
    throw std::runtime_error("entity has invalid connectivity");
  }

  std::vector<double> coordinates(static_cast<std::size_t>(3 * num_vertices));

  const auto coordinate_error =
      core.get_coords(connectivity, num_vertices, coordinates.data());

  checkMoab(coordinate_error, "MOAB failed to query vertex coordinates");

  std::vector<Vec3> result(static_cast<std::size_t>(num_vertices));

  for (int i = 0; i < num_vertices; ++i) {

    result[static_cast<std::size_t>(i)] = {
        coordinates[static_cast<std::size_t>(3 * i + 0)],

        coordinates[static_cast<std::size_t>(3 * i + 1)],

        coordinates[static_cast<std::size_t>(3 * i + 2)]};
  }

  return result;
}

// ============================================================
// Ordered cell -> vertex connectivity and dense vertex geometry
// ============================================================

void MoabMesh::buildVertexTopology() {
  const std::size_t num_cells = cell_handles_.size();

  cell_vertex_offsets_.resize(num_cells + 1);
  cell_vertices_.clear();

  for (std::size_t c = 0; c < num_cells; ++c) {
    cell_vertex_offsets_[c] = checkedUint32(cell_vertices_.size());

    const moab::EntityHandle* connectivity = nullptr;
    int num_vertices = 0;
    const auto error = core_->get_connectivity(
        cell_handles_[c], connectivity, num_vertices, true);

    checkMoab(error, "MOAB failed to query cell vertex connectivity");

    if (connectivity == nullptr || num_vertices < 3) {
      throw std::runtime_error("cell has invalid vertex connectivity");
    }

    for (int i = 0; i < num_vertices; ++i) {
      const auto vertex_it = vertex_id_.find(connectivity[i]);
      if (vertex_it == vertex_id_.end()) {
        throw std::runtime_error("cell vertex missing from dense map");
      }

      cell_vertices_.push_back(vertex_it->second);
    }
  }

  cell_vertex_offsets_[num_cells] = checkedUint32(cell_vertices_.size());

  vertex_coordinates_.resize(vertex_handles_.size());
  if (!vertex_handles_.empty()) {
    std::vector<double> coordinates(3 * vertex_handles_.size());
    checkMoab(core_->get_coords(vertex_handles_.data(), vertex_handles_.size(),
                                coordinates.data()),
              "MOAB failed to query dense vertex coordinates");

    for (std::size_t i = 0; i < vertex_handles_.size(); ++i) {
      vertex_coordinates_[i] = {coordinates[3 * i], coordinates[3 * i + 1],
                                coordinates[3 * i + 2]};
    }
  }
}

// ============================================================
// Geometry
//
// Current contract:
//
// 2D planar XY finite-volume mesh.
//
// cellVolume(cell):
//     geometric cell area
//
// faceArea(face):
//     edge length
//
// faceNormal(face):
//     unit normal
//
// Normal orientation invariant:
//
// internal face:
//     owner -> neighbor
//
// boundary face:
//     owner -> outside
// ============================================================

void MoabMesh::buildGeometry() {
  const std::size_t num_cells = cell_handles_.size();

  const std::size_t num_faces = face_handles_.size();

  cell_centers_.resize(num_cells);
  cell_volumes_.resize(num_cells);

  face_centers_.resize(num_faces);
  face_normals_.resize(num_faces);
  face_areas_.resize(num_faces);

  // --------------------------------------------------------
  // Cells
  // --------------------------------------------------------

  for (std::size_t c = 0; c < num_cells; ++c) {

    const auto vertices = getEntityVertices(*core_, cell_handles_[c]);

    const auto geometry = computePolygonGeometry(vertices);

    cell_centers_[c] = geometry.centroid;

    cell_volumes_[c] = geometry.area;
  }

  // --------------------------------------------------------
  // Faces
  //
  // In 2D each face is an edge with two vertices.
  // --------------------------------------------------------

  for (std::size_t f = 0; f < num_faces; ++f) {

    const auto vertices = getEntityVertices(*core_, face_handles_[f]);

    if (vertices.size() != 2) {
      throw std::runtime_error("2D face must contain exactly two vertices");
    }

    const Vec3& a = vertices[0];

    const Vec3& b = vertices[1];

    const Vec3 edge = b - a;

    const double length = norm(edge);

    if (length <= std::numeric_limits<double>::epsilon()) {
      throw std::runtime_error("zero-length mesh face");
    }

    const FaceId face_id = static_cast<FaceId>(f);

    //
    // Midpoint.
    //
    const Vec3 center = (a + b) / 2.0;

    face_centers_[face_id] = center;

    face_areas_[face_id] = length;

    //
    // Candidate normal for XY plane.
    //
    // Edge:
    //
    //     (dx, dy)
    //
    // perpendicular:
    //
    //     (dy, -dx)
    //
    Vec3 normal{edge.y, -edge.x, 0.0};

    normal = normalize(normal);

    //
    // Enforce pemu orientation convention.
    //
    const CellId owner = face_owner_[face_id];

    const CellId neighbor = face_neighbor_[face_id];

    Vec3 desired_direction{};

    if (neighbor != invalid_cell) {

      //
      // internal:
      //
      // owner -> neighbor
      //
      desired_direction = cell_centers_[neighbor] - cell_centers_[owner];

    } else {

      //
      // boundary:
      //
      // owner center -> boundary face center
      //
      // This points outward for a valid convex local cell.
      //
      desired_direction = center - cell_centers_[owner];
    }

    if (dot(normal, desired_direction) < 0.0) {

      normal = normal * -1.0;
    }

    face_normals_[face_id] = normal;
  }
}

// ============================================================
// Boundary metadata
//
// v0.1:
//
// all boundary faces are assigned BoundaryId = 0.
//
// Interior faces:
//
//     invalid_boundary
//
// Later this function can map MOAB entity sets / tags:
//
//     wall
//     electrode
//     symmetry
//     dielectric
//
// into pemu BoundaryId values.
// ============================================================

void MoabMesh::buildBoundaryMetadata() {
  //
  // Default:
  //
  // internal face -> invalid_boundary
  //
  // boundary face with no physical tag
  //               -> invalid_boundary
  //
  boundary_ids_.assign(numFaces(), invalid_boundary);

  // --------------------------------------------------------
  // Gmsh physical groups are imported by MOAB as
  // MATERIAL_SET entity sets.
  //
  // Example:
  //
  // physical id = 1 -> left
  // physical id = 2 -> right
  //
  // The numeric id is stored in MATERIAL_SET.
  // --------------------------------------------------------

  moab::Tag material_tag{};

  const auto tag_error = core_->tag_get_handle(
      MATERIAL_SET_TAG_NAME, 1, moab::MB_TYPE_INTEGER, material_tag);

  //
  // No MATERIAL_SET tag simply means:
  //
  // mesh contains no physical groups.
  //
  if (tag_error == moab::MB_TAG_NOT_FOUND) {
    return;
  }

  checkMoab(tag_error, "MOAB failed to query MATERIAL_SET tag");

  // --------------------------------------------------------
  // Find all entity sets carrying MATERIAL_SET.
  // --------------------------------------------------------

  moab::Range material_sets;

  const auto set_error = core_->get_entities_by_type_and_tag(
      0, moab::MBENTITYSET, &material_tag, nullptr, 1, material_sets);

  checkMoab(set_error, "MOAB failed to query MATERIAL_SET entity sets");

  // --------------------------------------------------------
  // Iterate physical groups.
  // --------------------------------------------------------

  for (const auto set : material_sets) {

    int physical_id = 0;

    const auto id_error =
        core_->tag_get_data(material_tag, &set, 1, &physical_id);

    checkMoab(id_error, "MOAB failed to read MATERIAL_SET id");

    //
    // pemu BoundaryId is unsigned.
    //
    if (physical_id < 0) {
      throw std::runtime_error("negative physical group id is not supported");
    }

    // ----------------------------------------------------
    // Query entities contained in the physical group.
    //
    // Recursive = true is deliberate:
    //
    // a set may contain other sets instead of directly
    // containing edges.
    // ----------------------------------------------------

    moab::Range entities;

    const auto entity_error =
        core_->get_entities_by_handle(set, entities, true);

    checkMoab(entity_error, "MOAB failed to query MATERIAL_SET contents");

    // ----------------------------------------------------
    // Only dimension-1 entities can represent FVM
    // boundary faces in a 2D mesh.
    //
    // domain physical group (dimension 2) is ignored here.
    // ----------------------------------------------------

    for (const auto entity : entities) {

      const auto type = core_->type_from_handle(entity);

      const int entity_dimension = moab::CN::Dimension(type);

      if (entity_dimension != dimension_ - 1) {

        //
        // e.g. physical group "domain"
        //
        continue;
      }

      // ------------------------------------------------
      // Map MOAB edge entity -> pemu FaceId.
      // ------------------------------------------------

      const auto face_it = face_id_.find(entity);

      if (face_it == face_id_.end()) {

        //
        // This edge may not belong to the active cell
        // topology represented by pemu.
        //
        continue;
      }

      const FaceId face = face_it->second;

      // ------------------------------------------------
      // Physical groups of dimension-1 are meaningful as
      // boundary groups only if the corresponding face
      // is actually topologically on the boundary.
      //
      // Do not silently accept an internal interface as
      // external boundary.
      // ------------------------------------------------

      if (!isBoundary(face)) {
        continue;
      }

      const BoundaryId id = static_cast<BoundaryId>(physical_id);

      // ------------------------------------------------
      // A face should normally belong to exactly one
      // physical boundary group.
      //
      // Catch conflicting assignment early.
      // ------------------------------------------------

      if (boundary_ids_[face] != invalid_boundary &&
          boundary_ids_[face] != id) {

        throw std::runtime_error(
            "boundary face belongs to multiple "
            "physical groups");
      }

      boundary_ids_[face] = id;
    }
  }
}

// ============================================================
// Basic information
// ============================================================

int MoabMesh::dimension() const noexcept {
  return dimension_;
}

std::size_t MoabMesh::numCells() const noexcept {
  return cell_handles_.size();
}

std::size_t MoabMesh::numFaces() const noexcept {
  return face_handles_.size();
}

std::size_t MoabMesh::numVertices() const noexcept {
  return vertex_handles_.size();
}

// ============================================================
// Topology query
// ============================================================

CellId MoabMesh::owner(const FaceId face) const {
  if (face >= face_owner_.size()) {
    throw std::out_of_range("FaceId out of range");
  }

  return face_owner_[face];
}

CellId MoabMesh::neighbor(const FaceId face) const {
  if (face >= face_neighbor_.size()) {
    throw std::out_of_range("FaceId out of range");
  }

  return face_neighbor_[face];
}

bool MoabMesh::isBoundary(const FaceId face) const {
  if (face >= face_neighbor_.size()) {
    throw std::out_of_range("FaceId out of range");
  }

  return face_neighbor_[face] == invalid_cell;
}

std::span<const FaceId> MoabMesh::cellFaces(const CellId cell) const {
  if (cell >= cell_handles_.size()) {
    throw std::out_of_range("CellId out of range");
  }

  const std::uint32_t begin = cell_face_offsets_[cell];

  const std::uint32_t end = cell_face_offsets_[cell + 1];

  return std::span<const FaceId>{cell_faces_.data() + begin,
                                 static_cast<std::size_t>(end - begin)};
}

std::span<const VertexId> MoabMesh::cellVertices(const CellId cell) const {
  if (cell >= cell_handles_.size()) {
    throw std::out_of_range("CellId out of range");
  }

  const std::uint32_t begin = cell_vertex_offsets_[cell];
  const std::uint32_t end = cell_vertex_offsets_[cell + 1];

  return std::span<const VertexId>{cell_vertices_.data() + begin,
                                   static_cast<std::size_t>(end - begin)};
}

// ============================================================
// Geometry query
// ============================================================

Vec3 MoabMesh::cellCenter(const CellId cell) const {
  if (cell >= cell_centers_.size()) {
    throw std::out_of_range("CellId out of range");
  }

  return cell_centers_[cell];
}

Vec3 MoabMesh::vertex(const VertexId vertex) const {
  if (vertex >= vertex_coordinates_.size()) {
    throw std::out_of_range("VertexId out of range");
  }

  return vertex_coordinates_[vertex];
}

double MoabMesh::cellVolume(const CellId cell) const {
  if (cell >= cell_volumes_.size()) {
    throw std::out_of_range("CellId out of range");
  }

  return cell_volumes_[cell];
}

Vec3 MoabMesh::faceCenter(const FaceId face) const {
  if (face >= face_centers_.size()) {
    throw std::out_of_range("FaceId out of range");
  }

  return face_centers_[face];
}

double MoabMesh::faceArea(const FaceId face) const {
  if (face >= face_areas_.size()) {
    throw std::out_of_range("FaceId out of range");
  }

  return face_areas_[face];
}

Vec3 MoabMesh::faceNormal(const FaceId face) const {
  if (face >= face_normals_.size()) {
    throw std::out_of_range("FaceId out of range");
  }

  return face_normals_[face];
}

// ============================================================
// Boundary
// ============================================================

BoundaryId MoabMesh::boundaryId(const FaceId face) const {
  if (face >= boundary_ids_.size()) {
    throw std::out_of_range("FaceId out of range");
  }

  return boundary_ids_[face];
}

}  // namespace pemu::mesh
