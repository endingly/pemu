#include <pemu/output/common/mesh_adapter.hpp>

#include <vtkCellType.h>
#include <vtkIdList.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>

#include <limits>
#include <stdexcept>

namespace pemu::output::common {

vtkSmartPointer<vtkUnstructuredGrid> toVtkUnstructuredGrid(
    const mesh::IMesh& mesh) {
  if (mesh.numVertices() >
          static_cast<std::size_t>(std::numeric_limits<vtkIdType>::max()) ||
      mesh.numCells() >
          static_cast<std::size_t>(std::numeric_limits<vtkIdType>::max()) ||
      mesh.numFaces() >
          static_cast<std::size_t>(std::numeric_limits<vtkIdType>::max())) {
    throw std::overflow_error("mesh exceeds VTK id range");
  }

  vtkNew<vtkPoints> points;
  points->SetDataTypeToDouble();
  points->SetNumberOfPoints(static_cast<vtkIdType>(mesh.numVertices()));

  for (mesh::VertexId vertex = 0; vertex < mesh.numVertices(); ++vertex) {
    const auto point = mesh.vertex(vertex);
    points->SetPoint(static_cast<vtkIdType>(vertex), point.x, point.y, point.z);
  }

  vtkNew<vtkUnstructuredGrid> grid;
  grid->SetPoints(points);
  grid->Allocate(static_cast<vtkIdType>(mesh.numCells()));

  for (mesh::CellId cell = 0; cell < mesh.numCells(); ++cell) {
    const auto vertices = mesh.cellVertices(cell);
    if (vertices.size() < 3) {
      throw std::invalid_argument("VTK polygon cell needs at least 3 vertices");
    }

    vtkNew<vtkIdList> ids;
    ids->SetNumberOfIds(static_cast<vtkIdType>(vertices.size()));
    for (std::size_t i = 0; i < vertices.size(); ++i) {
      ids->SetId(static_cast<vtkIdType>(i),
                 static_cast<vtkIdType>(vertices[i]));
    }

    grid->InsertNextCell(VTK_POLYGON, ids);
  }

  return grid;
}

}  // namespace pemu::output::common
