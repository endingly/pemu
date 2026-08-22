#pragma once

#include <pemu/mesh/i_mesh.hpp>

#include <vtkSmartPointer.h>

class vtkUnstructuredGrid;

namespace pemu::output::common {

[[nodiscard]] vtkSmartPointer<vtkUnstructuredGrid> toVtkUnstructuredGrid(
    const mesh::IMesh& mesh);

}  // namespace pemu::output::common
