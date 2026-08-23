#pragma once

#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>
#include <pemu/field/utils.hpp>
#include <pemu/mesh/geometry.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace pemu::discretization::operators {
namespace detail {

/** @brief Returns one Cartesian component by runtime index. */
[[nodiscard]] inline double vectorComponent(const mesh::Vec3& value,
                                            int index) noexcept {
  return index == 0 ? value.x : (index == 1 ? value.y : value.z);
}

/** @brief Solves a normal-equation system of dimension one to three. */
[[nodiscard]] inline std::array<double, 3> solveVectorNormalEquations(
    std::array<std::array<double, 4>, 3> augmented, int dimension) {
  double coefficient_scale = 0.0;
  for (int row = 0; row < dimension; ++row) {
    for (int column = 0; column < dimension; ++column) {
      coefficient_scale =
          std::max(coefficient_scale, std::abs(augmented[row][column]));
    }
  }
  if (!(coefficient_scale > 0.0) || !std::isfinite(coefficient_scale)) {
    throw std::runtime_error(
        "cell face normals cannot reconstruct a cell vector");
  }

  const double pivot_tolerance =
      64.0 * std::numeric_limits<double>::epsilon() * coefficient_scale;
  for (int pivot = 0; pivot < dimension; ++pivot) {
    int best_row = pivot;
    for (int row = pivot + 1; row < dimension; ++row) {
      if (std::abs(augmented[row][pivot]) >
          std::abs(augmented[best_row][pivot])) {
        best_row = row;
      }
    }
    if (std::abs(augmented[best_row][pivot]) <= pivot_tolerance) {
      throw std::runtime_error(
          "cell face normals are rank deficient for vector reconstruction");
    }
    std::swap(augmented[pivot], augmented[best_row]);

    for (int row = pivot + 1; row < dimension; ++row) {
      const double factor = augmented[row][pivot] / augmented[pivot][pivot];
      for (int column = pivot; column <= dimension; ++column) {
        augmented[row][column] -= factor * augmented[pivot][column];
      }
    }
  }

  std::array<double, 3> result{};
  for (int row = dimension - 1; row >= 0; --row) {
    double value = augmented[row][dimension];
    for (int column = row + 1; column < dimension; ++column) {
      value -= augmented[row][column] * result[column];
    }
    result[row] = value / augmented[row][row];
    if (!std::isfinite(result[row])) {
      throw std::runtime_error("cell-vector reconstruction overflowed");
    }
  }
  return result;
}

}  // namespace detail

/**
 * @brief Reconstructs a cell-vector magnitude from face-normal components.
 *
 * For each cell this solves the face-area-weighted least-squares problem
 * min_v sum_f A_f (n_f dot v - v_n,f)^2. It is exact for a constant Cartesian
 * vector when the cell normals span the mesh dimension.
 */
inline void reconstructCellVectorMagnitudeFromFaceNormal(
    const field::FaceField<double>& face_normal_component,
    field::CellField<double>& cell_magnitude) {
  field::ensureSameMesh(face_normal_component, cell_magnitude);
  const auto& mesh = face_normal_component.mesh();
  const int dimension = mesh.dimension();
  if (dimension < 1 || dimension > 3) {
    throw std::invalid_argument(
        "cell-vector reconstruction requires mesh dimension 1, 2, or 3");
  }
  for (const double value : face_normal_component) {
    if (!std::isfinite(value)) {
      throw std::invalid_argument(
          "face-normal component must contain only finite values");
    }
  }

  for (mesh::CellId cell = 0; cell < mesh.numCells(); ++cell) {
    std::array<std::array<double, 4>, 3> normal_equations{};
    for (const mesh::FaceId face : mesh.cellFaces(cell)) {
      const double area = mesh.faceArea(face);
      const auto normal = mesh.faceNormal(face);
      if (!std::isfinite(area) || area <= 0.0 || !std::isfinite(normal.x) ||
          !std::isfinite(normal.y) || !std::isfinite(normal.z)) {
        throw std::runtime_error(
            "cell-vector reconstruction requires finite face geometry");
      }

      for (int row = 0; row < dimension; ++row) {
        const double normal_row = detail::vectorComponent(normal, row);
        for (int column = 0; column < dimension; ++column) {
          normal_equations[row][column] +=
              area * normal_row * detail::vectorComponent(normal, column);
        }
        normal_equations[row][dimension] +=
            area * normal_row * face_normal_component[face];
      }
    }

    const auto reconstructed =
        detail::solveVectorNormalEquations(normal_equations, dimension);
    double magnitude_squared = 0.0;
    for (int component = 0; component < dimension; ++component) {
      magnitude_squared += reconstructed[component] * reconstructed[component];
    }
    cell_magnitude[cell] = std::sqrt(magnitude_squared);
    if (!std::isfinite(cell_magnitude[cell])) {
      throw std::runtime_error("cell-vector magnitude overflowed");
    }
  }
}

}  // namespace pemu::discretization::operators
