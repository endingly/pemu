#pragma once

#include <pemu/field/types.hpp>
#include <pemu/field/utils.hpp>
#include <pemu/mesh/i_mesh.hpp>

#include <concepts>
#include <type_traits>

namespace pemu::discretization::operators {

template <field::FaceFieldLike FluxField, field::CellFieldLike ResultField>
  requires std::same_as<typename FluxField::value_type,
                        typename ResultField::value_type> &&
           std::floating_point<typename FluxField::value_type>
void divergence(const FluxField& normal_flux, ResultField& result) {
  field::ensureSameMesh(normal_flux, result);

  const auto& mesh = normal_flux.mesh();

  result.fill(typename ResultField::value_type{0});

  // --------------------------------------------------------
  // FaceField stores:
  //
  //     Gamma_f · n_f
  //
  // where n_f is:
  //
  // internal:
  //     owner -> neighbor
  //
  // boundary:
  //     owner -> outside
  //
  // Therefore:
  //
  // owner:
  //     + Gamma_n A_f
  //
  // neighbor:
  //     - Gamma_n A_f
  // --------------------------------------------------------

  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {

    const auto owner = mesh.owner(face);

    const auto integrated_flux = normal_flux[face] * mesh.faceArea(face);

    result[owner] += integrated_flux;

    if (!mesh.isBoundary(face)) {

      const auto neighbor = mesh.neighbor(face);

      result[neighbor] -= integrated_flux;
    }
  }

  // --------------------------------------------------------
  // Convert integrated flux to divergence:
  //
  //        1
  // div = --- sum_f Gamma_f · n_f A_f
  //        V
  // --------------------------------------------------------

  for (mesh::CellId cell = 0; cell < mesh.numCells(); ++cell) {

    result[cell] /= mesh.cellVolume(cell);
  }
}

}  // namespace pemu::discretization::operators