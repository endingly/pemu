#pragma once

#include <pemu/field/cell_field.hpp>
#include <pemu/field/face_field.hpp>
#include <pemu/field/utils.hpp>
#include <pemu/mesh/geometry.hpp>
#include <pemu/mesh/i_mesh.hpp>

#include <cmath>
#include <stdexcept>

namespace pemu::physics {

/**
 * @brief Reconstructs cell-centred electron field power from face-normal data.
 *
 * The particle flux and electric field use the mesh face-normal orientation.
 * Their product is orientation independent. For an orthogonal finite-volume
 * cell, each face contribution is weighted by its centre-to-face distance:
 *
 *   S_E,P = -1/V_P sum_f Gamma_e,f E_f A_f d_P,f.
 *
 * With particle flux expressed in 1/(area*time) and electric field in V/length,
 * the returned numerical value is in eV/(volume*time). No explicit elementary
 * charge factor is needed because the target energy unit is electronvolt.
 */
inline void computeElectronFieldPowerDensity(
    const field::FaceField<double>& electron_particle_flux_normal,
    const field::FaceField<double>& electric_field_normal,
    field::CellField<double>& field_power_density) {
  field::ensureSameMesh(electron_particle_flux_normal, electric_field_normal);
  field::ensureSameMesh(electron_particle_flux_normal, field_power_density);

  const auto& mesh = electron_particle_flux_normal.mesh();
  field_power_density.fill(0.0);

  for (mesh::FaceId face = 0; face < mesh.numFaces(); ++face) {
    const auto owner = mesh.owner(face);
    const auto normal = mesh.faceNormal(face);
    const double owner_distance =
        mesh::dot(mesh.faceCenter(face) - mesh.cellCenter(owner), normal);
    if (!std::isfinite(owner_distance) || owner_distance <= 0.0) {
      throw std::runtime_error(
          "electron field power requires positive owner-face distance");
    }

    const double face_area = mesh.faceArea(face);
    const double owner_volume = mesh.cellVolume(owner);
    if (!std::isfinite(face_area) || face_area <= 0.0 ||
        !std::isfinite(owner_volume) || owner_volume <= 0.0) {
      throw std::runtime_error(
          "electron field power requires positive finite geometry");
    }
    if (!std::isfinite(electron_particle_flux_normal[face]) ||
        !std::isfinite(electric_field_normal[face])) {
      throw std::invalid_argument(
          "electron field power requires finite face fields");
    }

    const double face_power =
        -electron_particle_flux_normal[face] * electric_field_normal[face] *
        face_area;
    field_power_density[owner] +=
        face_power * owner_distance / owner_volume;

    if (mesh.isBoundary(face)) {
      continue;
    }

    const auto neighbor = mesh.neighbor(face);
    const double neighbor_distance =
        mesh::dot(mesh.cellCenter(neighbor) - mesh.faceCenter(face), normal);
    if (!std::isfinite(neighbor_distance) || neighbor_distance <= 0.0) {
      throw std::runtime_error(
          "electron field power requires positive neighbor-face distance");
    }
    const double neighbor_volume = mesh.cellVolume(neighbor);
    if (!std::isfinite(neighbor_volume) || neighbor_volume <= 0.0) {
      throw std::runtime_error(
          "electron field power requires positive finite geometry");
    }
    field_power_density[neighbor] +=
        face_power * neighbor_distance / neighbor_volume;
  }
}

/**
 * @brief Computes electron energy density from number density and mean energy.
 *
 * All inputs are numerical values in one consistent unit system, for example
 * 1/cm^3 and eV, producing eV/cm^3.
 */
inline void computeElectronEnergyDensity(
    const field::CellField<double>& number_density,
    const field::CellField<double>& mean_energy,
    field::CellField<double>& energy_density) {
  field::ensureSameMesh(number_density, mean_energy);
  field::ensureSameMesh(number_density, energy_density);

  for (mesh::CellId cell = 0; cell < number_density.mesh().numCells(); ++cell) {
    if (!std::isfinite(number_density[cell]) || number_density[cell] < 0.0) {
      throw std::invalid_argument(
          "electron number density must be finite and non-negative");
    }
    if (!std::isfinite(mean_energy[cell]) || mean_energy[cell] < 0.0) {
      throw std::invalid_argument(
          "electron mean energy must be finite and non-negative");
    }
    if (!std::isfinite(number_density[cell] * mean_energy[cell])) {
      throw std::overflow_error("electron energy density overflowed");
    }
  }

  for (mesh::CellId cell = 0; cell < number_density.mesh().numCells(); ++cell) {
    energy_density[cell] = number_density[cell] * mean_energy[cell];
  }
}

/**
 * @brief Recovers mean electron energy while defining vacuum cells as zero.
 *
 * Cells whose number density does not exceed density_floor receive zero mean
 * energy, avoiding division by a vanishing density.
 */
inline void computeElectronMeanEnergy(
    const field::CellField<double>& energy_density,
    const field::CellField<double>& number_density, double density_floor,
    field::CellField<double>& mean_energy) {
  field::ensureSameMesh(energy_density, number_density);
  field::ensureSameMesh(energy_density, mean_energy);

  if (!std::isfinite(density_floor) || density_floor < 0.0) {
    throw std::invalid_argument(
        "electron density floor must be finite and non-negative");
  }

  for (mesh::CellId cell = 0; cell < energy_density.mesh().numCells(); ++cell) {
    if (!std::isfinite(energy_density[cell]) || energy_density[cell] < 0.0) {
      throw std::invalid_argument(
          "electron energy density must be finite and non-negative");
    }
    if (!std::isfinite(number_density[cell]) || number_density[cell] < 0.0) {
      throw std::invalid_argument(
          "electron number density must be finite and non-negative");
    }
    if (number_density[cell] > density_floor &&
        !std::isfinite(energy_density[cell] / number_density[cell])) {
      throw std::overflow_error("electron mean energy overflowed");
    }
  }

  for (mesh::CellId cell = 0; cell < energy_density.mesh().numCells(); ++cell) {
    mean_energy[cell] = number_density[cell] > density_floor
                            ? energy_density[cell] / number_density[cell]
                            : 0.0;
  }
}

/**
 * @brief Converts mean electron energy in eV to equivalent Maxwellian Te in eV.
 *
 * Plasma notation reports k_B*T_e as an energy. For an isotropic Maxwellian
 * distribution, mean_energy = 3/2*k_B*T_e.
 */
[[nodiscard]] inline double electronTemperatureEv(double mean_energy_ev) {
  if (!std::isfinite(mean_energy_ev) || mean_energy_ev < 0.0) {
    throw std::invalid_argument(
        "electron mean energy must be finite and non-negative");
  }
  return (2.0 / 3.0) * mean_energy_ev;
}

/** @brief Converts a mean-energy field in eV to equivalent Te values in eV. */
inline void computeElectronTemperatureEv(
    const field::CellField<double>& mean_energy_ev,
    field::CellField<double>& electron_temperature_ev) {
  field::ensureSameMesh(mean_energy_ev, electron_temperature_ev);
  for (const double value : mean_energy_ev) {
    static_cast<void>(electronTemperatureEv(value));
  }
  for (mesh::CellId cell = 0; cell < mean_energy_ev.mesh().numCells(); ++cell) {
    electron_temperature_ev[cell] = electronTemperatureEv(mean_energy_ev[cell]);
  }
}

}  // namespace pemu::physics
