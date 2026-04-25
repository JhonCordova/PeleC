#include "ODTThermoBridge.H"

#include <algorithm>
#include <cmath>

#include "PelePhysics.H"

namespace pelec::odtles
{

ODTThermoBridge::RecoveryControls
ODTThermoBridge::defaultRecoveryControls()
{
  return RecoveryControls{};
}

const char*
ODTThermoBridge::failureModeString(FailureMode mode)
{
  switch (mode) {
  case FailureMode::None:
    return "none";
  case FailureMode::CellIndexOutOfRange:
    return "cell index out of range";
  case FailureMode::DensityBelowFloor:
    return "density below floor";
  case FailureMode::SpeciesSizeMismatch:
    return "species container size mismatch";
  case FailureMode::NonPositiveSpeciesSum:
    return "non-positive species conservative sum";
  case FailureMode::InternalEnergyBelowFloor:
    return "specific internal energy below floor";
  case FailureMode::EOSNonPhysicalTemperature:
    return "EOS returned non-physical temperature";
  default:
    return "unknown thermodynamic recovery failure";
  }
}

ODTThermoBridge::RecoveredState
ODTThermoBridge::recoverCellThermoState(
  const ODTLineState& line_state,
  int i,
  const RecoveryControls& ctrl)
{
  if (i < 0 || i >= line_state.numCells()) {
    RecoveredState out;
    out.failure_mode = FailureMode::CellIndexOutOfRange;
    return out;
  }
  return recoverCellThermoState(line_state.cell(i), ctrl);
}

ODTThermoBridge::RecoveredState
ODTThermoBridge::recoverCellThermoState(
  const ConservativeCell& cell,
  const RecoveryControls& ctrl)
{
  RecoveredState out;
  out.rho = cell.rho;

  if (out.rho <= ctrl.rho_floor) {
    out.failure_mode = FailureMode::DensityBelowFloor;
    return out;
  }

  if (cell.rhoY.size() != static_cast<std::size_t>(NUM_SPECIES)) {
    out.failure_mode = FailureMode::SpeciesSizeMismatch;
    return out;
  }

  const amrex::Real rho_inv = 1.0 / out.rho;
  out.velocity = {
    AMREX_D_DECL(cell.rhou * rho_inv, cell.rhov * rho_inv, cell.rhow * rho_inv)};
  out.kinetic_energy = 0.5 * AMREX_D_TERM(
                               out.velocity[0] * out.velocity[0], +
                               out.velocity[1] * out.velocity[1], +
                               out.velocity[2] * out.velocity[2]);

  amrex::Real species_sum = 0.0;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    const amrex::Real rhoY_n =
      std::max<amrex::Real>(cell.rhoY[static_cast<std::size_t>(n)], 0.0);
    out.rhoY[static_cast<std::size_t>(n)] = rhoY_n;
    species_sum += rhoY_n;
  }

  if (species_sum <= 0.0) {
    out.failure_mode = FailureMode::NonPositiveSpeciesSum;
    return out;
  }

  const amrex::Real species_scale = out.rho / species_sum;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    out.rhoY[static_cast<std::size_t>(n)] *= species_scale;
    out.Y[static_cast<std::size_t>(n)] = out.rhoY[static_cast<std::size_t>(n)] * rho_inv;
  }

  out.specific_internal_energy = cell.rhoE * rho_inv - out.kinetic_energy;
  if (
    out.specific_internal_energy <= ctrl.e_floor ||
    !std::isfinite(out.specific_internal_energy)) {
    out.failure_mode = FailureMode::InternalEnergyBelowFloor;
    return out;
  }

  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real T = 0.0;
  eos.REY2T(out.rho, out.specific_internal_energy, out.Y.data(), T);

  if (!std::isfinite(T) || T <= ctrl.T_floor) {
    out.failure_mode = FailureMode::EOSNonPhysicalTemperature;
    return out;
  }

  out.temperature = T;
  out.success = true;
  out.failure_mode = FailureMode::None;
  return out;
}

} // namespace pelec::odtles
