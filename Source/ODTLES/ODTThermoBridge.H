#ifndef PELEC_ODTLES_ODTTHERMOBRIDGE_H
#define PELEC_ODTLES_ODTTHERMOBRIDGE_H

#include <array>
#include <vector>

#include <AMReX_SPACE.H>

#include "ODTLineState.H"
#include "mechanism.H"

namespace pelec::odtles
{

class ODTThermoBridge
{
public:
  using ConservativeCell = ODTLineState::ConservativeCell;

  struct RecoveryControls
  {
    amrex::Real rho_floor = 1.0e-14;
    amrex::Real e_floor = 1.0e-12;
    amrex::Real T_floor = 1.0e-12;
  };

  enum class FailureMode : int
  {
    None = 0,
    CellIndexOutOfRange = 1,
    DensityBelowFloor = 2,
    SpeciesSizeMismatch = 3,
    NonPositiveSpeciesSum = 4,
    InternalEnergyBelowFloor = 5,
    EOSNonPhysicalTemperature = 6
  };

  struct RecoveredState
  {
    bool success = false;
    FailureMode failure_mode = FailureMode::None;
    std::array<amrex::Real, AMREX_SPACEDIM> velocity{
      AMREX_D_DECL(0.0, 0.0, 0.0)};
    std::vector<amrex::Real> rhoY;
    std::vector<amrex::Real> Y;
    amrex::Real rho = 0.0;
    amrex::Real kinetic_energy = 0.0;
    amrex::Real specific_internal_energy = 0.0;
    amrex::Real temperature = 0.0;

    RecoveredState()
      : rhoY(static_cast<std::size_t>(NUM_SPECIES), 0.0),
        Y(static_cast<std::size_t>(NUM_SPECIES), 0.0)
    {
    }
  };

  static RecoveryControls defaultRecoveryControls();
  static const char* failureModeString(FailureMode mode);

  static RecoveredState recoverCellThermoState(
    const ODTLineState& line_state,
    int i,
    const RecoveryControls& ctrl = defaultRecoveryControls());
  static RecoveredState recoverCellThermoState(
    const ConservativeCell& cell,
    const RecoveryControls& ctrl = defaultRecoveryControls());
};

} // namespace pelec::odtles

#endif
