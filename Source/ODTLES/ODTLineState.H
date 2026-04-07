#ifndef PELEC_ODTLES_ODTLINESTATE_H
#define PELEC_ODTLES_ODTLINESTATE_H

#include <vector>

#include <AMReX_REAL.H>

#include "ODTLineGeometry.H"

namespace pelec::odtles
{

class ODTLineState
{
public:
  struct ConservativeCell
  {
    amrex::Real rho = 0.0;
    amrex::Real rhou = 0.0;
    amrex::Real rhov = 0.0;
    amrex::Real rhow = 0.0;
    amrex::Real rhoE = 0.0;
  };

  struct ThermoInput
  {
    amrex::Real rho = 0.0;
    amrex::Real e_int = 0.0;
  };

  ODTLineState() = default;

  void initialize(const ODTLineGeometry& geom);
  // Validation-only utility for local MVP hooks.
  // This intentionally bypasses geometry-derived support construction and
  // should not be used by runtime ODT integration paths.
  void initializeForValidation(int n_cells);
  void reset();
  void setValid(bool valid) noexcept { m_valid = valid; }

  bool initialized() const noexcept { return m_initialized; }
  bool valid() const noexcept { return m_valid; }
  int numCells() const noexcept { return static_cast<int>(m_cells.size()); }

  const ConservativeCell& cell(int i) const;
  ConservativeCell& cell(int i);
  void setCell(int i, const ConservativeCell& state);

  // Velocity helpers derived on demand from conservative state.
  amrex::Real velocityU(int i) const;
  amrex::Real velocityV(int i) const;
  amrex::Real velocityW(int i) const;

  // Internal-energy helpers derived on demand from conservative state.
  amrex::Real rhoInternalEnergy(int i) const;
  amrex::Real specificInternalEnergy(int i) const;

  // Hook payload for future EOS integration (no persistent primitives stored).
  ThermoInput thermoInput(int i) const;

private:
  std::vector<ConservativeCell> m_cells;
  bool m_valid = false;
  bool m_initialized = false;
};

} // namespace pelec::odtles

#endif
