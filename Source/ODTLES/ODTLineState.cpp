#include "ODTLineState.H"

#include <array>
#include <stdexcept>

#include <AMReX.H>

namespace pelec::odtles
{

namespace
{

void
checkCellIndex(int i, int n)
{
  if (i < 0 || i >= n) {
    throw std::out_of_range("ODTLineState cell index out of range");
  }
}

void
checkSpeciesIndex(int n)
{
  if (n < 0 || n >= NUM_SPECIES) {
    throw std::out_of_range("ODTLineState species index out of range");
  }
}

} // namespace

void
ODTLineState::initialize(const ODTLineGeometry& geom)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    geom.isDefined(), "ODTLineState::initialize requires a defined geometry");

  const int n = geom.supportCellCount();
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    n > 0, "ODTLineState::initialize requires non-empty support");

  m_cells.assign(static_cast<std::size_t>(n), ConservativeCell{});
  m_valid = false;
  m_initialized = true;
}

void
ODTLineState::initializeForValidation(int n_cells)
{
  // Test-hook path only: allows validation code to exercise mapping logic on
  // controlled synthetic line lengths.
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    n_cells > 0,
    "ODTLineState::initializeForValidation requires a positive cell count");

  m_cells.assign(static_cast<std::size_t>(n_cells), ConservativeCell{});
  m_valid = false;
  m_initialized = true;
}

void
ODTLineState::reset()
{
  m_cells.clear();
  m_valid = false;
  m_initialized = false;
}

const ODTLineState::ConservativeCell&
ODTLineState::cell(int i) const
{
  checkCellIndex(i, numCells());
  return m_cells[static_cast<std::size_t>(i)];
}

ODTLineState::ConservativeCell&
ODTLineState::cell(int i)
{
  checkCellIndex(i, numCells());
  return m_cells[static_cast<std::size_t>(i)];
}

void
ODTLineState::setCell(int i, const ConservativeCell& state)
{
  checkCellIndex(i, numCells());
  m_cells[static_cast<std::size_t>(i)] = state;
}

amrex::Real
ODTLineState::velocityU(int i) const
{
  const auto& s = cell(i);
  return (s.rho != 0.0) ? (s.rhou / s.rho) : 0.0;
}

amrex::Real
ODTLineState::velocityV(int i) const
{
  const auto& s = cell(i);
  return (s.rho != 0.0) ? (s.rhov / s.rho) : 0.0;
}

amrex::Real
ODTLineState::velocityW(int i) const
{
  const auto& s = cell(i);
  return (s.rho != 0.0) ? (s.rhow / s.rho) : 0.0;
}

std::array<amrex::Real, AMREX_SPACEDIM>
ODTLineState::velocities(int i) const
{
  return {AMREX_D_DECL(velocityU(i), velocityV(i), velocityW(i))};
}

amrex::Real
ODTLineState::speciesMassFraction(int i, int n) const
{
  checkSpeciesIndex(n);
  const auto& s = cell(i);
  if (s.rho == 0.0 || s.rhoY.size() != static_cast<std::size_t>(NUM_SPECIES)) {
    return 0.0;
  }
  return s.rhoY[static_cast<std::size_t>(n)] / s.rho;
}

std::vector<amrex::Real>
ODTLineState::speciesMassFractions(int i) const
{
  std::vector<amrex::Real> Y(static_cast<std::size_t>(NUM_SPECIES), 0.0);
  if (!hasExpectedSpeciesContainerSize(i)) {
    return Y;
  }
  for (int n = 0; n < NUM_SPECIES; ++n) {
    Y[static_cast<std::size_t>(n)] = speciesMassFraction(i, n);
  }
  return Y;
}

bool
ODTLineState::hasExpectedSpeciesContainerSize(int i) const
{
  const auto& s = cell(i);
  return s.rhoY.size() == static_cast<std::size_t>(NUM_SPECIES);
}

bool
ODTLineState::isCellAdmissible(int i) const
{
  const auto& s = cell(i);
  if (s.rho <= 0.0) {
    return false;
  }
  if (!hasExpectedSpeciesContainerSize(i)) {
    return false;
  }
  for (int n = 0; n < NUM_SPECIES; ++n) {
    if (s.rhoY[static_cast<std::size_t>(n)] < 0.0) {
      return false;
    }
  }
  return true;
}

amrex::Real
ODTLineState::rhoInternalEnergy(int i) const
{
  const auto& s = cell(i);
  if (s.rho == 0.0) {
    return 0.0;
  }

  const amrex::Real inv_rho = 1.0 / s.rho;
  const amrex::Real kinetic =
    0.5 * inv_rho * (s.rhou * s.rhou + s.rhov * s.rhov + s.rhow * s.rhow);
  return s.rhoE - kinetic;
}

amrex::Real
ODTLineState::specificInternalEnergy(int i) const
{
  const auto& s = cell(i);
  return (s.rho != 0.0) ? (rhoInternalEnergy(i) / s.rho) : 0.0;
}

ODTLineState::ThermoInput
ODTLineState::thermoInput(int i) const
{
  const auto& s = cell(i);
  return ThermoInput{s.rho, specificInternalEnergy(i)};
}

} // namespace pelec::odtles
