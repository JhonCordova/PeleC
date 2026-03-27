#ifndef PELEC_ODTLES_ODTMANAGER_H
#define PELEC_ODTLES_ODTMANAGER_H

#include <cstddef>
#include <map>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_IntVect.H>

#include "ODTLineState.H"
#include "ODTParams.H"

namespace pelec::odtles
{

class ODTManager
{
public:
  ODTManager() = default;
  explicit ODTManager(const ODTParams& params);

  void setParams(const ODTParams& params);
  const ODTParams& params() const noexcept { return m_params; }

  void initializeLevel(
    int level,
    const amrex::BoxArray& grids,
    const amrex::DistributionMapping& dmap);
  bool isInitialized() const noexcept { return m_initialized; }
  int level() const noexcept { return m_level; }

  void clear();

  ODTLineState& getOrCreateLineState(const amrex::IntVect& owner_cell, int dir);
  ODTLineState* findLineState(const amrex::IntVect& owner_cell, int dir);
  const ODTLineState* findLineState(
    const amrex::IntVect& owner_cell,
    int dir) const;

  std::size_t numLineStates() const noexcept { return m_line_states.size(); }

private:
  struct LineKey
  {
    amrex::IntVect owner_cell{AMREX_D_DECL(0, 0, 0)};
    int dir = 0;

    bool operator<(const LineKey& other) const noexcept;
  };

  ODTParams m_params{};
  bool m_initialized = false;
  int m_level = -1;
  std::map<LineKey, ODTLineState> m_line_states;
};

} // namespace pelec::odtles

#endif
