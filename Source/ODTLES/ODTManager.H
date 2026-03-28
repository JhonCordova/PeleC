#ifndef PELEC_ODTLES_ODTMANAGER_H
#define PELEC_ODTLES_ODTMANAGER_H

#include <cstddef>
#include <map>
#include <vector>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

#include "ODTLineGeometry.H"
#include "ODTLineState.H"
#include "ODTParams.H"

namespace pelec::odtles
{

class ODTManager
{
public:
  struct LineEntry
  {
    ODTLineGeometry geometry;
    ODTLineState state;
  };

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

  bool hasLineEntry(int level, const amrex::IntVect& owner_cell, int dir) const;
  LineEntry*
  findLineEntry(int level, const amrex::IntVect& owner_cell, int dir);
  const LineEntry*
  findLineEntry(int level, const amrex::IntVect& owner_cell, int dir) const;
  LineEntry& getOrCreateLineEntry(
    int level,
    const amrex::IntVect& owner_cell,
    int dir,
    const amrex::Geometry& geom);

  std::size_t numLineEntries() const noexcept { return m_line_entries.size(); }

  ODTLineGeometry& getOrCreateLineGeometry(
    int level,
    const amrex::IntVect& owner_cell,
    int dir,
    const amrex::Geometry& geom);
  ODTLineGeometry*
  findLineGeometry(int level, const amrex::IntVect& owner_cell, int dir);
  const ODTLineGeometry* findLineGeometry(
    int level,
    const amrex::IntVect& owner_cell,
    int dir) const;

  std::size_t numLineGeometries() const noexcept
  {
    return m_line_entries.size();
  }

  ODTLineState&
  getOrCreateLineState(
    int level,
    const amrex::IntVect& owner_cell,
    int dir,
    const amrex::Geometry& geom);
  // Host-side input-contract path: assemble local support-cell conservative
  // means directly from valid LES cells in state and initialize first-time
  // line state.
  ODTLineState& initializeLineStateFromLESState(
    int level,
    const amrex::IntVect& owner_cell,
    int dir,
    const amrex::Geometry& geom,
    const amrex::MultiFab& state);
  ODTLineState* findLineState(
    int level,
    const amrex::IntVect& owner_cell,
    int dir);
  const ODTLineState* findLineState(
    int level,
    const amrex::IntVect& owner_cell,
    int dir) const;

  std::size_t numLineStates() const noexcept { return m_line_entries.size(); }

private:
  struct LineKey
  {
    int level = -1;
    amrex::IntVect owner_cell{AMREX_D_DECL(0, 0, 0)};
    int dir = 0;

    bool operator<(const LineKey& other) const noexcept;
  };

  void assertLevelScope(int level) const;
  std::vector<ODTLineState::ConservativeCell> assembleLocalSupportCellAverages(
    const ODTLineGeometry& line_geom,
    const amrex::MultiFab& state) const;
  ODTLineState& initializeLineStateFromLESSupportAverages(
    int level,
    const amrex::IntVect& owner_cell,
    int dir,
    const amrex::Geometry& geom,
    const std::vector<ODTLineState::ConservativeCell>&
      local_support_cell_averages);

  ODTParams m_params{};
  bool m_initialized = false;
  int m_level = -1;
  std::map<LineKey, LineEntry> m_line_entries;
};

} // namespace pelec::odtles

#endif
