#ifndef PELEC_ODTLES_ODTRECONCILE_H
#define PELEC_ODTLES_ODTRECONCILE_H

#include <vector>

#include "ODTLineGeometry.H"
#include "ODTLineState.H"

namespace pelec::odtles
{

class ODTReconcile
{
public:
  using ConservativeCell = ODTLineState::ConservativeCell;

  // First-time LES -> ODT initialization on local support.
  //
  // Contract for input averages:
  // - One conservative LES cell-average per support segment in geometry order.
  // - Values must correspond to valid local LES support cells associated with
  //   the ODT line support, not ghost-derived physical means.
  // - This layer performs no ghost-cell physical reconstruction.
  static void initializeLineStateFromLESSupportAverages(
    const ODTLineGeometry& geom,
    const std::vector<ConservativeCell>& local_support_cell_averages,
    ODTLineState& line_state);
};

} // namespace pelec::odtles

#endif
