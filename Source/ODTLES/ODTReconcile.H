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
  struct ReconcileControls
  {
    amrex::Real rho_floor = 1.0e-14;
    amrex::Real e_floor = 1.0e-12;
    amrex::Real alpha_min = 1.0e-6;
    int max_alpha_reductions = 40;
  };

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

  // Reconcile an existing valid line by enforcing the current owner-cell LES
  // conservative mean while preserving old subgrid residual content when
  // admissible:
  //   trial = target_mean + alpha * old_residual
  // with alpha starting at 1 and reduced uniformly if needed to satisfy
  // admissibility floors.
  // This path is distinct from first-time initialization and does not
  // reconstruct/reset the line from scratch.
  static ReconcileControls defaultReconcileControls();
  static void reconcileExistingLineStateToOwnerAverage(
    const ODTLineGeometry& geom,
    const ConservativeCell& owner_cell_average,
    ODTLineState& line_state);
  static void reconcileExistingLineStateToOwnerAverage(
    const ODTLineGeometry& geom,
    const ConservativeCell& owner_cell_average,
    ODTLineState& line_state,
    const ReconcileControls& ctrl);
};

} // namespace pelec::odtles

#endif
