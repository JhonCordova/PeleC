#ifndef PELEC_ODTLES_ODTPARAMS_H
#define PELEC_ODTLES_ODTPARAMS_H

#include <AMReX_REAL.H>

namespace pelec::odtles
{

struct ODTParams
{
  // Runtime enable of the ODT LES model path in host coupling.
  bool enabled = false;
  // Max local ODT substeps per host LES step.
  int max_local_substeps = 1;
  // Number of ODT subsegments per host LES cell on local support lines.
  int subsegments_per_host_cell = 3;
  // Poisson candidate sampling rate for local ODT eddy trials inside dt_LES.
  // Set to zero to keep diffusion-only local stepping.
  amrex::Real event_rate = 0.0;
  // Optional runtime extension: allow support entries sourced from host
  // coarse-fine fills when no boundary/invalid provenance is present.
  // Default remains strict same-level-only acceptance.
  bool allow_amr_coarse_fine_support = false;
};

} // namespace pelec::odtles

#endif
