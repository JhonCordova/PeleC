#ifndef PELEC_ODTLES_ODTPARAMS_H
#define PELEC_ODTLES_ODTPARAMS_H

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
};

} // namespace pelec::odtles

#endif
