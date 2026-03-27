#ifndef PELEC_ODTLES_ODTPARAMS_H
#define PELEC_ODTLES_ODTPARAMS_H

namespace pelec::odtles
{

struct ODTParams
{
  bool enabled = false;
  int max_local_substeps = 1;
};

} // namespace pelec::odtles

#endif
