#ifndef PELEC_ODTLES_ODTLINESTATE_H
#define PELEC_ODTLES_ODTLINESTATE_H

#include <AMReX_IntVect.H>

namespace pelec::odtles
{

class ODTLineState
{
public:
  ODTLineState() = default;
  ODTLineState(const amrex::IntVect& owner_cell, int line_dir);

  void reset();

  const amrex::IntVect& ownerCell() const noexcept { return m_owner_cell; }
  int lineDir() const noexcept { return m_line_dir; }
  bool valid() const noexcept { return m_valid; }

private:
  amrex::IntVect m_owner_cell{AMREX_D_DECL(0, 0, 0)};
  int m_line_dir = 0;
  bool m_valid = false;
};

} // namespace pelec::odtles

#endif
