#include "ODTLineState.H"

namespace pelec::odtles
{

ODTLineState::ODTLineState(const amrex::IntVect& owner_cell, int line_dir)
  : m_owner_cell(owner_cell), m_line_dir(line_dir), m_valid(true)
{
}

void
ODTLineState::reset()
{
  m_valid = false;
}

} // namespace pelec::odtles
