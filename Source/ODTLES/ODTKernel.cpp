#include "ODTKernel.H"

#include <vector>

#include <AMReX.H>

namespace pelec::odtles
{

bool
ODTKernel::applyVelocityKernelZeroNet(
  ODTLineState& line_state,
  int i_lo,
  int i_hi,
  const VelocityKernelSpec& kernel)
{
  if (!kernel.enabled) {
    return false;
  }

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    line_state.initialized() && line_state.valid(),
    "ODTKernel requires initialized, valid line state");
  AMREX_ALWAYS_ASSERT(i_lo >= 0 && i_hi < line_state.numCells() && i_lo <= i_hi);

  const int n = i_hi - i_lo + 1;
  std::vector<amrex::Real> shape(static_cast<std::size_t>(n), 1.0);
  if (!kernel.shape.empty()) {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      static_cast<int>(kernel.shape.size()) == n,
      "ODTKernel shape size must match mapped interval size");
    shape = kernel.shape;
  }

  amrex::Real rho_sum = 0.0;
  amrex::Real rho_shape_sum = 0.0;
  for (int i = 0; i < n; ++i) {
    const auto& c = line_state.cell(i_lo + i);
    rho_sum += c.rho;
    rho_shape_sum += c.rho * shape[static_cast<std::size_t>(i)];
  }
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    rho_sum > 0.0, "ODTKernel requires positive total density in interval");

  const amrex::Real shape_mean_rho = rho_shape_sum / rho_sum;

  for (int i = 0; i < n; ++i) {
    auto c = line_state.cell(i_lo + i);
    const amrex::Real centered_shape =
      shape[static_cast<std::size_t>(i)] - shape_mean_rho;

    // Velocity-kernel hook: adds rho*delta(u,v,w), zero-net momentum by
    // construction over the mapped interval; rho and rhoE untouched.
    c.rhou += c.rho * kernel.amp_u * centered_shape;
    c.rhov += c.rho * kernel.amp_v * centered_shape;
    c.rhow += c.rho * kernel.amp_w * centered_shape;
    line_state.setCell(i_lo + i, c);
  }

  return true;
}

} // namespace pelec::odtles

