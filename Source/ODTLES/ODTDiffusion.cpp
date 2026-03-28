#include "ODTDiffusion.H"

#include <array>
#include <cmath>
#include <vector>

#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_RealBox.H>

namespace pelec::odtles
{

namespace
{

using ConservativeCell = ODTLineState::ConservativeCell;
constexpr int NComp = static_cast<int>(ODTDiffusion::Component::NumComponents);

amrex::Real
getComponent(const ConservativeCell& c, int comp)
{
  AMREX_ALWAYS_ASSERT(comp >= 0 && comp < NComp);
  if (comp == 0) {
    return c.rho;
  }
  if (comp == 1) {
    return c.rhou;
  }
  if (comp == 2) {
    return c.rhov;
  }
  if (comp == 3) {
    return c.rhow;
  }
  return c.rhoE;
}

void
setComponent(ConservativeCell& c, int comp, amrex::Real val)
{
  AMREX_ALWAYS_ASSERT(comp >= 0 && comp < NComp);
  if (comp == 0) {
    c.rho = val;
  } else if (comp == 1) {
    c.rhou = val;
  } else if (comp == 2) {
    c.rhov = val;
  } else if (comp == 3) {
    c.rhow = val;
  } else {
    c.rhoE = val;
  }
}

void
solveTridiagonalInPlace(
  std::vector<amrex::Real>& a,
  std::vector<amrex::Real>& b,
  std::vector<amrex::Real>& c,
  std::vector<amrex::Real>& rhs,
  std::vector<amrex::Real>& x)
{
  const int n = static_cast<int>(rhs.size());
  AMREX_ALWAYS_ASSERT(n > 0);
  AMREX_ALWAYS_ASSERT(
    static_cast<int>(a.size()) == n && static_cast<int>(b.size()) == n &&
    static_cast<int>(c.size()) == n && static_cast<int>(x.size()) == n);

  for (int i = 1; i < n; ++i) {
    AMREX_ALWAYS_ASSERT(b[i - 1] != 0.0);
    const amrex::Real m = a[i] / b[i - 1];
    b[i] -= m * c[i - 1];
    rhs[i] -= m * rhs[i - 1];
  }

  AMREX_ALWAYS_ASSERT(b[n - 1] != 0.0);
  x[n - 1] = rhs[n - 1] / b[n - 1];
  for (int i = n - 2; i >= 0; --i) {
    AMREX_ALWAYS_ASSERT(b[i] != 0.0);
    x[i] = (rhs[i] - c[i] * x[i + 1]) / b[i];
  }
}

amrex::Real
maxAbsStateChange(const ODTLineState& a, const ODTLineState& b)
{
  AMREX_ALWAYS_ASSERT(a.numCells() == b.numCells());
  amrex::Real max_abs = 0.0;
  for (int i = 0; i < a.numCells(); ++i) {
    for (int comp = 0; comp < NComp; ++comp) {
      const amrex::Real da = getComponent(a.cell(i), comp);
      const amrex::Real db = getComponent(b.cell(i), comp);
      max_abs = std::max(max_abs, std::abs(db - da));
    }
  }
  return max_abs;
}

} // namespace

void
ODTDiffusion::applyImplicitUniform(
  const ODTLineGeometry& geom,
  ODTLineState& line_state,
  amrex::Real dt,
  const Controls& controls)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    geom.isDefined(), "ODTDiffusion requires defined line geometry");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    line_state.initialized() && line_state.valid(),
    "ODTDiffusion requires an initialized and valid line state");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    line_state.numCells() == geom.supportCellCount(),
    "ODTDiffusion requires state/geometry support-size consistency");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    dt >= 0.0, "ODTDiffusion requires non-negative dt");

  const int n = line_state.numCells();
  if (n <= 1 || dt == 0.0) {
    return;
  }

  const amrex::Real dx = geom.deltaS();
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dx > 0.0, "ODTDiffusion requires deltaS > 0");
  const amrex::Real inv_dx2 = 1.0 / (dx * dx);

  std::vector<amrex::Real> q_old(static_cast<std::size_t>(n), 0.0);
  std::vector<amrex::Real> a(static_cast<std::size_t>(n), 0.0);
  std::vector<amrex::Real> b(static_cast<std::size_t>(n), 0.0);
  std::vector<amrex::Real> c(static_cast<std::size_t>(n), 0.0);
  std::vector<amrex::Real> rhs(static_cast<std::size_t>(n), 0.0);
  std::vector<amrex::Real> q_new(static_cast<std::size_t>(n), 0.0);

  for (int comp = 0; comp < NComp; ++comp) {
    if (!controls.diffuse_component[static_cast<std::size_t>(comp)]) {
      continue;
    }

    const amrex::Real nu = controls.diffusivity[static_cast<std::size_t>(comp)];
    if (nu <= 0.0) {
      continue;
    }

    const amrex::Real r = dt * nu * inv_dx2;
    for (int i = 0; i < n; ++i) {
      q_old[static_cast<std::size_t>(i)] = getComponent(line_state.cell(i), comp);
      rhs[static_cast<std::size_t>(i)] = q_old[static_cast<std::size_t>(i)];

      // Zero-gradient (Neumann) closure at line ends:
      // ghost values mirror boundary-adjacent interior values.
      if (i == 0) {
        a[0] = 0.0;
        b[0] = 1.0 + r;
        c[0] = -r;
      } else if (i == n - 1) {
        a[static_cast<std::size_t>(i)] = -r;
        b[static_cast<std::size_t>(i)] = 1.0 + r;
        c[static_cast<std::size_t>(i)] = 0.0;
      } else {
        a[static_cast<std::size_t>(i)] = -r;
        b[static_cast<std::size_t>(i)] = 1.0 + 2.0 * r;
        c[static_cast<std::size_t>(i)] = -r;
      }
    }

    solveTridiagonalInPlace(a, b, c, rhs, q_new);

    for (int i = 0; i < n; ++i) {
      auto cell = line_state.cell(i);
      setComponent(cell, comp, q_new[static_cast<std::size_t>(i)]);
      line_state.setCell(i, cell);
    }
  }
}

ODTDiffusion::ValidationReport
ODTDiffusion::runMVPValidationHook()
{
  ValidationReport report{};

  const amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(2, 0, 0)));
  amrex::RealBox rb(
    AMREX_D_DECL(0.0, 0.0, 0.0), AMREX_D_DECL(3.0, 1.0, 1.0));
  int is_per[AMREX_SPACEDIM] = {AMREX_D_DECL(0, 0, 0)};
  const amrex::Geometry geom(domain, &rb, 0, is_per);
  const amrex::IntVect owner(AMREX_D_DECL(1, 0, 0));

  ODTLineGeometry line_geom(geom, 0, owner, 0);

  Controls ctrl{};
  ctrl.diffuse_component.fill(false);
  ctrl.diffuse_component[static_cast<int>(Component::RhoU)] = true;
  ctrl.diffusivity.fill(0.0);
  ctrl.diffusivity[static_cast<int>(Component::RhoU)] = 0.1;

  constexpr amrex::Real dt = 1.0e-2;
  constexpr amrex::Real tol = 1.0e-12;

  // Uniform-profile invariance check.
  ODTLineState uniform_state;
  uniform_state.initialize(line_geom);
  for (int i = 0; i < uniform_state.numCells(); ++i) {
    ODTLineState::ConservativeCell c{};
    c.rho = 1.0;
    c.rhou = 2.5;
    c.rhov = -0.4;
    c.rhow = 0.1;
    c.rhoE = 10.0;
    uniform_state.setCell(i, c);
  }
  uniform_state.setValid(true);
  ODTLineState uniform_before = uniform_state;
  applyImplicitUniform(line_geom, uniform_state, dt, ctrl);
  report.uniform_max_abs_change = maxAbsStateChange(uniform_before, uniform_state);
  report.uniform_invariant = (report.uniform_max_abs_change <= tol);

  // Nonuniform-profile evolution check.
  ODTLineState nonuniform_state;
  nonuniform_state.initialize(line_geom);
  for (int i = 0; i < nonuniform_state.numCells(); ++i) {
    ODTLineState::ConservativeCell c{};
    c.rho = 1.0;
    c.rhou = static_cast<amrex::Real>(i + 1);
    c.rhov = 0.0;
    c.rhow = 0.0;
    c.rhoE = 5.0;
    nonuniform_state.setCell(i, c);
  }
  nonuniform_state.setValid(true);
  ODTLineState nonuniform_before = nonuniform_state;
  applyImplicitUniform(line_geom, nonuniform_state, dt, ctrl);
  report.nonuniform_max_abs_change =
    maxAbsStateChange(nonuniform_before, nonuniform_state);
  report.nonuniform_evolved = (report.nonuniform_max_abs_change > tol);

  return report;
}

} // namespace pelec::odtles
