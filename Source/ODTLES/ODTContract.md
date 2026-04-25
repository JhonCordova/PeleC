# ODTLES Embedded SGS Contract (PeleC)

## Purpose
This document is the T1 software-and-coupling contract for embedded ODT SGS in
PeleC. It defines what is active now, what is deferred by phase, and what is
explicitly excluded.

## Host-solver ownership (non-negotiable)
- PeleC is the only host LES solver.
- PeleC owns global timestep control (`dt_LES`), CFL control, AMR sequencing,
  conservative update assembly, and resolved-flow evolution.
- ODT is a subordinate local SGS module invoked only from the existing PeleC
  LES source path.
- ODT does not own a global advance loop and does not maintain a second
  resolved-field state.

## Software-contract touchpoints in PeleC
The embedded ODT contract is anchored to these host files and locations:
- `Source/Sources.cpp`: source selection/dispatch path using `les_src`.
- `Source/LES.cpp`: `getLESTerm()` dispatch and ODT LES-term assembly path
  (`getODTLESTerm()`).
- `Source/LES.H`: LES model interfaces and flux/source helpers used by LES-term
  assembly.
- `Source/PeleC.H`: LES model enum, LES APIs, and ODT manager sidecar
  ownership (`odt_manager`).
- `Source/ODTLES`: ODT local module implementation
  (manager/geometry/state/reconcile/stepper/moment extraction).

## ODT state ownership and identity
- ODT persistent state is sidecar-owned by `ODTManager` and is outside the main
  conserved `MultiFab`.
- Persistent identity is `(level, owner_cell, dir)`.
- One owner cell can own up to three directional local lines per level.

## LES -> ODT -> LES contract (active path)
1. LES -> ODT line reconciliation:
   same-level accepted support data are used to initialize/reconcile each local
   directional line for the owner cell.
2. ODT local advancement:
   each line advances only within host-provided `dt_LES`; local subcycling is
   subordinate and closes `dt_LES` exactly.
3. ODT -> LES SGS return:
   ODT returns directional SGS momentum contributions derived from moments over
   the owner-cell central interval, not corrected LES fields.
4. Host conservative assembly:
   PeleC maps the directional returns to LES fluxes and applies the existing
   flux-divergence route in the standard LES source assembly path.

## Owner-cell central interval and Favre moments
- SGS moments are extracted over the owner-cell central interval of each local
  directional line.
- Favre-consistent moments are used to define
  `tau_ij = <rho u_i u_j> - <rho u_i><rho u_j>/<rho>`.
- Directional momentum SGS return uses `tau_ij` and is mapped through the
  existing PeleC LES flux/divergence machinery.

## Phase scope freeze (MVP and extensions)
- Phase-1 (active): momentum SGS only via `tau_ij`.
- Phase-1 (active): `UEDEN` contribution is explicitly kept zero.
- Phase-2 (deferred): add total-energy SGS work contribution.
- Phase-3 (deferred): add SGS heat flux `Q_j`.
- `J_j` and `D_j` are explicitly deferred (not in Phase-1/2/3 MVP closure
  package).

## Runtime ODT-Line Molecular Model (production standard)
- Runtime ODT line momentum diffusion is unambiguous in production:
  - thermochemical admissibility/EOS recovery is always active for transport
    queries,
  - host molecular viscosity is always used,
  - momentum diffusion form is
    `d(rho*u_i)/dt = d/ds(mu * d(u_i)/ds)`, `u_i=(rho*u_i)/rho`.
- User-facing toggles for enabling/disabling or selecting alternate molecular
  momentum forms are deprecated and ignored in production runtime.
- The only runtime ODT local-step control retained is:
  - `pelec.odt_max_local_substeps`
- Thermochemical recovery remains a local admissibility/EOS bridge only;
  persistent ODT line storage remains conservative (`rho`, momentum, `rhoE`,
  `rhoY_k`) with no persistent temperature variable.

Minimal input example:
```ini
pelec.do_les = 1
pelec.les_model = 4
pelec.odt_max_local_substeps = 1
```

## Explicit exclusions for this contract
- No patch-based McDermott/XLES-style coupling in this package.
- No direct overwrite/correction of resolved LES fields from ODT.
- No second global solver/driver architecture ported from standalone ODT.
- No data-layout migration that embeds ODT state inside the main conserved
  `MultiFab`.

## T1 closure statement
For T1 documentation closure, this repository defines ODT as an embedded,
host-owned, conservative SGS path through existing PeleC LES source plumbing,
with Phase-1 restricted to momentum SGS and `UEDEN = 0`.
