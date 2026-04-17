# ODTLES Porting Boundary Plan (PeleC)

## Scope
This plan defines what can be reused conceptually from standalone ODT and what must be implemented natively for embedded PeleC ODTLES.

## May be adapted conceptually from standalone ODT
- 1D line-evolution concepts for local directional lines.
- Eddy-event concepts and triplet-map algorithmic ideas.
- 1D conservative-state handling patterns.
- Diffusion/evolution algorithm references.
- Favre-moment extraction concepts for SGS stress construction.

## Must be reimplemented natively in PeleC
- Data ownership and lifecycle under PeleC level/cell/direction identity `(level, owner_cell, dir)`.
- Conservative coupling through PeleC LES source/flux-divergence pathway.
- LES-controlled timing and call order (host-controlled sequencing).
- Manager/geometry/state APIs aligned with PeleC host responsibilities.
- Sidecar ownership of ODT state via `Source/ODTLES` and host integration through:
  - `Source/Sources.cpp`,
  - `Source/LES.cpp`,
  - `Source/LES.H`,
  - `Source/PeleC.H`.

## Must NOT be copied from standalone ODT
- Standalone global solver loop or time-integration ownership.
- Standalone domain/driver/parser/processor architecture.
- Any second resolved-flow architecture parallel to PeleC.
- Any coupling pattern that replaces PeleC resolved fields directly instead of returning SGS-consistent conservative contributions.
- Any patch-based McDermott/XLES-style coupling pattern for this embedded path.

## MVP mathematical scope freeze
- Phase-1 (active): momentum SGS closure through `tau_ij` only.
- Phase-1 (active): returned energy-slot contribution remains zero (`UEDEN = 0`).
- Phase-2 (deferred): add total-energy SGS work contribution.
- Phase-3 (deferred): add SGS heat flux `Q_j`.
- Explicitly deferred terms: `J_j` and `D_j`.

## Runtime contract status for current repository
The active embedded path includes:
- LES -> ODT reconciliation on owner-cell directional lines,
- local ODT stepping subordinate to `dt_LES` with exact local closure of `dt_LES`,
- owner-central interval Favre-moment extraction and directional `tau_ij`-based return,
- conservative return through existing PeleC LES flux/divergence source assembly.

## Implementation discipline
- Keep all production changes inside `ODTLES/PeleC`.
- Use `ODTLES/ODT` as reference input only.
- Preserve PeleC host-solver ownership and conservative consistency.
