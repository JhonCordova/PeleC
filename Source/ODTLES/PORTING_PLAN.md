# ODTLES Porting Boundary Plan (PeleC)

## Scope
This plan defines what can be reused conceptually from standalone ODT and what must be implemented natively for embedded PeleC ODTLES.

## May be adapted conceptually from standalone ODT
- 1D line-evolution concepts for local directional lines.
- Eddy-event concepts and triplet-map algorithmic ideas.
- 1D conservative-state handling patterns.
- Diffusion/evolution algorithm references.

## Must be reimplemented natively in PeleC
- Data ownership and lifecycle under PeleC level/cell/direction identity `(level, owner_cell, dir)`.
- Conservative coupling through PeleC LES source/flux-divergence pathway.
- LES-controlled timing and call order (host-controlled sequencing).
- Manager/geometry/state APIs aligned with PeleC host responsibilities.

## Must NOT be copied from standalone ODT
- Standalone global solver loop or time-integration ownership.
- Standalone domain/driver/parser/processor architecture.
- Any second resolved-flow architecture parallel to PeleC.
- Any coupling pattern that replaces PeleC resolved fields directly instead of returning SGS-consistent conservative contributions.

## Current phase boundary
Current repository state (through T3) is infrastructure-first only:
- persistent ownership manager,
- local line geometry,
- conservative line state.

Not yet implemented:
- LES -> ODT reconciliation,
- ODT stepping/event execution,
- SGS extraction return path,
- restart/regrid migration.

## Implementation discipline
- Keep all production changes inside `ODTLES/PeleC`.
- Use `ODTLES/ODT` as reference input only.
- Preserve PeleC host-solver ownership and conservative consistency.
