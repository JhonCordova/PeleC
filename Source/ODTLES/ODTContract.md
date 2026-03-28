# ODTLES Embedded SGS Contract (PeleC)

## Purpose
This document defines the repository-local engineering contract for embedded ODT SGS development in PeleC.

## Host-solver ownership
- PeleC is the host LES solver.
- PeleC owns timestep control, CFL control, AMR sequencing, conservative updates, source-term assembly, and resolved-flow evolution.
- ODT is subordinate to the LES timestep and participates only as an embedded SGS closure path.

## Architecture boundary
- Standalone `ODTLES/ODT` is algorithmic reference only.
- The target architecture is native PeleC integration, not a second standalone solver architecture.
- Conservative compatibility with PeleC flux/source plumbing is required.

## Current local ownership identity
- Persistent local ODT entry identity is:
  - `(level, owner_cell, dir)`
- One owner cell may have up to three directional entries per level.

## Current responsibility split
- `ODTLineGeometry`: local 1D support geometry metadata and owner-cell interval metadata.
- `ODTLineState`: conservative 1D line state storage (`rho`, `rhou`, `rhov`, `rhow`, `rhoE`) plus derived helper accessors.
- `ODTManager`: persistent ownership/lifecycle and lookup/create APIs for `(level, owner_cell, dir)` entries.

## MVP execution priority
1. Persistent local infrastructure (identity, geometry, conservative line state, manager lifecycle).
2. LES -> ODT reconciliation.
3. ODT stepping/subordinate local advancement.
4. SGS moment/flux extraction and host-side conservative coupling.

## Explicitly out of current T1-T3 scope
- LES -> ODT reconciliation logic.
- ODT stepping/event logic.
- SGS extraction/closure return logic.
- Restart/regrid migration of persistent entries.
- Any standalone-style global domain/solver/driver architecture.

## Contributor rule
All ODTLES production changes must remain in `ODTLES/PeleC` and preserve PeleC host ownership and conservative integration semantics.
