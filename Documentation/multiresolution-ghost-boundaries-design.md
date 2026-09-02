# Ghost boundary particles at multiresolution interfaces — design

Status: agreed design (pending implementation approval)
Scope: `SolverMultiSetBlockMultiresolution` cases (2D 2-zone / 3-zone dam break, 3D dam break)

## 1. Problem

The multiresolution decomposition partitions the wall boundary over subdomains: each
subdomain owns the wall segments inside its box at its own resolution, and
`MultiresolutionBoundary` patches represent only the **neighboring fluid** (buffer
particles). A fluid particle near a subdomain corner (e.g. left-bottom edge of the L1/L2
zone) reaches with its kernel past the subdomain edge, but the wall continuing on the
other side belongs to the neighbor's boundary set and is invisible to it. Result:
one-sided wall support at the seam — asymmetric wall force, density dip/leak and corner
artifacts. The missing interaction is the **neighbor's wall** in the interface band.

Goal: make each fluid set interact with the complete wall next to its subdomain, while
keeping the structure clean and the model kernels untouched.

## 2. Variants

The whole pipeline (band selection, merge into the boundary set, index mapping, per-step
value refresh, unchanged fluid interaction) is shared; only *who manufactures the ghost
wall data* differs.

| stage | Variant B (implement first) | Variant A (prepared for, later) |
|---|---|---|
| ghost resolution | **neighbor's** resolution | **own** resolution |
| ghost geometry (r, n, ghostNodes, elementSize) | copied once from the neighbor's boundary particles in the band | own-dp lattice generated **in the init script** (`init_multi_resolution.py`) and loaded from VTK |
| ghost values (rho, gamma) | **copied** from the source boundary's fresh values, every step | **interpolated** from the source boundary's fresh values (MFD, same machinery as `interpolateVariables`) |
| location of construction | solver init (runtime) | init script (geometry) + solver init (mapping/interpolation) |
| wall representation near seam | mixed-dp (fine fluid over coarser wall) | uniform own-dp wall |

B delivers the physically decisive part (the wall is complete and continuous); A only
smooths wall sampling in the corner.

## 3. Geometry and band definition

- Ghost band = the **buffer band** of the directed interface: cells in
  `frameBack ⊖ frameFront` (outer patch: 1 own cell *outside* the frame; inner patch:
  1 own cell *inside* the frame, i.e. the shrunk ring). Same geometric band where the
  fluid buffer lives, by construction.
- For variant B, the band lies inside the **overlap ring** of the boundary set's own
  cell list (`dimsWithOverlap`, overlap = 1 cell), so binning/finding ghosts works on the
  existing grids with **no grid changes**.
- Band selection at solver init: cell-coordinate box test on the neighbor's boundary
  (`isInsideBox(gc, frameBack) != isInsideBox(gc, frameFront)` signalling band
  membership, matching the patch's `inner_overlap`/`outer_overlap` orientation).

## 4. Agreed data flow

1. **Init (once, solver side):**
   - For each directed interface (own → neighbor), select the neighbor's boundary
     particles in the band → append copies (`r`, `n`, `ghostNodes`, `elementSize`) at the
     **tail of the own boundary set** (`boundary[own]`).
   - Record a mapping per interface (`boundaryGhostParticles`):
     `{ ownIdx, neighborIdx, ghostBegin, ghostEnd, srcIndexList (device) }`.
     `srcIndexList` maps each ghost slot to its source particle index in
     `boundary[neighbor]` — used by the value refresh; also anchors the future variant-A
     interpolation stage.
   - Allocation: boundary allocations (`boundary_n_allocated`, factor ~2×) absorb the
     ghost band; hard assert otherwise.
2. **Every step, inside `interact()` (`SolverMultiSetBlockMultiresolution`):**
   - **Phase 1 — boundary updates:** `model.updateSolidBoundary(fluid_i, boundary_i)`
     for all sets (unchanged; its `forAll` skips overlap-ring particles, so ghost slots
     are *never* touched by it — intended).
   - **Phase 2 — ghost value refresh:** for each `boundaryGhostParticles` entry, set
     `rho`/`gamma` of ghost slots from the *source* boundary's fresh values:
     copy kernel (B) / interpolation kernel (A). This is the single seam where A will
     later replace B's copy with MFD interpolation using the same `srcIndexList`.
   - **Phase 3 — interactions (unchanged):** `model.interaction(fluid_i, boundary_i)`
     and `model.interactionWithOpenBoundary(fluid_i, patches...)`. The merged boundary
     set (own wall + ghosts) is treated uniformly — no extra fluid sweep, no zone
     plumbing.
3. **`multiresolutionUpdate()`** is *not* involved with ghost values (it is the
   fluid–fluid coupling, after `integrate`); ghost value refresh must happen in the same
   step's `interact()` or the wall forcing would be one step stale.

## 5. Why these semantics

- **Boundary values are properties of the neighboring subdomain's regular wall
  treatment**: the ghost's `rho/gamma` must come from the neighbor's boundary state
  (fresh from that subdomain's own `updateSolidBoundary`), never self-evaluated against
  the own fluid (circular, and it would erase exactly the shared boundary state).
- For BC types that require no per-step boundary update, variant B reduces to pure
  geometry — copied statics suffice; the refresh pass is a no-op.

## 6. Files to touch (variant B)

| file | change |
|---|---|
| `include/SPH/solvers/SolverMultiSetBlockMultiresolution.h` | `BoundaryGhostParticles` struct + per-interface list member |
| `include/SPH/solvers/SolverMultiSetBlockMultiresolution.hpp` | init: band selection + merge into boundary sets (+ assert on allocation); `interact()`: phase separation (update all → refresh ghosts → interactions); `save()`: unchanged (ghosts ride along in boundary output — cosmetics only) |
| `include/SPH/MultiresolutionRectangleBuffer.h` | none (band cells taken from the patch's public frame getters; optionally a small helper exposing the band test) |
| model kernels (`WCSPH_BI/*`) | none |
| example init scripts | none (unchanged for B) |

## 7. Variant A preparation (what changes later)

1. **Init script** (`init_multi_resolution.py`): emit, per subdomain, an own-dp ghost
   wall lattice for the band outside its box (walls + corners at the band ends, own
   ghost-node convention, own `elementSize`) — a clipped sub-generation exercising the
   existing boundary generator; save as an extra VTK set and load it at solver init as
   the ghost geometry (replaces B's copy-from-neighbor geometry).
2. **Solver**: same `boundaryGhostParticles` mapping (built by spatial match of the
   loaded ghost points against the source boundary, e.g. per-cell nearest search);
   refresh stage swaps the copy kernel for MFD interpolation from the source boundary
   (`rho`, `gamma`), mirroring `interpolateVariables` of the multiresolution buffer.
3. No changes to phases, merging, interactions, verification.

## 8. Safety notes verified in code

- `SolverMultiSetBase::removeParticlesOutOfDomain()` loops `fluidSets` only — boundary
  ghosts are not at risk in the block (non-MPI) path. `resetOverlaps()` calls it on
  `boundarySets` too, but that is the MPI/distributed path, unused by these cases
  (re-guard if enabled later).
- `updateSolidBoundary`'s `Particles::forAll` skips overlap-ring particles — ghosts in
  the ring are automatically excluded from the own-fluid evaluation ✅ (needed by design).
- Boundary sets are not re-sorted per step (boundary search runs at step 0 / on demand),
  so `srcIndexList` stays valid; add a re-derivation hook if boundary searches are ever
  enabled.

## 9. Limitations / watch points

- **Moving walls** are out of scope (ghosts are static at init); re-validation needed if
  a case with moving boundaries is used.
- Source-boundary values are computed from the *source subdomain's fluid*; where that
  fluid is absent (e.g. L1's wall under the L2 region while L1 is empty), the source
  wall value is weak there — self-consistent (the region is covered by the deeper
  refined set's own wall once populated), but corner values during transients deserve a
  look in verification.
- VTK output of `boundary_*` files will include ghost particles (duplicates along seams)
  — cosmetic.

## 10. Verification

1. **2D 3-zone dam break (t = 1.0):** corner behavior at L1/L2 left-bottom corners —
   density distribution continuous across the seam (no dip/leak line), pressure-sensor
   traces near corners vs. pre-change run, mass conservation at the previous level
   (≈ −1 %).
2. **3D 2-zone regression (t = 2.0):** unchanged away from walls; floor seam under the
   refined box gets identical coverage automatically.
3. Allocation asserts never fire; no CUDA errors; no per-step cost regression beyond a
   few hundred-copy kernel per step.
