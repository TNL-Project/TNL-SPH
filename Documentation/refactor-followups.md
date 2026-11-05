<!-- Moved from the refactor session's .omo registry on 2026-11-05; entries were deferred from 2026-10-22 onward. -->

# Refactor follow-ups (deferred work registry)

Items parked during the `tnl-sph-refactor-tasks` walkthrough (`.omo/plans/tnl-sph-refactor-tasks.md`),
to be executed as separate follow-up work. Each entry is self-contained: verified inventory + approach.

WORKER NOTE: on first execution of any follow-up item, move/merge this file's content into
`Documentation/refactor-followups.md` (repo location; T24 of the plan may already have created it —
append there, keeping this section) and use that repo file as the single registry from then on.

---

## DEFERRED-T15 — Adopt TNL's assert/exception error-handling idiom

Deferred 2026-10-22 at user request during the plan walkthrough. Findings stay valid (verified counts
below); line numbers are pre-`tnl-sph-refactor-tasks` anchors — re-locate by grep content.

### Finding

TNL's own AGENTS.md (`build/_deps/tnl-src/AGENTS.md`, error-handling section) prescribes:
`TNL_ASSERT_*` macros for debug assertions, `TNL_BACKEND_SAFE_CALL` for CUDA/HIP API calls, and
`TNL::Exceptions::BackendRuntimeError` for runtime errors. TNL-SPH currently uses **zero** `TNL_ASSERT_*`
and instead mixes four idioms: bare `assert()` (inert in Release), `throw std::runtime_error`,
`std::cerr` diagnostics, and `std::exit`. Verified inventory:

bare `assert()` (4 sites, all overlap-state guards in MR buffers):
- include/SPH/MultiresolutionRectangleBufferLocalTimeStepping.h:125,141
- include/SPH/MultiresolutionRectangleBuffer.h:160,176

`throw` sites (7):
- include/SPH/parseConfigFile.h:114,195,202,225 (`std::runtime_error`)
- include/SPH/InitReport.hpp:73 (`std::runtime_error`)
- include/SPH/solvers/SolverMultiSetBlockMultiresolution.hpp:344 (`std::runtime_error`)
- include/SPH/DecompositionTopology.h:125 (`std::out_of_range` — acceptable as-is or align)

`std::cerr` per file (36 hits / 16 files):
- include/SPH/configSetup.h:9, configInit.h:4, SimulationMonitor.h:3, shared/debugging.h:3,
  SPHMultiSet_CFD.hpp:2 (was SPHMultiset_CFD.hpp), legacy/SolverMultiSet.hpp:2 (post-quarantine path),
  shared/PeriodicBoundaryConditions.h:2, Models/WCSPH_DBC/OpenBoundaryConfig.h:2,
  Models/WCSPH_BI/OpenBoundaryConfig.h:2, solvers/SolverMultiSetRemeshed.hpp:1,
  solvers/SolverMultiSet.hpp:1, solvers/SolverMultiSetBlockMultiresolutionLocalTimestepping.hpp:1
  (post-T9 path ends with LocalTimeStepping), solvers/SolverMultiSetBlockMultiresolution.hpp:1,
  solvers/SolverMultiSetBase.hpp:1, Models/WCSPH_DBC/OpenBoundaryConditions.hpp:1,
  Models/WCSPH_BI/OpenBoundaryConditions.hpp:1

`std::exit` (2): include/SPH/configInit.h:52,59 — after `--config-help` printing; CLI semantics, keep
behavior but route output through `TNL::Logger`.

### Proposed approach (former plan task T15)

1. Replace the 4 bare `assert(` with the matching `TNL_ASSERT_*` macro + message (TNL/Assert.h).
2. Triage the 36 `std::cerr` sites: genuine user-facing diagnostics → keep (documented keep-list in the
   change's evidence) or route through `TNL::Logger`; true fatal errors → `TNL::Exceptions::ConfigError`
   / `BackendRuntimeError`.
3. Type the 6 `std::runtime_error` throws from `TNL/Exceptions/*` (ConfigError for parseConfigFile.h +
   InitReport.hpp; BackendRuntimeError for the MR solver check).
4. Keep `DecompositionTopology.h:125` `std::out_of_range` or align; document the choice.
5. Keep `configInit.h` exit-after-help semantics; logger for output. Do NOT touch the known
   uncaught-ConfigError abort on bare `--config-help` (pre-existing issue listed in AGENTS.md).
6. Guardrails: no exceptions inside hot GPU/device code paths; no silent return-contract changes.

Acceptance when executed: `TNL_ASSERT_*` usage ≥4, bare `assert(` == 0, `std::cerr` count ≤ the
documented keep-list, builds + unit tier green.

---

## DEFERRED-T17/T18 — TNL-style Doxygen floor (both batches)

Deferred 2026-10-22 at user request during the plan walkthrough ("documentation resolved in a separate
session"). Two micro-fixes that were embedded in T18 were RELOCATED to plan task T14 (wrong RSPH CLI
delimiter label + stale `} // BCTypes` closing comments) so nothing is lost by this deferral.

### Finding

Header documentation is essentially absent (verified census): 0 of 116 headers carry `\file`/`@file`,
0 `\class`, `\brief` appears in only 30 of 116 files, `\param` ×17 and `\tparam` ×9 library-wide —
while TNL heads every file with a full Doxygen block (`\brief`, `\tparam`, `\ref`, `\par`;
idiom: build/_deps/tnl-src/src/TNL/Containers/Array.h:14-25). Examples of critical undocumented API
files: include/SPH/ParticleSet.h (base class of all particle sets — zero doc blocks),
VariablesBase.h, SPHTraits.h, SPHMultiset_CFD.h driver API, Models/*/Interactions.h model classes.

### Proposed approach (former plan tasks T17 + T18)

1. Every header gets a `\file <basename>` + `\brief` file block in TNL style.
2. Every public class gets a `\class`/`\brief` block (+ `\tparam` on templates, `\param` on the
   public methods already partially covered). Describe existing behavior — do not invent physics.
3. Batch split (for independent execution): batch A = top level + solvers/ + shared/ (all headers
   except Models/; exclude legacy/ — quarantine-covered), batch B = include/SPH/Models/** incl.
   policy files (DiffusiveTerms, ViscousTerms(post-T6), EquationOfState, PressureGradient,
   BoundaryViscousTerms, RiemannSolvers) + per-model BoundaryConditions*/IntegrationSchemes/
   OpenBoundary*; module-level page documenting the compile-time policy-bundle philosophy.
4. Acceptance when executed: every header in batch scope carries `\file`; doc-only commits keep
   builds byte-identical.

---

## DEFERRED-T19 — Extract shared open-boundary implementation (BI/DBC twins)

Deferred 2026-10-22 at user request during the plan walkthrough.

### Finding

Two model files are a maintenance copy-paste (verified): the open-boundary data-extrapolation headers
include/SPH/Models/WCSPH_DBC/OpenBoundaryConditionsDataExtrapolation.h and
include/SPH/Models/WCSPH_BI/OpenBoundaryConditionsDataExtrapolation.h are **707 lines each with ~74
differing lines (~92% identical)**; their OpenBoundaryConditions.hpp twins differ by only 13 lines.
Additionally the BI file reaches across model boundaries:
`#include "../WCSPH_DBC/details.h"` (its line 2) to obtain `matrixCorrection2D`, which lives as a
free function directly in `TNL::SPH` inside WCSPH_DBC/details.h.

### Proposed approach (former plan task T19)

1. Hoist shared code into one canonical template `Models/OpenBoundaryConditionsDataExtrapolation.{h,hpp}`
   (+ shared conditions base); expose the ~10% model-specific behavior as CRTP/policy hooks consumed by
   each model's thin per-model file; same treatment for OpenBoundaryConditions.hpp (13-line diff).
2. Move `matrixCorrection2D` into a shared `include/SPH/Models/Details.h` in `namespace TNL::SPH::details`
   (mirroring SHTC/Details.h's pattern — interplay with executor task T20b, which only adds
   `namespace details` inside each per-model Details.h and remains executable independently).
3. Acceptance when executed: exactly one shared implementation; each model-specific remainder < 60 diff
   lines; `grep -rni 'WCSPH_DBC/details.h' include/SPH/Models/WCSPH_BI` == 0; DBC+BI example builds green;
   one ≤10-step open-channel boot of each matches pre-change state.

NOTE: RSPH/SHTC keep using the EmptyOpenBCPlaceholder path (post-T11), unchanged by this task.

---

## DEFERRED-T20 — Normalize the model-directory skeleton

Deferred 2026-10-22 at user request during the plan walkthrough.

### Finding

Of the 4 production models, only 2 match the documented `include/SPH/Models/<name>/` layout:
- SHTC names its entry `SHTC.{h,hpp}` instead of `Interactions.{h,hpp}` (class `SHTC` inside) and adds a model-specific `stress.h`;
- `details.h` semantics differ: WCSPH_DBC/details.h puts `matrixCorrection2D` free into `TNL::SPH`, while SHTC/Details.h uses `namespace details`;
- scheme-variable classes follow 3 naming patterns: DBC/RSPH generic `IntegrationSchemeVariables`,
  BI inside one directory mixes `MidpointIntegrationSchemeVariables` + `RK4IntegrationSchemeVariables`
  with `SymplecticVerletSchemeVariables`, SHTC uses `ExplicitEulerSchemeVariables`;
- BC extras: BI keeps model-specific `schemeTypes.h` (lowercase-start, no peer);
- per-model copies of `IntegrationSchemes/VerletScheme.h` exist in DBC, BI, RSPH (drift unrendezvoused).

### Proposed approach (former plan task T20)

1. `git mv include/SPH/Models/SHTC/SHTC.{h,hpp}` → `Interactions.{h,hpp}` (class untouched), fix bottom-pair includes + SHTC example config includes.
2. `namespace details` inside every per-model `Details.h` (mirror SHTC/Details.h).
3. Unify scheme-variable class names to `<Scheme>SchemeVariables` (`VerletSchemeVariables`, `MidpointSchemeVariables`, `RK4SchemeVariables`); keep the per-model `using IntegrationSchemeVariables = ...` alias consumed by SPHMultiSet_CFD.h's drivers.
4. Diff the three `VerletScheme.h`; hoist alias-level commonalities into `Models/IntegrationSchemes/VerletSchemeBase.h`, else document divergence.
5. Write the sanctioned skeleton to `Documentation/model-skeleton.md` (+ AGENTS.md pointer); keep BI's BC struct names (`BI_*`, `BIConsistent_*` = user config tokens) untouched.
6. Acceptance when executed: all 4 model dirs match the documented skeleton; all representative builds green (incl. RSPH + lidDrivenCavity2D_SHTC_cuda); unit tier baseline-consistent.

---

## DEFERRED-T21 — Multiresolution next-audit bundle (buffers, LTS internals, MR marker backlog)

Deferred 2026-10-22 at user request during the plan walkthrough: "skip everything related to
multiresolution algorithms — handled separately during the next audit." Deliberately NOT part of this
bundle (kept in the active plan): T13 legacy-quarantine move and the W2 MR-adjacent renames
(T9 solver spelling, T11 MR config-helper rename) — user confirmed keeping all of them.

### Findings

- Three MR buffer generations coexisted pre-quarantine: DEAD `include/SPH/MultiresolutionBuffer.h`
  (770 L, legacy-moved by T13) plus the two LIVE monoliths `MultiresolutionRectangleBuffer.h`
  (1308 L) and `MultiresolutionRectangleBufferLocalTimeStepping.h` (1433 L, LTS deltas interleaved;
  defines `MultiresolutionBoundaryLTS`). None use the project's .h/.hpp split idiom.
- MR-marker backlog carried inside those files (verified counts): MultiresolutionRectangleBuffer.h 34,
  MultiresolutionRectangleBufferLocalTimeStepping.h 33, dead MultiresolutionBuffer.h 29
  (the dead one's markers leave with T13's quarantine), plus the legacy solver's 95 (quarantined).
  Plan task T24 skips these clusters — the bundle owns their triage.
- The fixtures' pre-existing `dummy2D-lts-center` diff (AGENTS.md expected-on-branch list) lives in
  this feature area — reference it before any MR behavior-affecting rework.

### Proposed approach (former plan task T21, expanded)

1. Read Documentation/multiresolution-reference.md, multiresolution-ghost-boundaries-design.md and
   multiresolution-local-timestepping.md FIRST (project rule).
2. Merge the two live generations into a decl/impl-split family: shared MultiresolutionRectangleBuffer.{h,hpp}
   + ...LocalTimeStepping.{h,hpp} carrying only LTS deltas; same mechanical .h/.hpp split for
   ParticleSet.h (469 L, no .hpp) IF the bundle session chooses to include it — ParticleSet is not
   MR-internal, it was bundled here only as a "monolith" peer.
3. Triage the MR marker backlog listed above (fix-mechanical vs park vs known-broken label) in the
   same session while the files are already being reworked.
4. Acceptance when executed: LTS-only delta isolated (record before/after line counts); MR + MR-LTS
   example builds green; init-level dummy-MR fixture run shows no new diffs vs the pre-work capture
   (incl. the pre-existing dummy2D-lts-center diff).

---

## DEFERRED-T22 — Solver-stack audit (SPHMultiSet_CFD decomposition etc.)

Deferred 2026-10-22 at user request during the plan walkthrough: "we slowly moved from
SPHMultiset_CFD to solvers — we audit the solver in detail next time." The driver/solver
architecture is the next audit's own session; only its comment-level noise stays in the
active plan (T24 triages the 12 markers inside SPHMultiSet_CFD.hpp without touching architecture).

### Findings

- The top-level driver `include/SPH/SPHMultiset_CFD.{h,hpp}` (→ `SPHMultiSet_CFD` post-T8;
  376 + 1039 lines) does ~5 jobs in one class: init (8 methods), config (SIX Config/ParameterContainer
  pairs), the time-loop, I/O, and monitoring; ~40 public methods; self-flagged
  comments: `IntegratorPointer integrator; // I hate this.` (current line 314) and
  `// TEMP: And btw the names are AWFUL` (:338).
- The broken never-instantiated `operator<<` (:364-373, calls writeProlog(Logger&) vs
  writeProlog(bool)) is REMOVED by in-plan task T14 — the split session need not care about it.
- The migration question (driver vs solvers/ hierarchy) is exactly what this audit session
  must decide — gather: who consumes SPHMultiSet_CFD today vs who already instantiates
  SolverMultiSet* directly (grep: examples all call SPHMultiset_CFD<Model> via case.h;
  fixture configs use the two solver generations directly).

### Proposed approach (former plan task T22)

1. Keep the public driver contract byte-stable for examples during any refactor
   (init/performNeighborSearch/interact/computeTimeStep/integrate*/makeSnapshot/measure/updateTime).
2. Candidate splits (evaluate during the audit): SimulationInitializer (init/read) +
   SimulationConfigState (the 6 config pairs) + SimulationMonitor delegation; .hpp split
   driver/IO; resolve the 2 self-flagged comments by fixing the names themselves.
3. Decide the SPHMultiSet_CFD ↔ solvers/ end-state (the documented migration).
4. Acceptance when executed: no case.h edits; builds green; .hpp < 500 lines; one-config struct;
   unit tier + one DBC boot unchanged.

---

## DEFERRED-T24 — TODO/FIXME marker backlog triage (owner: user, personally)

Deferred 2026-10-22 at user request during the plan walkthrough: "keep them, I want to resolve it
myself." The active plan does NOT touch the marker backlog. This registry entry just preserves the
verified inventory + the suggested triage rule, so the user's own session can start from data.

### Verified inventory (~322 TODO/FIXME lines; exact counts audited)

- Quarantined with T13's legacy move (leave or resolve on that branch): legacy SolverMultiSet.hpp ~95,
  dead MultiresolutionBuffer.h 29.
- DEFERRED-T21 bundle (MR markers — resolve during the next MR audit):
  MultiresolutionRectangleBuffer.h 34, MultiresolutionRectangleBufferLocalTimeStepping.h 33.
- Live, in-plan scope (WAS T24's): SPHMultiSet_CFD.hpp 12 (comment-level only; architecture is
  DEFERRED-T22), solvers/SolverMultiSetBase.hpp 8 (incl. SIX `//FIXME: This function is broken`
  methods, lines 167-239), ParticleSet.h 8, SHTC/Control.h 7, OBC twins 6+6 (still duplicated —
  DEFERRED-T19), WCSPH_BI/Interactions.hpp 6, WCSPH_BI/IntegrationSchemes/SymplecticVerletScheme.h 6;
  plus `//temp` ×22, `is broken` ×7, `NOTE:` ×11.

### Suggested triage rule (former plan task T24)

For each marker: (a) fix if mechanical (typo, completed TODO, stale comment); (b) park as one line in
`Documentation/refactor-followups.md` (this registry's repo copy) if physics judgment is needed;
(c) label `FIXME(known-broken)` explicitly for the 6 behavioral gaps — do not silently "fix" them.

TARGET end-state when executed: `grep -rn 'FIXME\|TODO' include/SPH` count equals the parked-list
length; every parked item has one registry line; builds green.

---

<!-- Future triage items will be parked below this line. -->

## FINDING (from T10 execution, 2026-11-05): SensorInterpolation invalid write with zero-sized monitor outputs
- Booting a case with finalTime below the monitor output periods triggers CUDA error 700: invalid write in SensorInterpolation::interpolate (zero-sized monitor arrays). Pre-existing, NOT introduced by the refactor (verified by T10 worker via compute-sanitizer; worked around in the external test config only).
- Repro: run any example with finalTime < min(monitor period) and monitors enabled. Candidate fix belongs to the measuretool/monitor owners — park for a dedicated session; do not bundle into style refactors.

## [T16] draft/broken marker census

Executed 2026-11-05 (T16, comments-only; user decision: drafts kept as documented drafts).
Command: `grep -rni "not implement\|is broken" include/SPH` — 9 sites before, 7 sites after T16 edits.

Resolved by T16 (converted to explicit `\todo` DRAFT Doxygen notes, zero behavior change):
- include/SPH/Models/WCSPH_DBC/BoundaryConditionsTypes.h:27 — MDBCVelocity (was `//not implement yet`)
- include/SPH/Models/WCSPH_DBC/BoundaryConditionsTypes.h:36 — GWBC (was `//not implemented`)

Remaining sites (listed-but-NOT-edited; user-owned backlog per DEFERRED-T24 triage rule):
- include/SPH/solvers/SolverMultiSetBase.hpp:176 — `//FIXME: This function is broken`
- include/SPH/solvers/SolverMultiSetBase.hpp:190 — `//FIXME: This function is broken`
- include/SPH/solvers/SolverMultiSetBase.hpp:207 — `//FIXME: This function is broken`
- include/SPH/solvers/SolverMultiSetBase.hpp:223 — `//FIXME: This function is broken`
- include/SPH/solvers/SolverMultiSetBase.hpp:237 — `//FIXME: This function is broken`
- include/SPH/solvers/SolverMultiSetBase.hpp:252 — `//FIXME: This function is broken`
  (the six known-broken methods cited under DEFERRED-T24; label FIXME(known-broken) at triage, do not silently fix)
- include/SPH/Models/WCSPH_DBC/Interactions.h:130 — `* TODO: Not implemented.` (doc TODO on an unimplemented detail; T24 triage)

Census note: the new T16 draft-stub comments use the wording "unimplemented", so they are deliberately absent from the census grep above.
