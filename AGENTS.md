# TNL-SPH PROJECT KNOWLEDGE BASE

**Updated:** 2026-10-01
**Branch:** multiresolution-patch-lts

## OVERVIEW

TNL-SPH is a C++20/CUDA header-only Smoothed Particle Hydrodynamics (SPH) framework built on top of the Template Numerical Library (TNL).
It provides pluggable SPH formulations (δ-WCSPH, WCSPH with boundary integrals, RSPH, SHTC), interchangeable diffusive/viscous/EOS terms,
open (inlet/outlet) and periodic boundary conditions, block-based multiresolution particle refinement with local time stepping,
and distributed multi-GPU execution through CUDA-aware MPI with 1D domain decomposition.

## STRUCTURE

```
.
├── include/
│   └── SPH/                            # Core header-only library
│       ├── Models/                     # Pluggable SPH formulations
│       │   ├── WCSPH_DBC/              # δ-WCSPH with dynamic boundary conditions
│       │   ├── WCSPH_BI/               # WCSPH with boundary integrals (has experimental/ subdir)
│       │   ├── RSPH/                   # Riemann SPH
│       │   ├── SHTC/                   # SHTC formulation
│       │   ├── WCSPH_MFD/              # archived/incomplete (only Interactions.hpp)
│       │   ├── SPHMultisetSolverTemplate.h   # base class for models
│       │   └── DiffusiveTerms.h, VisousTerms.h, EquationOfState.h,
│       │       DensityFilters.h, VelocityFilters.h, RiemannSolvers.h,
│       │       PressureGradient.h, BoundaryViscousTerms.h   # interchangeable policy terms
│       ├── solvers/                    # SolverMultiSet* family (single-res, remeshed, MR, MR+LTS)
│       ├── shared/                     # Cross-model utilities (Measuretool, PST, interpolation, energy, ...)
│       ├── SPHMultiset_CFD.h/.hpp      # Main simulation orchestrator
│       ├── ParticleSet.h, Fluid.h, Boundary.h   # particle-set base + fluid/boundary sets
│       ├── OpenBoundaryBuffers.h, OpenBoundaryConfig.h, PeriodicBoundaryBuffers.h
│       ├── MultiresolutionBuffer.h, MultiresolutionRectangleBuffer.h,
│       │   MultiresolutionRectangleBufferLocalTimeStepping.h
│       ├── GhostZone.h, BoundaryGhostParticles.h   # multiresolution ghost handling
│       ├── DecompositionTopology.h     # 1D MPI domain decomposition
│       ├── SPHTraits.h, Kernels.h, TimeStep.h, SimulationMonitor.h, TimeMeasurement.h
│       └── configInit.h, configSetup.h, parseConfigFile.h   # JSONC runtime configuration
├── examples/                           # Stand-alone case problems, one directory per case
│   ├── WCSPH-DBC/    WCSPH-BI/    RSPH/    SHTC/
│   └── extra/  resources/              # auxiliary cases + shared images/data
├── src/
│   ├── tools/                          # Python utilities (multirunner.py, splitToSubdomains.py,
│   │                                   #  testSolver.py, plotting, VTK writers) — imported by init.py scripts
│   ├── UnitTests/                      # CTest-registered test executables
│   │   ├── multiresolution/            # dummy MR simulations used as test fixtures
│   │   └── DistributedParticles/
│   ├── Benchmarks/                     # optional (TNL-SPH_BUILD_BENCHMARKS)
│   └── extra/  archivedParts/          # auxiliary and legacy code
├── Documentation/                      # Doxygen config + multiresolution design documents
├── CMakeLists.txt                      # Root config: C++20+CUDA, -Werror, FetchContent deps
└── .gitlab-ci.yml                      # Release + Debug builds (Ninja, CMAKE_BUILD_PARALLEL_LEVEL=4)
```

Example layout convention — each example directory contains:

```
<case>.cu                  # solver entry, compiled to <case>_cuda
case.h                     # thin main() driver including template/config.h
template/                  # config.h (compile-time policy bundle), config_template.jsonc, ...
init.py                    # generates particles/config into sources/ from template/
run.py                     # init (if needed) + run driver; resolves binary in build/ mirror path
sources/                   # GENERATED: particles (.vtk) + config.jsonc — do not edit by hand
results/                   # GENERATED: simulation output — never commit
postpro.py                 # post-processing
```

## WHERE TO LOOK

| Task | Location | Notes |
|------|----------|-------|
| Add a new SPH formulation | `include/SPH/Models/<NewModel>/` | Mirror `WCSPH_DBC/` layout: `control.h` (config class), `Interactions.h/.hpp` (model class deriving from `SPHMultisetSolverTemplate`), `Variables.h`, `BoundaryConditionsTypes.h`, `IntegrationSchemes/` |
| Add a diffusive term | `include/SPH/Models/DiffusiveTerms.h` | Existing: `None`, `MolteniDiffusiveTerm`, `FourtakasDiffusiveTerm`; select in the example's `SPHParams` |
| Add a viscous term | `include/SPH/Models/VisousTerms.h` (fluid), `BoundaryViscousTerms.h` (boundary) | Existing: `ArtificialViscosity`, `PhysicalViscosity`, `PhysicalViscosity_MVT/_MGVT`, `CombinedViscosity` |
| Add an equation of state | `include/SPH/Models/EquationOfState.h` | Existing: `TaitWeaklyCompressibleEOS` (+ linearized) |
| Add a pressure gradient formulation | `include/SPH/Models/PressureGradient.h` | Existing: `Symmetric`, `TIC` |
| Add a boundary condition type (DBC) | `include/SPH/Models/WCSPH_DBC/BoundaryConditionsTypes.h` + `BoundaryConditions/` | Enum-like structs (`DBC`, `MDBC`, `OpenBoundary`), each with `integrateInTime()`; `MDBCVelocity`/`GWBC` are stubs |
| Add a boundary condition type (BI) | `include/SPH/Models/WCSPH_BI/BoundaryConditionsTypes.h` | Numerical variants: `BI_*`, `BIConsistent_*`, `BIConservative_numeric` |
| Add an integration scheme | `include/SPH/Models/<model>/IntegrationSchemes/` | Verlet, SymplecticVerlet, RK4, Midpoint (BI), ExplicitEuler (SHTC); each file defines scheme + scheme-variables classes |
| Change time stepping | `include/SPH/TimeStep.h` | `ConstantTimeStep`, `VariableTimeStep`, `VariableTimeStepWithReduction`; selected via `SPHParams::TimeStepping` |
| Change solver strategy (MR / remeshing / LTS) | `include/SPH/solvers/SolverMultiSet*.h/.hpp` | `SolverMultiSet` (single-resolution), `SolverMultiSetRemeshed`, `SolverMultiSetBlockMultiresolution`, `SolverMultiSetBlockMultiresolutionLocalTimestepping` |
| Multiresolution buffers / level overlays | `include/SPH/MultiresolutionRectangleBuffer{,LocalTimeStepping}.h`, `GhostZone.h`, `BoundaryGhostParticles.h` | Design docs: `Documentation/multiresolution-*.md` — read before touching |
| Open boundaries (inlet/outlet) | `include/SPH/OpenBoundary*.h` + per-model `OpenBoundaryConditions*` / `OpenBoundaryConfig.h` | Buffers store free-zone/buffer-zone particles; per-model extrapolation in `OpenBoundaryConditionsDataExtrapolation.h` |
| Periodic boundaries | `include/SPH/PeriodicBoundaryBuffers.h`, `include/SPH/shared/PeriodicBoundaryConditions.h` | `SPHConfig::numberOfPeriodicBuffers` selects copies per side |
| MPI decomposition | `include/SPH/DecompositionTopology.h`, `src/tools/splitToSubdomains.py` | 1D split with domain overlaps; `SubdomainDescriptor`/`InterfaceDescriptor` |
| Measurement: probes, sensors, grid output | `include/SPH/shared/Measuretool.h` (`InterpolateToGrid`, `SensorInterpolation`, `SensorWaterLevel`), `measureVolumetricFlowRate.h`, `massConservationMonitor.h` | Driven by `SimulationMonitor.h` timers; runtime-configured via `.jsonc` |
| Create a new example | copy e.g. `examples/WCSPH-DBC/damBreak2D_WCSPH-DBC/` | Keep the layout convention; CMakeLists adds a `<case>_cuda` executable linked to `TNL::TNL` |
| Unit tests | `src/UnitTests/{multiresolution,DistributedParticles}/` | Plain test executables registered via `add_test`; MR fixtures are dummy simulations with their own `run.py` |

## CODE MAP

| Symbol | Type | Location | Role |
|--------|------|----------|------|
| `SPHMultiset_CFD<Model>` | class | `include/SPH/SPHMultiset_CFD.h/.hpp` | Top-level simulation: owns `Fluid`, `Boundary`, `OpenBoundary` particle sets, time stepping, VTK reader/writer, `SimulationMonitor`; API used by `case.h` (`init`, `performNeighborSearch`, `interact`, `computeTimeStep`, `integrateVerletStep`, `makeSnapshot`, `measure`, `updateTime`) |
| `SolverMultiSetBase<Model>` | class | `include/SPH/solvers/SolverMultiSetBase.h/.hpp` | Base of the modular solver hierarchy |
| `SolverMultiSet` / `SolverMultiSetRemeshed` | class | `include/SPH/solvers/SolverMultiSet{,Remeshed}.h/.hpp` | Single-resolution and remeshed solver variants |
| `SolverMultiSetBlockMultiresolution` | class | `include/SPH/solvers/SolverMultiSetBlockMultiresolution.h/.hpp` | Block-based multiresolution solver (level overlays, coarse→fine interaction) |
| `SolverMultiSetBlockMultiresolutionLocalTimestepping` | class | `include/SPH/solvers/SolverMultiSetBlockMultiresolutionLocalTimestepping.h/.hpp` | MR solver with per-level local time stepping — focus of the current branch |
| `ParticleSet<...>` | class | `include/SPH/ParticleSet.h` | Base class for every particle collection (fluid, boundary, open-boundary, MR buffers) |
| `Fluid` / `Boundary` / `OpenBoundary` | class | `include/SPH/{Fluid,Boundary,OpenBoundaryBuffers}.h` | Particle sets for fluid / boundary / open-boundary zones |
| `VariablesBase<Derived>` | class | `include/SPH/VariablesBase.h` | CRTP base for per-model variable sets (pairs with `VariablesGhostInit.h`) |
| `SPHFluidTraits<SPHConfig>` | class | `include/SPH/SPHTraits.h` | Type-traits bundle (index, real, vector, array types) |
| `WCSPH_DBC` / `WCSPH_BI` / `RSPH` / `SHTC` | class | `include/SPH/Models/*/Interactions.h(.hpp)`, `SHTC.h` | SPH model implementations, composed of the policy terms below |
| `SPHMultisetSolverTemplate<Particles,ModelConfig>` | class | `include/SPH/Models/SPHMultisetSolverTemplate.h` | Base class shared by models (RSPH, SHTC derive publicly) |
| `WendlandKernel` | class | `include/SPH/Kernels.h` | Wendland C2 kernel, 2D/3D specializations |
| `DiffusiveTerms::*`, `ViscousTerms::*`, `EquationsOfState::*` | classes | `include/SPH/Models/{DiffusiveTerms,VisousTerms,EquationOfState}.h` | Interchangeable compile-time policy terms |
| `ConstantTimeStep` / `VariableTimeStep` / `VariableTimeStepWithReduction` | class | `include/SPH/TimeStep.h` | Time-step computation; variable variants from CFL-type conditions |
| `MultiresolutionBoundary` / `MassNodes` | class | `include/SPH/MultiresolutionRectangleBuffer{,LocalTimeStepping}.h` | Per-level boundary particle buffers for MR; `MultiresolutionBoundaryLTS` adds local time stepping |
| `BoundaryGhostParticles` | class | `include/SPH/BoundaryGhostParticles.h` | Boundary ghost values on MR block interfaces (updated from the neighbor fluid — model-dependent) |
| `ParticleZone` | class | `include/SPH/GhostZone.h` | Ghost-zone classification of MR blocks |
| `DecompositionTopology(Flat)` | class | `include/SPH/DecompositionTopology.h` | 1D MPI subdomain layout and 1-ring interface graph |
| `ParticleSystemConfig` / `SPHConfig` / `SPHParams` (“SPHDefs”) | struct | per-example `template/config.h` | Compile-time policy bundle: particle-system types, kernel, EOS, diffusive/viscous terms, BC type, time stepping, integration scheme |
| `SimulationMonitor` | class | `include/SPH/SimulationMonitor.h` | Timed snapshots, measurements, energy evaluation, mass conservation |
| `InterpolateToGrid` / `SensorInterpolation` / `SensorWaterLevel` | class | `include/SPH/shared/Measuretool.h` | Grid interpolation and point-sensor measurements used by the monitor |
| `MeasureVolumetricFlowRate` / `MassConservationMonitor` | class | `include/SPH/shared/measureVolumetricFlowRate.h`, `massConservationMonitor.h` | Open-boundary flow-rate and mass tracking |
| `TimerMeasurement` | class | `include/SPH/TimeMeasurement.h` | Per-stage wall-time accounting |
| `ElasticBounce`, `PST`, `MFD`, `TaylorMonomials` | class | `include/SPH/shared/` | Cross-model helpers: boundary-bounce modes, particle shifting, MLS/Taylor interpolation |

## CONVENTIONS

- **Commit messages**: short sentence-case, imperative mood, no trailing period; WIP may be prefixed with `Draft - ` (see `git log -10`).
- **Assisted-by trailer**: When AI tools contribute, add `Assisted-by: AGENT:MODEL` (e.g. `Assisted-by: Opencode:kimi-k3`).
- **Header-only library**: `TNL_SPH` is a CMake `INTERFACE` target (alias `TNL::TNL_SPH`); all compilation cost sits in the executables.
- **C++ source suffixes**: headers `.h` hold declarations; `.hpp` files are template implementations included from the corresponding `.h`.
- **Formatting** (`.editorconfig`): 3 spaces for C++/CUDA and CMake, 4 spaces for Python; LF line endings; final newline; no trailing whitespace.
- **TNL idioms** (mirror the parent library): `__cuda_callable__` instead of separate `__host__`/`__device__`, `TNL_ASSERT_*` macros, `#pragma once` include guards.
- **Git discipline**: never `git add -A`/`git add .` — stage files explicitly (`results/`, `sources/`, VTK, logs must not be committed).
- **-Werror builds**: warnings are compile errors in both C++ and CUDA host code (few `-Wno-error=` escapes in root `CMakeLists.txt`); fix warnings, never suppress them locally.
- **Executable naming**: every example target is suffixed `_cuda` (`<case>_cuda`) and links `TNL::TNL`; OpenMP is enabled per-target with `-DHAVE_OPENMP` (for nvcc: `-Xcompiler=-fopenmp`).
- **Runtime configuration**: JSON-with-comments `.jsonc` files (comments stripped by `parseConfigFile`); numerical/runtime knobs belong there, compile-time model policy belongs in `template/config.h`.
- **Dependencies** are fetched with `FetchContent`: TNL from the fork `gitlab.com/tomashalada/tnl`, branch `TH/particles` (NOT upstream `tnl-project/tnl`), and nlohmann/json pinned to v3.11.3.
- **MPI is optional** (`TNL_SPH_ENABLE_MPI`, default ON) but the library targets CUDA; `Device` is fixed to `TNL::Devices::Cuda` in configs.
- **CUDA architecture**: defaults to `"native"`; use `75` when configuring on a GPU-less machine (same fallback as CI).

## ANTI-PATTERNS (THIS PROJECT)

- **Parallel builds above 2 jobs**: C++20/CUDA template instantiation is memory-hungry; highly parallel `make -jN`/`ninja -jN` builds get `cicc`/`cc1plus` SIGKILLed by the OOM killer. **Always build with two cores** (`--parallel 2` / `-j2`). CI itself caps at 4.
- **Building all targets for a small change**: rebuild only the affected case with `--target <case>_cuda`; full builds waste a lot of time.
- **Editing `sources/` by hand**: `sources/` is regenerated by `init.py` from `template/` (and `run.py` auto-runs init when it is missing). Make lasting changes in `template/config.h` / `template/config_template.jsonc` and re-run `init.py`.
- **Editing `build/_deps/tnl-src/`**: that tree is a FetchContent mirror of the forked TNL; changes there are in a separate repository (with its own `AGENTS.md`) and belong to the TNL project, not here.
- **Confusing the two `SolverMultiSet` headers**: `include/SPH/SolverMultiSet.h/.hpp` (top level, unitary legacy-style solver) and `include/SPH/solvers/SolverMultiSet.h/.hpp` (part of the `SolverMultiSetBase` hierarchy) are different classes; check includes before editing.
- **Including both MR buffer implementations**: `MultiresolutionBuffer.h` and `MultiresolutionRectangleBuffer.h` define the same class names (`MassNodes`, `MultiresolutionBoundary`) — never include both in one translation unit; the LTS variant lives in `MultiresolutionRectangleBufferLocalTimeStepping.h`.
- **Committing generated content**: `results/`, VTK/VTU output, `sources/`, `__pycache__`, and simulation logs are git-ignored helpers; adding them with `-A`/force pollutes the repo.
- **Host-only builds**: TNL's `HeadersOnly` mode and non-CUDA devices are unsupported here; CUDA is required and configs pin `TNL::Devices::Cuda`.
- **CUDA error-chasing bottom-up**: `nvcc` errors cascade; fix the first error, rebuild with `--parallel 2` (or the single failing TU), repeat. Runtime CUDA failures: `compute-sanitizer`; segfaults: rebuild `RelWithDebInfo`/`Debug`, then `gdb`.

## UNIQUE STYLES

- **Compile-time policy bundle per simulation**: every case defines `ParticleSystemConfig`, `SPHConfig`, and `SPHParams` (aliased to `SPHDefs`) in `template/config.h` — kernel, EOS, diffusive/viscous terms, BC type, time stepping, and integration scheme are all `using` typedefs swapped per case; `Model = <Formulation><ParticlesSys, SPHParams>` and `Simulation = SPHMultiset_CFD<Model>` complete the assembly.
- **Uniform `case.h` driver loop**: every example's `case.h` is the same thin `main()`: `init → while(runTheSimulation()){ performNeighborSearch; interact; computeTimeStep; integrate{Scheme}Step(BCType::integrateInTime()); makeSnapshot; measure; updateTime; }` — model differences live in the headers, not the loop.
- **Everything is a `ParticleSet`**: `Fluid`, `Boundary`, `OpenBoundary`, and the multiresolution buffers all derive `ParticleSet` with their own variable sets, so ghost/boundary logic shares the particle-machinery API.
- **CRTP variable sets**: `VariablesBase<Derived>` + per-model `FluidVariables`/`BoundaryVariables`/`OpenBoundaryVariables` (two layers: `FluidVariablesBase` + final class); ghost initialization via `VariablesGhostInit.h`.
- **Model subpackages mirror each other**: each `Models/<name>/` carries `control.h` (typed config class), `Interactions.h/.hpp`, `Variables.h`, `BoundaryConditionsTypes.h`, `IntegrationSchemes/`, and (where relevant) `OpenBoundary*` files — adding a model means reproducing this layout.
- **Python-driven case management**: `init.py` generates particles + config from `template/` into `sources/` (imports helpers from `src/tools`), `run.py` resolves the binary in the build mirror `build/<same relative path>/<case>_cuda` and streams solver output; distributed examples run under `mpirun` after `splitToSubdomains.py`.
- **Measurement as monitor plug-ins**: grid interpolation, point sensors, water level, mass conservation, and flow rate are assembled in `SimulationMonitor` from `shared/Measuretool*.h` and timer-driven from the runtime `.jsonc` config.
- **Suffix conventions**: executables `_cuda`; config templates `config_template.{h,jsonc}` generated into concrete `sources/config.jsonc`.

## COMMANDS

```bash
# Hard rule: never exceed 2 parallel build jobs (template-heavy → OOM kills cicc/cc1plus)
export CMAKE_BUILD_PARALLEL_LEVEL=2

# Configure (Release / Debug / RelWithDebInfo; use CMAKE_CUDA_ARCHITECTURES=75 on GPU-less machines)
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release

# Build everything — or much faster, only what you touch
cmake --build build --parallel 2
cmake --build build --parallel 2 --target damBreak2D_WCSPH-DBC_cuda

# Initialize + run an example (init runs automatically when sources/ is missing)
./examples/WCSPH-DBC/damBreak2D_WCSPH-DBC/run.py
./examples/WCSPH-DBC/damBreak2D_WCSPH-DBC/init.py            # force re-generation of particles
./examples/WCSPH-DBC/damBreak2D_WCSPH-DBC/run.py --config sources/custom.jsonc

# Distributed examples (1D domain split first, then MPI)
python3 src/tools/splitToSubdomains.py --help
mpirun -np 2 ./build/examples/WCSPH-DBC/damBreak2D_WCSPH-DBC_distributed/damBreak2D_WCSPH-DBC_distributed_cuda --config sources/config.jsonc

# Tests
ctest --test-dir build
ctest --test-dir build -R DistributedParticleSystemTest

# LSP: build/compile_commands.json is produced at configure time
```

## MULTIRESOLUTION & LOCAL TIME STEPPING

- Design rationale before coding: `Documentation/multiresolution-reference.md`,
  `multiresolution-ghost-boundaries-design.md`, and `multiresolution-local-timestepping.md`
  are the of-record documents — consult them before touching MR buffers, ghosts, or LTS.
- Solvers: `SolverMultiSetBlockMultiresolution` (uniform step across levels) and
  `SolverMultiSetBlockMultiresolutionLocalTimestepping` (per-level steps), both in the
  `SolverMultiSetBase` hierarchy.
- Buffers: `MultiresolutionBoundary` provides per-level boundary particle overlays
  (with `MassNodes` for density/mass redistribution across levels);
  `MultiresolutionBoundaryLTS` is the LTS-specific buffer variant.
- Ghost interfaces: `GhostZone.h` classifies blocks at level interfaces;
  `BoundaryGhostParticles` holds ghost values of boundary particles on interfaces.
  On the current branch, **boundary-ghost updates moved into the WCSPH-BI model**
  (updating from the neighbor fluid is a selectable method; recent commits reworked
  `ghostBoundaryIndices` handling and split the code into `BoundaryGhostParticles.h`).
- Working examples: `examples/WCSPH-BI/*_multiresolution*/`; unit fixtures under
  `src/UnitTests/multiresolution/testConfigurations/` (dummy 2D/3D MR and MR-LTS
  simulations, each with its own `run.py` and JSONC template).
- `MultiresolutionBuffer.h` and `MultiresolutionRectangleBuffer.h` are two generations
  of the same classes — do not mix them (see anti-patterns).

## NOTES

- `include/SPH/Models/WCSPH_MFD/` is archived/incomplete (only `Interactions.hpp`).
- `include/SPH/Models/WCSPH_BI/experimental/` holds non-production BI variants.
- Release builds compile with `-O3 -DNDEBUG -use_fast_math`; `RelWithDebInfo` adds `--generate-line-info` (nvcc) for `gdb`-friendly line tables.
- `build/` may contain ad-hoc experiment directories (`mr2d-*`, `*-smoke`, `ghost-test`, `opencode-tests`) and scratch configs — they are not build targets.
- Root-level scratch files (e.g. `task`, stray smoke-test directories) are personal notes, not project files.
- Executables land in `build/` mirroring the source tree (e.g. `build/examples/WCSPH-DBC/damBreak2D_WCSPH-DBC/damBreak2D_WCSPH-DBC_cuda`); `run.py` computes that path from its own location.
- TNL's `HeadersOnly` mode is not the target configuration: TNL-SPH requires CUDA, a CUDA-aware MPI (OpenMPI) for distributed runs, and Python 3 with NumPy/VTK for the init/run tooling.
- Generated content (`sources/`, `results/`, VTK, logs, `__pycache__`) is covered by `.gitignore`; the CI image (`archlinux-tnl-cuda`) assumes Ninja, but Makefiles work locally.
