# SPH Code Generator

Generates C++ SPH model headers (`Variables.h`, `control.h`, `Interactions.h`, `Interactions.hpp`) from a Python scheme definition. The generated code plugs directly into `SPHMultiset_CFD<Model>`.

## Quick Start

```bash
cd codegen
python sph_codegen/run.py scheme-delta-wcsph.py --output ../../include/SPH/Models/generated
```

Rebuild the example:

```bash
cmake --build build --target damBreak2D_generated_cuda
```

## How It Works

The generator has three layers:

1. **DSL** (`dsl2.py`) — Expr dataclasses representing the symbolic tree (`BinOp`, `Dot`, `Ternary`, `FieldRef`, `ConstRef`, etc.). Pure data, no side effects.

2. **Parser** (`scheme.py`) — Parses Python function bodies via `ast.parse()` into Expr trees. No execution, no monkey-patching. Two decorators:
   - `@interactions` — parses neighbor-pair functions (`fluid_fluid`, `fluid_boundary`, `boundary_fluid`)
   - `@integration` — parses time-stepping rules (`fluid_update`, `boundary_update`)

3. **Emitter** (`emit2.py` + `gen2.py`) — Walks Expr trees to C++ strings and assembles the four output files.

## Writing a Scheme

A scheme file defines five classes: `Fluid`, `Boundary`, `Constants`, `Interactions`, `Integration`.

### Fields

```python
class Fluid(FluidSet):
    rho  = ScalarField(t=2, read=True, write=True)
    drho = ScalarField()
    p    = ScalarField(write=True, eos=True)
    v    = VectorField(t=2, read=True, write=True)
    a    = VectorField()
```

| Parameter | Meaning |
|---|---|
| `t=2` | Needs swap array for sorting (time-coupled field) |
| `read=True` | Loaded from input VTK files |
| `write=True` | Written to output VTK files |
| `eos=True` | Computed from EOS (`DensityToPressure`), not read from view |

### Constants

```python
class Constants(ConstantSet):
    Scalars = "h dp mass delta alpha speedOfSound rho0"
    Vectors = "gravity"
```

Space-separated names. These become `modelParams` members in C++.

### Interactions

Function names encode the particle-set pair: `fluid_fluid`, `fluid_boundary`, `boundary_fluid`.

```python
@interactions(Fluid, Boundary, Constants)
class Interactions:

    def fluid_fluid(i, j):
        eps_h2 = h * h * eps                          # local temp (plain assignment)
        v_ij   = v[i] - v[j]                          # field read via [i] / [j]
        drdv   = dot(r_ij, v_ij)                      # dot product
        visco  = -2 * alpha * mu / (rho[i] + rho[j]) if drdv < 0 else 0  # ternary

        drho[i] += dot(v_ij, gradW) * mass             # accumulation
        a[i]    += -(p_term + visco) * gradW * mass
```

**Syntax rules:**

| Construct | Syntax | Notes |
|---|---|---|
| Field read | `v[i]`, `rho[j]` | `i` = current particle, `j` = neighbor |
| Field write | `drho[i] += expr` | Only `+=` (accumulation over neighbors) |
| Local temp | `name = expr` | Plain Python assignment |
| Dot product | `dot(a, b)` | |
| Conditional | `yes if cond else no` | Native Python conditional expression |
| Constants | `h`, `mass`, `alpha`, ... | Used by bare name |
| Builtins | `r_ij`, `drs`, `gradW` | Predefined per-neighbor quantities |

### Integration

Function names: `fluid_update`, `boundary_update`.

```python
@integration(Fluid, Boundary, Constants)
class Integration:

    def fluid_update(i):
        if step % 40 == 0:
            v[i](t+dt)   << v[i](t)   + dt * a[i](t)
            rho[i](t+dt) << rho[i](t) + dt * drho[i](t)
        else:
            v[i](t+dt)   << v[i](t-dt)   + 2 * dt * a[i](t)
            rho[i](t+dt) << rho[i](t-dt) + 2 * dt * drho[i](t)
```

**Syntax rules:**

| Construct | Syntax | Notes |
|---|---|---|
| Time reference (LHS) | `v[i](t+dt)` | Future value — assignment target |
| Time reference (RHS) | `v[i](t)`, `v[i](t-dt)` | Current / previous value |
| Assignment operator | `<<` | Required — Python forbids `=` on call LHS |
| Conditionals | `if step % 40 == 0:` | Standard Python if/else blocks |

### Built-in Quantities

Available inside interaction functions without declaration:

| Name | Type | Meaning |
|---|---|---|
| `r_ij` | Vector | Distance vector `r_i - r_j` |
| `drs` | Scalar | Scalar distance `‖r_ij‖` |
| `gradW` | Vector | Kernel gradient `r_ij * F(drs, h)` |

### Available Helpers

```python
from sph_codegen.scheme import dot, norm, powf, gradW, r_ij, drs, t, dt
```

## Generated Output

| File | Contents |
|---|---|
| `Variables.h` | `FluidVariables`, `BoundaryVariables`, `OpenBoundaryVariables` — field arrays, setSize, sort, read/write VTK, MPI sync |
| `control.h` | `WCSPH_DBCConfig` — model parameters, config setup, init, `writePrologModel` |
| `Interactions.h` | Class declaration — type aliases, method signatures, empty stubs |
| `Interactions.hpp` | CUDA kernel implementations — `interaction()`, `updateSolidBoundary()`, `computePressureFromDensity()`, periodic patches |
