**WCSPH-BI** is TNL-SPH's flagship solution for free surface flow problems. It is built using the standard WCSPH model but includes renormalization and natively incorporates boundary integral formulations of solid walls (BI).

The WCSPH-BI model can solve three governing systems:

- An inviscid weakly compressible barotropic fluid in isothermal flow (**Euler equations**). The relation between density and pressure is given by the Cole equation of state. System (REF).
- A viscous compressible Newtonian barotropic fluid in isothermal flow (**Navier-Stokes equations** for compressible fluid). The relation between density and pressure is given by the Cole equation of state. System (REF).
- A viscous weakly compressible Newtonian barotropic fluid in isothermal flow (continuity equation for compressible flow with momentum conservation for inviscid fluid). The relation between density and pressure is given by the Cole equation of state. System (REF).

Table \ref{tab:wcsph_bi:variables} summarizes the particle variables defined for WCSPH-BI, including their type and software name. In addition to standard flow field variables, a renormalization factor $\gamma$ of scalar type is used, which by default stores the Shepard factor of each particle. When higher-order normalization terms (such as the Libersky matrix) are employed, $\gamma$ can be predefined as the required type. For boundary particles, their size $S$ and boundary normal $\mathbf{n}$ are included.

**Table: Variables of Fluid ($\mathcal{F}$) and Boundary ($\mathcal{B}$) particle sets in WCSPH-BI**

| Field                  | Name            | Type                | Description                                      | Units                      |
|------------------------|-----------------|---------------------|--------------------------------------------------|----------------------------|
| **Fluid particles** ($\mathcal{F}$) ||||
| $\rho$                 | `rho`           | `ScalarArrayType`   | fluid density                                    | kg·m⁻³                     |
| $\mathrm{D}\rho/\mathrm{D}t$ | `drhodt`  | `ScalarArrayType`   | time derivative of fluid density                 | kg·m⁻³·s⁻¹                 |
| $p$                    | `p`             | `ScalarArrayType`   | pressure *(this field is not actively used)*     | Pa                         |
| $\mathbf{v}$           | `v`             | `VectorArrayType`   | fluid/particle velocity                          | m·s⁻¹                      |
| $\mathrm{D}\mathbf{v}/\mathrm{D}t$ | `dvdt` | `VectorArrayType`   | time derivative of fluid velocity                | m·s⁻²                      |
| $\gamma$               | `gamma`         | `ScalarArrayType`   | particle concentration (Shepard factor)          | —                          |
| **Boundary particles** ($\mathcal{B}$) ||||
| $\rho$                 | `rho`           | `ScalarArrayType`   | fluid density                                    | kg·m⁻³                     |
| $\mathrm{D}\rho/\mathrm{D}t$ | `drhodt`  | `ScalarArrayType`   | time derivative of fluid density                 | kg·m⁻³·s⁻¹                 |
| $p$                    | `p`             | `ScalarArrayType`   | pressure *(this field is not actively used)*     | Pa                         |
| $\mathbf{v}$           | `v`             | `VectorArrayType`   | fluid/particle velocity                          | m·s⁻¹                      |
| $\mathrm{D}\mathbf{v}/\mathrm{D}t$ | `dvdt` | `VectorArrayType`   | time derivative of fluid velocity                | m·s⁻²                      |
| $\gamma$               | `gamma`         | `ScalarArrayType`   | particle concentration (Shepard factor)          | —                          |
| $\mathbf{n}$           | `n`             | `VectorArrayType`   | normal vector to boundary element                | —                          |
| $S$                    | `elementSize`   | `ScalarArrayType`   | length or area of boundary element               | —                          |

### Spatial discretization

Table below provides a list of all spatial SPH discretizations employed within the WCSPH-BI model. These discretizations are used to construct several built-in schemes:

- Inviscid WCSPH-BI Consistent scheme
- Inviscid WCSPH-BI Conservative scheme
- Weakly Compressible WCSPH-BI Consistent scheme
- Weakly Compressible WCSPH-BI Conservative scheme

which, in combination with implicit midpoint time integration, result in the **unconditionally stable SPH algorithm** proposed by Cercos-Pita (2024).

**Table: Spatial discretizations in WCSPH-BI**

| Term                                            | Bulk discretization $\langle \bullet \rangle^{\Omega}$                                                                 | Boundary discretization $\langle \bullet \rangle^{\partial\Omega}$                                                                 |
|-------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------|
| $\langle \nabla \cdot \mathbf{v} \rangle$       | $\displaystyle \langle \nabla \cdot \mathbf{v} \rangle_i^{\Omega, \text{asym.}} = \sum_{j \in \mathcal{F}} (\mathbf{v}_j - \mathbf{v}_i) \nabla_i W_{ij}^{\mathcal{L}} V_j$ | $\displaystyle \langle \nabla \cdot \mathbf{v} \rangle_i^{\partial\Omega, \text{asym.}} = \sum_{j \in \mathcal{B}} (\mathbf{v}_j - \mathbf{v}_i) \mathbf{n}_j W_{ij}^{\mathcal{L}} S_j$ <br><br> **or conservative form** <br> $\displaystyle \langle \nabla \cdot \mathbf{v} \rangle_i^{\partial\Omega, \text{conserv.}} = 2 \sum_{j \in \mathcal{B}} (\mathbf{v}_j - \mathbf{v}_i) \mathbf{n}_j W_{ij} S_j$ |
| $\langle \nabla p \rangle$                      | $\displaystyle \langle \nabla p \rangle_i^{\Omega, \text{sym.}} = \sum_{j \in \mathcal{F}} (p_i + p_j) \nabla_i W_{ij}^{\mathcal{L}} V_j$ | $\displaystyle \langle \nabla p \rangle_i^{\partial\Omega, \text{sym.}} = \sum_{j \in \mathcal{B}} (p_i + p_j) \mathbf{n}_j W_{ij}^{\mathcal{L}} S_j$ <br><br> **or conservative form** <br> $\displaystyle \langle \nabla p \rangle_i^{\partial\Omega, \text{conserv.}} = 2 p_i \sum_{j \in \mathcal{B}} \mathbf{n}_j W_{ij} S_j$ |
| $\langle \Delta \mathbf{v} \rangle$             | $\displaystyle \langle \Delta \mathbf{v} \rangle_i^{\Omega,\text{MVT}} = 2\sum_{j \in \mathcal{F}} \frac{ (\mathbf{v}_i - \mathbf{v}_j) \cdot (\mathbf{x}_i - \mathbf{x}_j) }{ \lVert \mathbf{x}_i - \mathbf{x}_j \rVert^2 } \nabla_i W_{ij} \, V_j$ | $\displaystyle \langle \Delta \mathbf{v} \rangle_i^{\partial\Omega,\text{MVT}} = 2\sum_{j \in \mathcal{B}} \frac{ \mathbf{v}_i - \mathbf{v}_j }{ (\mathbf{x}_i - \mathbf{x}_j) \cdot \mathbf{n}_j } W_{ji} S_j$ |
| $\langle \Delta \mathbf{v} + 2 \nabla (\nabla \cdot \mathbf{v}) \rangle$ | $\displaystyle \langle \Delta \mathbf{v} + 2 \nabla (\nabla \cdot \mathbf{v}) \rangle_i^{\Omega,\text{MGVT}} = 2(n+2)\mu \sum_{j \in \mathcal{F}} \frac{ (\mathbf{v}_i - \mathbf{v}_j) \cdot (\mathbf{x}_i - \mathbf{x}_j) }{ \lVert \mathbf{x}_i - \mathbf{x}_j \rVert^2 } \nabla_i W_{ij}^{\mathcal{L}} V_j$ | $\displaystyle \langle \Delta \mathbf{v} + 2 \nabla (\nabla \cdot \mathbf{v}) \rangle_i^{\partial\Omega,\text{MGVT}} = 2(n+2)\mu \sum_{j \in \mathcal{B}} \frac{ \mathbf{v}_i - \mathbf{v}_j }{ (\mathbf{x}_i - \mathbf{x}_j) \cdot \mathbf{n}_j } W_{ij}^{\mathcal{L}} S_j$ |

(Note: The viscosity terms use the **Morris** formulation for $\Delta \mathbf{v}$ and the **Monaghan-Gingold** formulation for the full compressible viscous operator $\Delta \mathbf{v} + 2 \nabla (\nabla \cdot \mathbf{v})$ respectively.)

The spatial WCSPH-BI schemes include, by default, two additional algorithms inspired by Cercos-Pita (2015).

- The **Conservative Velocity-Based Elastic Bounce** (EB-Conservative), which enforces the non-penetrability of solid walls.
- Alternatively, the simple **Elastic Bounce** (EB-ε) algorithm, which allows control of the relative energy of reflection through a parameter ε ∈ [0, 1].

The second is the **Near-Boundary Volume Correction Particle Shifting** (NBVC-PST), which adjusts the positions of particles near the boundary to prevent particles from getting too close to the wall and to ensure their geometrical volume does not penetrate the walls.

### Inviscid WCSPH-BI Conservative scheme

$$
\begin{aligned}
\frac{\mathrm{D} \mathbf{x}_i(t)}{\mathrm{D}t} &= \mathbf{v}_i(t) \\[1em]
\left\langle \frac{\mathrm{D} \rho}{\mathrm{D}t} \right\rangle_i(t) &= -\rho_i(t) \left\langle \nabla \cdot \mathbf{v} \right\rangle_i(t) + \delta_{\psi} h c_0 \mathcal{D}_i(t) \\[1em]
\left\langle \frac{\mathrm{D} \mathbf{v}}{\mathrm{D}t} \right\rangle_i(t) &= -\frac{1}{\rho_i(t)} \langle \nabla p \rangle_i(t) - \mathbf{g} \\[1em]
p_i(t) &= p_0 + c_0^2 (\rho_i(t) - \rho_0)
\end{aligned}
$$

- $\displaystyle \left\langle \nabla \cdot \mathbf{v} \right\rangle_i = \left\langle \nabla \cdot \mathbf{v} \right\rangle_i^{\Omega, \text{asym.}} + \left\langle \nabla \cdot \mathbf{v} \right\rangle_i^{\partial \Omega, \text{conserv.}}$
- $\displaystyle \left\langle \nabla p \right\rangle_i = \left\langle \nabla p \right\rangle_i^{\Omega, \text{sym.}} + \left\langle \nabla p \right\rangle_i^{\partial \Omega, \text{conserv.}}$
- Optional diffusive term $\mathcal{D}_i$
- Conservative elastic bounce (EB-Conservative)
- Near-boundary volume shifting correction (NBVC-PST)

**Caption:** Spatial scheme for (weakly) compressible inviscid barotropic fluid in isothermal flow (Euler equations), using conservative WCSPH-BI formulation.

### Weakly Compressible WCSPH-BI Conservative scheme

$$
\begin{aligned}
\frac{\mathrm{D} \mathbf{x}_i(t)}{\mathrm{D}t} &= \mathbf{v}_i(t) \\[1em]
\left\langle \frac{\mathrm{D} \rho}{\mathrm{D}t} \right\rangle_i(t) &= -\rho_i(t) \left\langle \nabla \cdot \mathbf{v} \right\rangle_i(t) + \delta_{\psi} h c_0 \mathcal{D}_i(t) \\[1em]
\left\langle \frac{\mathrm{D} \mathbf{v}}{\mathrm{D}t} \right\rangle_i(t) &= -\frac{1}{\rho_i(t)} \langle \nabla p \rangle_i(t) + \mu \left\langle \Delta \mathbf{v} + 2 \nabla (\nabla \cdot \mathbf{v}) \right\rangle_i(t) - \mathbf{g} \\[1em]
p_i(t) &= p_0 + c_0^2 (\rho_i(t) - \rho_0)
\end{aligned}
$$

- $\displaystyle \left\langle \nabla \cdot \mathbf{v} \right\rangle_i = \left\langle \nabla \cdot \mathbf{v} \right\rangle_i^{\Omega, \text{asym.}} + \left\langle \nabla \cdot \mathbf{v} \right\rangle_i^{\partial \Omega, \text{conserv.}}$
- $\displaystyle \left\langle \nabla p \right\rangle_i = \left\langle \nabla p \right\rangle_i^{\Omega, \text{sym.}} + \left\langle \nabla p \right\rangle_i^{\partial \Omega, \text{conserv.}}$
- $\displaystyle \left\langle \Delta \mathbf{v} + 2 \nabla (\nabla \cdot \mathbf{v}) \right\rangle_i = \left\langle \Delta \mathbf{v} + 2 \nabla (\nabla \cdot \mathbf{v}) \right\rangle_i^{\Omega, \text{MG}} + \left\langle \Delta \mathbf{v} + 2 \nabla (\nabla \cdot \mathbf{v}) \right\rangle_i^{\partial \Omega, \text{MG,no-slip}}$
- Optional diffusive term $\mathcal{D}_i$
- Conservative elastic bounce (EB-Conservative)
- Near-boundary volume shifting correction (NBVC-PST)

**Caption:** Spatial scheme for (weakly) compressible Newtonian barotropic fluid in isothermal flow (Navier-Stokes equations), using conservative WCSPH-BI formulation.

### Inviscid WCSPH-BI Consistent scheme

$$
\begin{aligned}
\frac{\mathrm{D} \mathbf{x}_i(t)}{\mathrm{D}t} &= \mathbf{v}_i(t) \\[1em]
\left\langle \gamma \right\rangle_i(t) &= \sum_{j \in \mathcal{F}} W_{ij} V_j \\[1em]
\left\langle \frac{\mathrm{D} \rho}{\mathrm{D}t} \right\rangle_i(t) &= \frac{1}{\left\langle \gamma \right\rangle_i(t)} \left( -\rho_i(t) \left\langle \nabla \cdot \mathbf{v} \right\rangle_i(t) + \delta_{\psi} h c_0 \mathcal{D}_i(t) \right) \\[1em]
\left\langle \frac{\mathrm{D} \mathbf{v}}{\mathrm{D}t} \right\rangle_i(t) &= -\frac{1}{\rho_i(t) \left\langle \gamma \right\rangle_i(t)} \langle \nabla p \rangle_i(t) - \mathbf{g} \\[1em]
p_i(t) &= p_0 + c_0^2 (\rho_i(t) - \rho_0)
\end{aligned}
$$

- $\displaystyle \left\langle \nabla \cdot \mathbf{v} \right\rangle_i = \left\langle \nabla \cdot \mathbf{v} \right\rangle_i^{\Omega, \text{asym.}} + \left\langle \nabla \cdot \mathbf{v} \right\rangle_i^{\partial \Omega, \text{asym.}}$
- $\displaystyle \left\langle \nabla p \right\rangle_i = \left\langle \nabla p \right\rangle_i^{\Omega, \text{sym.}} + \left\langle \nabla p \right\rangle_i^{\partial \Omega, \text{sym.}}$
- Optional diffusive term $\mathcal{D}_i$
- Conservative elastic bounce (EB-Conservative)
- Near-boundary volume shifting correction (NBVC-PST)

**Caption:** Consistent WCSPH-BI scheme for (weakly) compressible inviscid barotropic fluid.

### Weakly Compressible WCSPH-BI Consistent scheme

$$
\begin{aligned}
\frac{\mathrm{D} \mathbf{x}_i(t)}{\mathrm{D}t} &= \mathbf{v}_i(t) \\[1em]
\left\langle \gamma \right\rangle_i(t) &= \sum_{j \in \mathcal{F}} W_{ij} V_j \\[1em]
\left\langle \frac{\mathrm{D} \rho}{\mathrm{D}t} \right\rangle_i(t) &= \frac{1}{\left\langle \gamma \right\rangle_i(t)} \left( -\rho_i(t) \left\langle \nabla \cdot \mathbf{v} \right\rangle_i(t) + \delta_{\psi} h c_0 \mathcal{D}_i(t) \right) \\[1em]
\left\langle \frac{\mathrm{D} \mathbf{v}}{\mathrm{D}t} \right\rangle_i(t) &= \frac{1}{\left\langle \gamma \right\rangle_i(t)} \left( -\frac{1}{\rho_i(t)} \langle \nabla p \rangle_i(t) + \mu \left\langle \Delta \mathbf{v} + 2 \nabla (\nabla \cdot \mathbf{v}) \right\rangle_i(t) \right) - \mathbf{g} \\[1em]
p_i(t) &= p_0 + c_0^2 (\rho_i(t) - \rho_0)
\end{aligned}
$$

- $\displaystyle \left\langle \nabla \cdot \mathbf{v} \right\rangle_i = \left\langle \nabla \cdot \mathbf{v} \right\rangle_i^{\Omega, \text{asym.}} + \left\langle \nabla \cdot \mathbf{v} \right\rangle_i^{\partial \Omega, \text{asym.}}$
- $\displaystyle \left\langle \nabla p \right\rangle_i = \left\langle \nabla p \right\rangle_i^{\Omega, \text{sym.}} + \left\langle \nabla p \right\rangle_i^{\partial \Omega, \text{sym.}}$
- $\displaystyle \left\langle \Delta \mathbf{v} + 2 \nabla (\nabla \cdot \mathbf{v}) \right\rangle_i = \left\langle \Delta \mathbf{v} + 2 \nabla (\nabla \cdot \mathbf{v}) \right\rangle_i^{\Omega, \text{MG}} + \left\langle \Delta \mathbf{v} + 2 \nabla (\nabla \cdot \mathbf{v}) \right\rangle_i^{\partial \Omega, \text{MG,no-slip}}$
- Optional diffusive term $\mathcal{D}_i$
- Conservative elastic bounce (EB-Conservative)
- Near-boundary volume shifting correction (NBVC-PST)

**Caption:** Consistent WCSPH-BI scheme for (weakly) compressible Newtonian barotropic fluid.

### Time integration – Implicit Mid-Point scheme (IMP)

To compute the time derivatives at $t_{n+1/2}$, the implicit midpoint scheme includes an inner iteration loop (indexed by superscript $m$):

$$
\begin{aligned}
\hat{\rho}_i^m (t_{n+1/2}) &= \rho_i^n + \frac{\Delta t}{2} \widetilde{\rho}_i^m \\[1em]
\widetilde{\rho}_i^{m+1} &= f^m \left\langle \frac{\mathrm{D} \rho}{\mathrm{D}t} \right\rangle_i (t_{n+1/2}) + (1 - f^m) \widetilde{\rho}_i^m
\end{aligned}
$$

$$
\begin{aligned}
\hat{\mathbf{v}}_i^m (t_{n+1/2}) &= \mathbf{v}_i^n + \frac{\Delta t}{2} \widetilde{\mathbf{v}}_i^m \\[1em]
\widetilde{\mathbf{v}}_i^{m+1} &= f^m \left\langle \frac{\mathrm{D} \mathbf{v}}{\mathrm{D}t} \right\rangle_i (t_{n+1/2}) + (1 - f^m) \widetilde{\mathbf{v}}_i^m
\end{aligned}
$$

where $f^m \in [0,1]$ is a relaxation factor.

The iteration is controlled by the residual of energy conservation:

$$
R_{\Delta t}[E](t) = \sum_{i \in \mathcal{F}} \left| m_i \mathbf{v}_i(t) \cdot \left\langle \frac{\mathrm{D}\mathbf{v}}{\mathrm{D}t} \right\rangle_i(t) \right| + \sum_{i \in \mathcal{F}} \left| \frac{m_i p_i(t)}{\rho_i^2(t)} \left\langle \frac{\mathrm{D} \rho}{\mathrm{D}t} \right\rangle_i(t) \right|
$$

At the first iteration, relaxation is not applied ($f^0 = 1$). The relaxation factor is progressively increased as $f^{m+1} = \Delta f + (1 - \Delta f) f^m$ if the residual reduction rate does not meet a geometrical threshold $r_R$ (i.e., if $r_R > R_{\Delta t}^{m}[E] / R_{\Delta t}^{m-1}[E]$). Both $\Delta f$ and $r_R$ are tunable parameters.

### Time integration – Verlet scheme

$$
\begin{aligned}
\mathbf{r}_i(t_{n+1}) &= \mathbf{r}_i(t_{n-1}) + \Delta t \mathbf{v}_i(t_n) + \frac{1}{2} \Delta t^2 \left\langle \frac{\mathrm{D} \mathbf{v}}{\mathrm{D}t} \right\rangle_i(t_n) \\[1em]
\mathbf{v}_i(t_{n+1}) &= \mathbf{v}_i(t_{n-1}) + 2 \Delta t \left\langle \frac{\mathrm{D} \mathbf{v}}{\mathrm{D}t} \right\rangle_i(t_n) \\[1em]
\rho_i(t_{n+1}) &= \rho_i(t_{n-1}) + 2 \Delta t \left\langle \frac{\mathrm{D} \rho}{\mathrm{D}t} \right\rangle_i(t_n)
\end{aligned}
$$

Due to staggered integration intervals, the scheme can become unstable. To stabilize it, an intermediate explicit Euler step is performed every $f_{ES}$-th step (default $f_{ES} = 40$):

$$
\begin{aligned}
\mathbf{r}_i(t_{n+1}) &= \mathbf{r}_i(t_{n-1}) + \Delta t \mathbf{v}_i(t_n) + \frac{1}{2} \Delta t^2 \left\langle \frac{\mathrm{D} \mathbf{v}}{\mathrm{D}t} \right\rangle_i(t_n) \\[1em]
\mathbf{v}_i(t_{n+1}) &= \mathbf{v}_i(t_n) + \Delta t \left\langle \frac{\mathrm{D} \mathbf{v}}{\mathrm{D}t} \right\rangle_i(t_n) \\[1em]
\rho_i(t_{n+1}) &= \rho_i(t_n) + \Delta t \left\langle \frac{\mathrm{D} \rho}{\mathrm{D}t} \right\rangle_i(t_n)
\end{aligned}
$$

### Time integration – Symplectic Verlet scheme

The symplectic position Verlet method is a two-stage predictor-corrector scheme that is time-reversible and preserves geometric structure when applied to suitable systems.

**Predictor:**

$$
\begin{aligned}
\mathbf{r}_i(t_{n+1/2}) &= \mathbf{r}_i(t_n) + \frac{1}{2} \Delta t \mathbf{v}_i(t_n) \\[1em]
\mathbf{v}_i(t_{n+1/2}) &= \mathbf{v}_i(t_n) + \frac{1}{2} \Delta t \left\langle \frac{\mathrm{D} \mathbf{v}}{\mathrm{D}t} \right\rangle_i(t_n) \\[1em]
\rho_i(t_{n+1/2}) &= \rho_i(t_n) + \frac{1}{2} \Delta t \left\langle \frac{\mathrm{D} \rho}{\mathrm{D}t} \right\rangle_i(t_n)
\end{aligned}
$$

**Corrector:**

$$
\begin{aligned}
\mathbf{r}_i(t_{n+1}) &= \mathbf{r}_i(t_{n-1}) + \Delta t \frac{\mathbf{v}_i(t_n) + \mathbf{v}_i(t_{n+1})}{2} \\[1em]
\mathbf{v}_i(t_{n+1}) &= \mathbf{v}_i(t_n) + \Delta t \left\langle \frac{\mathrm{D} \mathbf{v}}{\mathrm{D}t} \right\rangle_i(t_{n+1/2}) \\[1em]
\rho_i(t_{n+1}) &= \rho_i(t_n) + \frac{2 - \epsilon}{2 + \epsilon}, \qquad \epsilon = -\Delta t \frac{1}{\rho_i(t_{n+1/2})} \left\langle \frac{\mathrm{D} \rho}{\mathrm{D}t} \right\rangle_i(t_{n+1/2})
\end{aligned}
$$

### Interaction functions

WCSPH-BI is designed as a model for a multi-set SPH solver, requiring three interaction functions:

- The **initialize interaction** function $\mathscr{F}_{\text{initialize}}$ precomputes variables required for the interaction, if necessary.
- The main **interact** function $\mathscr{F}$ handles the core particle interactions.
- The **finalize interaction** function $\mathscr{F}_{\text{finalize}}$ performs post-interaction routines.

**Algorithm: Initialize interaction** $\mathscr{F}_{\text{initialize}}$

- *empty function*

**Algorithm: Interact** $\mathscr{F}$

1. Loop over boundary particles: $\forall i \in \mathcal{B}$ *(skipped if not required by selected scheme)*
   - Loop over fluid neighbors: $\forall j \in \mathcal{N}_i(\mathcal{F})$ compute contribution to $\rho_i$

2. Loop over fluid particles: $\forall i \in \mathcal{F}$
   - Loop over fluid neighbors: $\forall j \in \mathcal{N}_i(\mathcal{F})$ compute contribution to
     $\left\langle \frac{\mathrm{D}\rho}{\mathrm{D}t} \right\rangle_i$, $\left\langle \frac{\mathrm{D}\mathbf{v}}{\mathrm{D}t} \right\rangle_i$ and $\left\langle \gamma \right\rangle_i$
   - Loop over boundary neighbors: $\forall j \in \mathcal{N}_i(\mathcal{B})$ compute contribution to
     $\left\langle \frac{\mathrm{D}\rho}{\mathrm{D}t} \right\rangle_i$ and $\left\langle \frac{\mathrm{D}\mathbf{v}}{\mathrm{D}t} \right\rangle_i$

**Algorithm: Finalize interaction** $\mathscr{F}_{\text{finalize}}$

1. Loop over fluid particles: $\forall i \in \mathcal{F}$
   - Loop over fluid neighbors: $\forall j \in \mathcal{N}_i(\mathcal{F})$ compute contribution to $\rho_i$

2. Loop over fluid particles: $\forall i \in \mathcal{F}$
   - Normalize derivatives $\left\langle \frac{\mathrm{D}\rho}{\mathrm{D}t} \right\rangle_i$ and $\left\langle \frac{\mathrm{D}\mathbf{v}}{\mathrm{D}t} \right\rangle_i$ with $\left\langle \gamma \right\rangle_i$
   - Add acceleration due to external forces $\mathbf{f}_i$

### Runtime parameters of WCSPH-BI model

| Var              | Variable Name                        | Type          | Commentary                                              | Default Value |
|------------------|--------------------------------------|---------------|---------------------------------------------------------|---------------|
| $\Delta x$       | `dp`                                 | `RealType`    | initial particle distance [m]                           | 0.f           |
| $h$              | `h`                                  | `RealType`    | smoothing length [m]                                    | 0.f           |
| $r_{\text{cut}}$ | `searchRadius`                       | `RealType`    | radius of kernel support [m]                            | 0.f           |
| $m_{\mathcal{F}}$| `mass`                               | `RealType`    | particle mass [kg]                                      | 0.f           |
| $m_{\mathcal{B}}$| `massBoundary`                       | `RealType`    | boundary particle mass [kg]                             | 0.f           |
| $S$              | `boundaryElementSize`                | `RealType`    | size of boundary element [m]                            | 0.f           |
| $\delta_{\psi}$  | `delta`                              | `RealType`    | coefficient of diffusive term (DT) [-]                  | 0.1f          |
| $\alpha$         | `alpha`                              | `RealType`    | parameter of artificial viscosity [-]                   | 0.02f         |
| $\mu$            | `dynamicViscosity`                   | `RealType`    | value of dynamic viscosity [Pa·s]                       | 0.f           |
|                  | `scaleBVTCoef`                       | `RealType`    | coefficient used to tune BVT effect to physical behavior [-] | 1.f     |
| $c_s$            | `speedOfSound`                       | `RealType`    | numerical speed of sound used in SPH calculations [m/s] | 0.f           |
|                  | `coefB`                              | `RealType`    | coefficient of the Tait equation of state               | 0.f           |
| $\rho_0$         | `rho0`                               | `RealType`    | referential density of the fluid [kg/m³]                | 0.f           |
| $\varepsilon$    | `enableElasticBounce`                | `bool`        | enable elastic bounce boundary correction               | true          |
| $\epsilon_{\text{bounce}}$ | `elasticFactor`                | `RealType`    | elastic boundary correction factor                      | 1.f           |
| $r_{\text{box}}$ | `r_boxFactor`                        | `RealType`    | $r_{\text{box}} = r_{\text{boxFactor}} \times dp$       | 1.5f          |
| $h_{\text{fac}}$ | `minimalDistanceFactor`              | `RealType`    | minimal distance factor                                 | 0.5f          |
| $\Delta t_0$     | `dtInit`                             | `RealType`    | initial time step [s]                                   | 0.f           |
| CFL              | `cfl`                                | `RealType`    | CFL number [-]                                          | 0.f           |
| $\Delta t_{\min}$| `dtMin`                              | `RealType`    | minimal allowed time step [s]                           | 0.f           |
| $n_{\max iter}$  | `midpointMaxInterations`             | `int`         | max number of midpoint iterations [-]                   | 50            |
| $R_{\Delta t,\text{tol.}}$ | `midpointResidualTolerance`    | `RealType`    | midpoint iteration tolerance to reach                   | 1.e-7         |
| $f$              | `midpointRelaxCoef`                  | `RealType`    | midpoint relaxation coefficient [-]                     | 0             |
| $f_0$            | `midpointRelaxCoef_0`                | `RealType`    | midpoint relaxation coefficient in first iteration [-]  | 0             |
| $g_R$            | `midpointResidualMinimalDecay`       | `RealType`    | midpoint minimal decay of residual [-]                  | 0.2f          |
| $\Delta f$       | `midpointRelaxCoefIncrement`         | `RealType`    | midpoint increment of relaxation coefficient [-]        | 0             |
| $\mathbf{g}$     | `gravity`                            | `VectorType`  | external forces [m/s²]                                  | 0.f           |
| $\varepsilon$    | `eps`                                | `RealType`    | constant to prevent zero in denominator [-]             | 0.001f        |

### Modular plug-in terms provided in WCSPH-BI model

| Parameter              | Option                              | Description                                                                 |
|------------------------|-------------------------------------|-----------------------------------------------------------------------------|
| `KernelFunction`       | `WendlandKernel`                    | Wendland C² kernel                                                          |
| `DiffusiveTerm`        | `None`                              | —                                                                           |
|                        | `MolteniDiffusiveTerm`              | Diffusive term proposed by Molteni                                          |
|                        | `FourtakasDiffusiveTerm`            | Diffusive term proposed by Fourtakas                                        |
| `ViscousTerm`          | `None`                              | —                                                                           |
|                        | `ArtificialViscosity`               | Artificial viscosity (Monaghan 1985)                                        |
|                        | `PhysicalViscosity_MVT`             | Morris viscous term                                                         |
|                        | `PhysicalViscosity_MGVT`            | Monaghan & Gingold viscous term                                             |
| `EOS`                  | `ColeEquationOfState`               | Weakly compressible Cole equation of state                                  |
|                        | `LinearizedColeEquationOfState`     | Linearized weakly compressible Cole equation of state                       |
| `BCType`               | `BIConservative_numeric`            | Conservative formulation of boundary integrals                              |
|                        | `BIConsistent_numeric`              | Numeric renormalized formulation of boundary integrals                      |
| `TimeStepping`         | `ConstantTimeStep`                  | —                                                                           |
|                        | `VariableTimeStep`                  | —                                                                           |
| `IntegrationScheme`    | `MidpointScheme`                    | Implicit midpoint scheme                                                    |
|                        | `VerletScheme`                      | Verlet scheme                                                               |
|                        | `SymplecticVerletScheme`            | Symplectic Verlet scheme                                                    |
|                        | `RK45`                              | Runge–Kutta 4(5) scheme                                                     |
