## WCSPH_DBC

Basics SPH model for the [simple fluid-stationary boundary free surface flow simulation](./sph_simple_fluid.md), using the simple multi-layer boundary condition.  It implements standard WCSPH variant as well as \f$\delta-\f$WCSPH scheme with chooseable viscous term formulation.

### Equations

The governig equations consist from continuity equation for compresible fluid, momentum equation for incompressible viscid fluid, together with the equation for particle motion and the system is closed by relation between density and pressure where Tait equation of state or linearized Tait equation of state can be used.

$$ \frac{\text{D} \mathbf{x}_i}{\text{D} t} = \mathbf{v}_i  $$

$$ \frac{\text{D} \rho_i}{\text{D} t} = - \big\langle \rho \text{div}(\mathbf{v}) \big\rangle_i = \rho_i \sum_{j=1}^N  ( \mathbf{v}_i - \mathbf{v}_j ) \cdot \nabla W_{ij} V_j + \delta_{\psi}hc_0 \mathscr{D} $$

$$ \frac{\text{D}\mathbf{v}_i }{\text{D}t} = - \frac{1 }{\rho_i}\big\langle \nabla p  \big\rangle_{i,1} + \big\langle \nabla \cdot \tau  \big\rangle_{i}  +\mathbf{f}_i = - \frac{ 1 }{ \rho_i }\sum_{j=1}^N ( p_i + p_j ) \nabla_i W_{ij} V_j	+ \frac{ 1 }{ \rho_i }\sum_{j=1}^N \Pi_{ij}\nabla_i W_{ij} V_j + \mathbf{f}_i $$

### Model variables
#### Constants

| Variable    | Type            | doubled     | Description |
| ----------- | -----------     | ----------- | ----------- |
| dp          | RealType        | yes         | Title       |
| h           | ScalarArrayType | no          | Text        |
| mass        | ScalarArrayType | no          | Text        |

#### Variable fields

| Variable    | Type            | doubled     | Description                                                                                  |
| ----------- | -----------     | ----------- | -----------                                                                                  |
| rho         | ScalarArrayType | **yes**     | \f$ \rho \f$, density field                                                                  |
| drho        | ScalarArrayType | no          | \f$ \text{D}\rho/\text{D}{t} \f$, time material derivative of density                        |
| p           | ScalarArrayType | *inactive*  | \f$ p \f$, pressure field (used only for output)                                             |
| v           | VectorArrayType | **yes**     | \f$ \mathbf{v} \f$, velocity field                                                           |
| a           | VectorArrayType | no          | \f$  \text{D}\mathbf{v}/\text{D}{t} \f$, time material derivative of velocity (acceleration) |





| Sub-modelmodels | Type            | doubled     | Description      |
| -----------     | -----------     | ----------- | -----------      |
| DiffusiveTerm   | RealType        | yes         | $ \mathscr{D} $  |
| ViscousTerm     | ScalarArrayType | no          | \f$  \Pi  \f$    |
| EOS             | ScalarArrayType | no          | Text             |
