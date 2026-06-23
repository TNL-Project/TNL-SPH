from sph_codegen.scheme import (
    FluidSet, BoundarySet, ConstantSet,
    ScalarField, VectorField,
    interactions, integration,
    dot,
    gradW, r_ij, drs, t, dt,
)

class Fluid(FluidSet):
    rho  = ScalarField(t=2, read=True, write=True)
    drho = ScalarField()
    p = ScalarField(write=True, eos=True)
    v = VectorField(t=2, read=True, write=True)
    a = VectorField()

class Boundary(BoundarySet):
    rho = ScalarField(t=2)
    drho = ScalarField()
    p = ScalarField(eos=True)
    v = VectorField(read=True)

class Constants(ConstantSet):
    Scalars = "h dp mass delta alpha dynamicViscosity speedOfSound rho0 dtInit cfl dtMin eps"
    Vectors = "gravity"

@interactions(Fluid, Boundary, Constants)
class Interactions:

    def fluid_fluid(i, j):
        eps_h2 = h * h * eps
        v_ij = v[i] - v[j]
        drdv = dot(r_ij, v_ij)
        mu = h * drdv / (drs * drs + eps_h2)
        visco = -2 * alpha * speedOfSound * mu / (rho[i] + rho[j]) if drdv < 0 else 0
        psi = 2 * h * delta * speedOfSound * (rho[j] - rho[i]) / (drs * drs + eps_h2)

        drho[i] += dot(v_ij, gradW) * mass - psi * dot(r_ij, gradW) * mass / rho[j]
        a[i] += -((p[i] + p[j]) / (rho[i] * rho[j]) + visco) * gradW * mass

    def fluid_boundary(i, j):
        eps_h2 = h * h * eps
        v_ij = v[i] - v[j]
        drdv = dot(r_ij, v_ij)
        mu = h * drdv / (drs * drs + eps_h2)
        visco = -2 * alpha * speedOfSound * mu / (rho[i] + rho[j]) if drdv < 0 else 0
        psi = 2 * h * delta * speedOfSound * (rho[j] - rho[i]) / (drs * drs + eps_h2)

        drho[i] += dot(v_ij, gradW) * mass - psi * dot(r_ij, gradW) * mass / rho[j]
        a[i]  += -((p[i] + p[j]) / (rho[i] * rho[j]) + visco) * gradW * mass

    def boundary_fluid(i, j):
        v_ij = v[i] - v[j]
        drho[i] += dot(v_ij, gradW) * mass * (rho[i] / rho[j])

@integration(Fluid, Boundary, Constants)
class Integration:

    def fluid_update(i):
        if step % 40 == 0:
            v[i](t+dt) << v[i](t) + dt * a[i](t)
            rho[i](t+dt) << rho[i](t) + dt * drho[i](t)
        else:
            v[i](t+dt) << v[i](t-dt) + 2 * dt * a[i](t)
            rho[i](t+dt) << rho[i](t-dt) + 2 * dt * drho[i](t)

    def boundary_update(i):
        rho[i](t+dt) << rho[i](t-dt) + 2 * dt * drho[i](t)
