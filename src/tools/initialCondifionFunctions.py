def hydrostatic_density(ry: float, fluid_height: float,
                        rho0: float, c: float) -> float:
    p = rho0 * 9.81 * (fluid_height - ry)
    return ((p / (c**2 * rho0 / 7) + 1) ** (1.0 / 7.0)) * rho0
