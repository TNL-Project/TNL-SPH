#! /usr/bin/env python3

import numpy as np
import sys
import os
import math
import argparse
from pprint import pprint

sys.path.append('../../../src/tools')
import saveParticlesVTK

def generate_moving_square_fluid_particles(setup):
    fluid_rx = []; fluid_ry = []
    dp = setup['dp']
    L = setup['channel_length']
    H = setup['channel_height']
    D = setup['square_side']
    x_square = setup['x_square']
    y_square = setup['y_square']
    U = setup['inlet_velocity']
    rho0 = setup['density']

    # Fluid domain: from x = 0 to x = L, y = 0 to y = H
    nx = int(round(L / dp))
    ny = int(round(H / dp))

    half_D = 0.5 * D

    for i in range(nx):
        x = (i + 0.5) * dp
        for j in range(ny):
            y = (j + 0.5) * dp

            # Skip particles inside the square obstacle
            if (abs(x - x_square) < half_D) and (abs(y - y_square) < half_D):
                continue

            fluid_rx.append(x)
            fluid_ry.append(y)

    fluid_n = len(fluid_rx)
    fluid_r = np.array((fluid_rx, fluid_ry, np.zeros(fluid_n)), dtype=float).T
    # Inlet velocity in x-direction
    fluid_v = np.zeros((fluid_n, 3))
    fluid_v[:, 0] = 0
    fluid_rho = rho0 * np.ones(fluid_n)
    fluid_p = np.zeros(fluid_n)
    fluid_ptype = np.zeros(fluid_n)  # fluid particles

    fluid_poly = saveParticlesVTK.create_pointcloud_polydata(
        fluid_r, fluid_v, fluid_rho, fluid_p, fluid_ptype)
    saveParticlesVTK.save_polydata(fluid_poly, "sources/moving_square_fluid.vtk")

    print(f"Generated {fluid_n} fluid particles")
    setup["fluid_n"] = fluid_n


def generate_moving_square_boundary_particles(setup):
    box_rx = []; box_ry = []
    normal_x = []; normal_y = []
    particle_type = []

    dp = setup['dp']
    L = setup['channel_length']
    H = setup['channel_height']
    D = setup['square_side']
    x_square = setup['x_square']
    y_square = setup['y_square']
    rho0 = setup['density']
    n_layers = setup['n_boundary_layers']

    half_D = 0.5 * D

    # Outer channel walls (bottom, top, inlet, outlet)
    # Bottom wall (y = 0)
    nx = int(round(L / dp))
    for i in range(0, nx):
        x = i * dp + dp / 2
        box_rx.append(x)
        box_ry.append(0.0)
        normal_x.append(0.0)
        normal_y.append(1.0)
        particle_type.append(0)

    # Top wall (y = H)
    for i in range(0, nx):
        x = i * dp + dp / 2
        box_rx.append(x)
        box_ry.append(H)
        normal_x.append(0.0)
        normal_y.append(-1.0)
        particle_type.append(0)

    # Inlet (left wall, x = 0) – only between bottom and top
    ny = int(round(H / dp))
    for j in range(0, ny):
        y = j * dp + dp / 2
        box_rx.append(0.0)
        box_ry.append(y)
        normal_x.append(1.0)
        normal_y.append(0.0)
        particle_type.append(0)

    # Outlet (right wall, x = L)
    for j in range(0, ny):
        y = j * dp + dp / 2
        box_rx.append(L)
        box_ry.append(y)
        normal_x.append(-1.0)
        normal_y.append(0.0)
        particle_type.append(0)

    # Moving square obstacle (four sides)
    # Bottom side of square
    nx_side = int(round(D / dp))
    #for i in range(nx_side + 1):
    for i in range(0, nx_side):
        x = x_square - half_D + i * dp + dp / 2
        y = y_square - half_D
        box_rx.append(x)
        box_ry.append(y)
        normal_x.append(0.0)
        normal_y.append(-1.0)
        particle_type.append(1)

    # Top side of square
    #for i in range(nx_side + 1):
    for i in range(0, nx_side):
        x = x_square - half_D + i * dp + dp / 2
        y = y_square + half_D
        box_rx.append(x)
        box_ry.append(y)
        normal_x.append(0.0)
        normal_y.append(1.0)
        particle_type.append(1)

    # Left side of square
    #for j in range(nx_side + 1):
    for j in range(0, nx_side):
        x = x_square - half_D
        y = y_square - half_D + j * dp + dp / 2
        box_rx.append(x)
        box_ry.append(y)
        normal_x.append(-1.0)
        normal_y.append(0.0)
        particle_type.append(1)

    # Right side of square
    #for j in range(nx_side + 1):
    for j in range(0, nx_side):
        x = x_square + half_D
        y = y_square - half_D + j * dp + dp / 2
        box_rx.append(x)
        box_ry.append(y)
        normal_x.append(1.0)
        normal_y.append(0.0)
        particle_type.append(1)

    boundary_n = len(box_rx)
    boundary_r = np.array((box_rx, box_ry, np.zeros(boundary_n)), dtype=float).T
    boundary_normal = np.array((normal_x, normal_y, np.zeros(boundary_n)), dtype=float).T
    boundary_v = np.zeros((boundary_n, 3))
    boundary_rho = rho0 * np.ones(boundary_n)
    boundary_p = np.zeros(boundary_n)
    boundary_elementSize = dp * np.ones(boundary_n)
    boundary_ptype = np.array(particle_type)  # boundary particles

    boundary_poly = saveParticlesVTK.create_pointcloud_polydata(
        boundary_r, boundary_v, boundary_rho, boundary_p, boundary_ptype,
        normals=boundary_normal, elementSize=boundary_elementSize)
    saveParticlesVTK.save_polydata(boundary_poly, "sources/moving_square_boundary.vtk")

    print(f"Generated {boundary_n} boundary particles")
    setup["boundary_n"] = boundary_n

    # Domain limits (for later expansion)
    setup["domain_origin_x"] = min(box_rx)
    setup["domain_origin_y"] = min(box_ry)
    setup["domain_end_x"] = max(box_rx)
    setup["domain_end_y"] = max(box_ry)


def compute_domain_size(setup):
    search_radius = setup["search_radius"]
    eps = 1.2  # safety margin
    domain_origin_x = setup["domain_origin_x"] - eps * search_radius
    domain_origin_y = setup["domain_origin_y"] - eps * search_radius
    domain_end_x = setup["domain_end_x"] + eps * search_radius
    domain_end_y = setup["domain_end_y"] + eps * search_radius

    domain_size_x = domain_end_x - domain_origin_x
    domain_size_y = domain_end_y - domain_origin_y

    extra = {
        "domain_origin_x": domain_origin_x,
        "domain_origin_y": domain_origin_y,
        "domain_size_x": domain_size_x,
        "domain_size_y": domain_size_y,
    }
    setup.update(extra)


def write_simulation_params(setup):
    with open('template/config_template.ini', 'r') as file:
        config_file = file.read()

    config_file = config_file.replace('placeholderSearchRadius', f'{setup["search_radius"]}')
    config_file = config_file.replace('placeholderDomainOrigin-x', f'{setup["domain_origin_x"]:.5f}')
    config_file = config_file.replace('placeholderDomainOrigin-y', f'{setup["domain_origin_y"]:.5f}')
    config_file = config_file.replace('placeholderDomainSize-x', f'{setup["domain_size_x"]:.5f}')
    config_file = config_file.replace('placeholderDomainSize-y', f'{setup["domain_size_y"]:.5f}')

    config_file = config_file.replace('placeholderInitParticleDistance', f'{setup["dp"]}')
    config_file = config_file.replace('placeholderSmoothingLength', f'{setup["smoothing_length"]}')
    config_file = config_file.replace('placeholderMass', f'{setup["particle_mass"]}')
    config_file = config_file.replace('placeholderBoundaryElementSize', f'{setup["dp"]}')
    config_file = config_file.replace('placeholderSpeedOfSound', f'{setup["speed_of_sound"]}')
    config_file = config_file.replace('placeholderDensity', f'{setup["density"]}')
    config_file = config_file.replace('placeholderTimeStep', f'{setup["time_step"]}')
    config_file = config_file.replace('placeholderCFL', f'{setup["cfl"]}')
    config_file = config_file.replace('placeholderAlpha', f'{setup["alpha"]}')
    config_file = config_file.replace('placeholderBackroundPressure', f'{setup["background_pressure"]}')
    config_file = config_file.replace('placeholderDynamicVicosity', f'{setup["dynamic_viscosity"]}')
    config_file = config_file.replace('placeholderFluidParticles', f'{setup["fluid_n"]}')
    config_file = config_file.replace('placeholderAllocatedFluidParticles', f'{setup["fluid_n"]}')
    config_file = config_file.replace('placeholderBoundaryParticles', f'{setup["boundary_n"]}')
    config_file = config_file.replace('placeholderAllocatedBoundaryParticles', f'{setup["boundary_n"]}')

    with open('sources/config.ini', 'w') as file:
        file.write(config_file)

    # Header file (if needed)
    with open('template/config_template.h', 'r') as file:
        config_file = file.read()

    config_file = config_file.replace('#placeholderBoundaryConditionsType', setup["bc_type"])
    config_file = config_file.replace('#placeholderBoundaryCorrection', setup["bc_correction"])
    config_file = config_file.replace('#placeholderDiffusiveTerm', setup["diffusive_term"])
    config_file = config_file.replace('#placeholderViscosTerm', setup["viscous_term"])
    config_file = config_file.replace('#placeholderTimeIntegration', setup["time_integration"])
    config_file = config_file.replace('#placeholderPst', setup["time_integration"])

    with open('template/config.h', 'w') as file:
        file.write(config_file)


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="Moving square in channel (SPHERIC benchmark 06) initial condition generator")

    g = argparser.add_argument_group("resolution parameters")
    g.add_argument("--dp", type=float, default=0.0125, help="initial particle spacing (default corresponds to ny≈400 in original script)")
    g.add_argument("--h-coef", type=float, default=4.0, help="smoothing length coefficient h = h_coef * dp")
    g = argparser.add_argument_group("geometry parameters")
    g.add_argument("--channel-length", type=float, default=10.0, help="channel length L (in units of D)")
    g.add_argument("--channel-height", type=float, default=5.0, help="channel height H (in units of D)")
    g.add_argument("--square-side", type=float, default=1.0, help="side length of the square obstacle D")
    g.add_argument("--x-square", type=float, default=1.5, help="x-position of square center (default 0.15*L with L=10*D)")
    g.add_argument("--y-square", type=float, default=2.5, help="y-position of square center (default 0.5*H)")
    g = argparser.add_argument_group("flow parameters")
    g.add_argument("--inlet-velocity", type=float, default=1.0, help="inlet velocity U")
    g.add_argument("--reynolds", type=float, default=100.0, help="Reynolds number Re = rho*U*D/nu")
    g.add_argument("--density", type=float, default=1.0, help="reference fluid density")
    g.add_argument("--speed-of-sound", type=float, default=50.0, help="speed of sound (multiples of U)")
    g.add_argument("--cfl", type=float, default=0.1, help="CFL number")
    g.add_argument("--alpha", type=float, default=0.0, help="artificial viscosity alpha")
    g.add_argument("--background-pressure", type=float, default=0.0, help="background pressure")
    g = argparser.add_argument_group("boundary and numerical options")
    g.add_argument("--bc-type", type=str, default="BIConsistent_numeric", help="boundary condition type")
    g.add_argument("--bc-correction", type=str, default="ElasticBounce", help="boundary correction")
    g.add_argument("--diffusive-term", type=str, default="MolteniDiffusiveTerm", help="density diffusion term")
    g.add_argument("--viscous-term", type=str, default="PhysicalViscosity_MVT", help="viscosity formulation")
    g.add_argument("--time-integration", type=str, default="VerletScheme", help="time integration scheme")
    g.add_argument("--pst", type=str, default="Simple", help="particles shifting scheme")

    args = argparser.parse_args()

    # Scale geometry with D
    D = args.square_side
    L = args.channel_length * D
    H = args.channel_height * D
    x_square = args.x_square * D if args.x_square > 1 else args.x_square * L
    y_square = args.y_square * D if args.y_square > 1 else args.y_square * H
    U = args.inlet_velocity

    # Dynamic viscosity from Reynolds number
    dynamic_viscosity = args.density * U * D / args.reynolds
    speed_of_sound = args.speed_of_sound * U

    setup = {
        "dp": args.dp,
        "h_coef": args.h_coef,
        "channel_length": L,
        "channel_height": H,
        "square_side": D,
        "x_square": x_square,
        "y_square": y_square,
        "inlet_velocity": U,
        "density": args.density,
        "speed_of_sound": speed_of_sound,
        "cfl": args.cfl,
        "alpha": args.alpha,
        "background_pressure": args.background_pressure,
        "dynamic_viscosity": dynamic_viscosity,
        "particle_mass": args.density * (args.dp ** 2),
        "smoothing_length": args.h_coef * args.dp,
        "search_radius": 2.0 * args.h_coef * args.dp,  # typical for quadratic kernel
        "time_step": args.cfl * (args.h_coef * args.dp) / speed_of_sound,
        "bc_type": args.bc_type,
        "bc_correction": args.bc_correction,
        "diffusive_term": args.diffusive_term,
        "viscous_term": args.viscous_term,
        "time_integration": args.time_integration,
        "pst": args.pst,
        "n_boundary_layers": 1,  # single layer boundary as in modern SPH
    }

    # Create output directories
    os.makedirs("results", exist_ok=True)
    os.makedirs("sources", exist_ok=True)

    # Generate particles
    generate_moving_square_fluid_particles(setup)
    generate_moving_square_boundary_particles(setup)

    # Compute domain
    compute_domain_size(setup)

    print("Final setup:")
    pprint(setup)

    # Write configuration files
    write_simulation_params(setup)

    print("Particle generation and configuration completed.")
