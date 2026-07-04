#! /usr/bin/env python3

import numpy as np
import sys
sys.path.append('../../../src/tools')
import saveParticlesVTK
from writeInitConfigFile import save_params_to_json

def compute_hydrostatic_density( ry, fluid_height, rho0, speed_of_sound ):
    hydrostaticPressure = rho0 * 9.81 * ( fluid_height - ry )
    hydrostaticDensity = ( ( hydrostaticPressure / ( speed_of_sound** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0
    #hydrostaticDensity = rho0 + hydrostaticPressure / speed_of_sound**2
    return hydrostaticDensity

def generate_dam_break_fluid_particles( setup, fluid_length, fluid_height ):
    fluid_rx = []; fluid_ry = []
    fluid_density = []
    dp = setup[ 'dp' ]
    fluid_lenght_n = round( fluid_length / dp )
    fluid_height_n = round( fluid_height / dp )
    rho0 = setup[ 'rho0' ]
    speed_of_sound = setup[ 'speedOfSound' ]

    for x in range( fluid_lenght_n ):
        for y in range( fluid_height_n ):
            fluid_rx.append( dp * ( x + 1 ) )
            fluid_ry.append( dp * ( y + 1 ) )
            fluid_density.append( compute_hydrostatic_density( fluid_ry[ -1 ], fluid_height, rho0, speed_of_sound ) )

    fluid_n = len( fluid_rx )
    fluid_r = np.array( ( fluid_rx, fluid_ry, np.zeros( fluid_n ) ), dtype=float ).T #!!
    fluid_v = np.zeros( ( fluid_n, 3 ) )
    fluid_rho = np.array( fluid_density, dtype=float )
    fluid_p = np.zeros( fluid_n )
    fluid_ptype = np.zeros( fluid_n )
    fluid_to_write = saveParticlesVTK.create_pointcloud_polydata( fluid_r, fluid_v, fluid_rho, fluid_p, fluid_ptype )
    saveParticlesVTK.save_polydata( fluid_to_write, "sources/dambreak_fluid.vtk" )

    # compute potential energy
    mass = rho0 * ( dp * dp )
    Epot0 = mass * 9.81 * np.sum( fluid_ry )
    print( f"Initial potential energy of fluid Epot0: {Epot0}" )
    setup[ "numberOfParticles" ] = fluid_n
    setup[ "numberOfAllocatedParticles" ] = fluid_n

def generate_dam_break_boundary_particles( setup, box_length, box_height, n_boundary_layers ):
    box_rx = []; box_ry = []
    ghost_rx = []; ghost_ry = []
    normal_x = []; normal_y = []
    dp = setup[ 'dp' ]
    box_length_n = round( box_length / dp )
    box_height_n = round( box_height / dp )
    rho0 = setup[ 'rho0' ]

    # left wall
    for layer in range( n_boundary_layers ):
        for y in range( box_height_n - 1 ):
            box_rx.append( 0. - layer * dp )
            box_ry.append( ( y + 1 ) * dp )
            ghost_rx.append( 0. + dp * ( layer + 1 ) )
            ghost_ry.append( ( y + 1 ) * dp )
            normal_x.append( 1. )
            normal_y.append( 0. )

    # bottom wall
    for layer in range( n_boundary_layers ):
        for x in range( box_length_n - n_boundary_layers + 1 ):
            box_rx.append( ( x + 1 ) * dp )
            box_ry.append( 0. - layer * dp )
            ghost_rx.append( ( x + 1 ) * dp )
            ghost_ry.append( 0. + dp * ( layer + 1 ) )
            normal_x.append( 0. )
            normal_y.append( 1. )

    x_last = box_rx[ -1 ] + dp #due to discretisation, we need to save last value of bottom wall

    # right wall
    for layer in range( n_boundary_layers ):
        for y in range( box_height_n - 1 ):
            box_rx.append( x_last + dp * layer )
            box_ry.append( ( y + 1 ) * dp )
            ghost_rx.append( x_last - dp * ( layer + 1 ) )
            ghost_ry.append( ( y + 1 ) * dp )
            normal_x.append( -1. )
            normal_y.append( 0. )

    # generate the corners
    def generate90degCorner( x, y, dirx, diry ):
      for layer in range( n_boundary_layers ):
        for k in range( n_boundary_layers ):
          box_rx.append( x + k * dp * dirx )
          box_ry.append( y + layer * dp * diry )
          ghost_rx.append( x + ( k + 1 ) * dp * dirx * ( -1 ) )
          ghost_ry.append( y + ( layer + 1 ) * dp * diry * ( -1 ) )
          drx = ghost_rx[ -1 ] - box_rx[ -1 ]
          dry = ghost_ry[ -1 ] - box_ry[ -1 ]
          n_norm = np.sqrt( drx**2 + dry**2 )
          normal_x.append( drx / n_norm )
          normal_y.append( dry / n_norm )

    generate90degCorner( 0, 0., -1, -1 )
    generate90degCorner( x_last, 0., +1, -1 )

    boundary_n = len( box_rx )
    boundary_r = np.array( ( box_rx, box_ry, np.zeros( boundary_n ) ), dtype=float ).T #!!
    boundary_ghostNodes = np.array( ( ghost_rx, ghost_ry, np.zeros( boundary_n ) ), dtype=float ).T #!!
    boundary_normals = np.array( ( normal_x, normal_y, np.zeros( boundary_n ) ), dtype=float ).T #!!
    boundary_v = np.zeros( ( boundary_n, 3 ) )
    boundary_rho = rho0 * np.ones( boundary_n )
    boundary_p = np.zeros( boundary_n )
    boundary_ptype = np.ones( boundary_n )
    box_to_write = saveParticlesVTK.create_pointcloud_polydata(
                    boundary_r,
                    boundary_v,
                    boundary_rho,
                    boundary_p,
                    boundary_ptype,
                    ghostNodes=boundary_ghostNodes,
                    normals=boundary_normals )
    saveParticlesVTK.save_polydata( box_to_write, "sources/dambreak_boundary.vtk" )

    setup[ "numberOfBoundaryParticles" ] = boundary_n
    setup[ "numberOfAllocatedBoundaryParticles" ] = boundary_n
    setup[ "domainOrigin-x" ] = min( box_rx )
    setup[ "domainOrigin-y" ] = min( box_ry )
    return max( box_rx ), max( box_ry )

def compute_domain_size( setup, domain_end_x, domain_end_y ):
    search_radius = setup[ "searchRadius" ]
    # Resize domain by one layer of cells
    eps = 1.005
    eps_sloshing = 1.2
    domain_origin_x = eps * ( setup[ "domainOrigin-x" ] - search_radius )
    domain_origin_y = eps * ( setup[ "domainOrigin-y" ] - search_radius )
    domain_end_x = eps * ( domain_end_x + search_radius )
    domain_end_y = eps_sloshing * ( domain_end_y + search_radius )
    domain_size_x = domain_end_x - domain_origin_x
    domain_size_y = domain_end_y - domain_origin_y

    setup[ "domainOrigin-x" ] = domain_origin_x
    setup[ "domainOrigin-y" ] = domain_origin_y
    setup[ "domainSize-x" ] = domain_size_x
    setup[ "domainSize-y" ] = domain_size_y

def write_config_header( compile_time_config ):
    with open( 'template/config_template.h', 'r' ) as file :
      config_file = file.read()

    config_file = config_file.replace( '#placeholderBoundaryConditionsType',  compile_time_config[ "bc_type" ] )
    config_file = config_file.replace( '#placeholderDiffusiveTerm', compile_time_config[ "diffusive_term" ] )
    config_file = config_file.replace( '#placeholderViscosTerm', compile_time_config[ "viscous_term" ] )

    with open( 'template/config.h', 'w' ) as file:
      file.write( config_file )

if __name__ == "__main__":
    import sys
    import argparse
    import os
    from pprint import pprint

    argparser = argparse.ArgumentParser(description="Dam break example initial condition and setup generator")
    g = argparser.add_argument_group("resolution parameters")
    g.add_argument("--dp", type=float, default=0.002, help="initial distance between particles")
    g.add_argument("--h-coef", type=float, default=2**0.5, help="smoothing length coefitient")
    g = argparser.add_argument_group("domain parameters")
    g.add_argument("--box-length", type=float, default=1.61, help="length of dam break box")
    g.add_argument("--box-height", type=float, default=1.0, help="height of dam break box")
    g.add_argument("--fluid-length", type=float, default=0.6, help="length of fluid block")
    g.add_argument("--fluid-height", type=float, default=0.3, help="height of fluid block")
    g = argparser.add_argument_group("simulation parameters")
    g.add_argument("--density", type=float, default=1000, help="referential density of the fluid")
    g.add_argument("--speed-of-sound", type=float, default=34.3, help="speed of sound")
    g.add_argument("--cfl", type=float, default=0.25, help="referential density of the fluid")
    g.add_argument("--bc-type", type=str, default="DBC", help="type of solid walls boundary conditions")
    g.add_argument("--diffusive-term", type=str, default="MolteniDiffusiveTerm", help="type of solid walls boundary conditions")
    g.add_argument("--viscous-term", type=str, default="ArtificialViscosity", help="type of solid walls boundary conditions")
    g.add_argument("--alpha", type=float, default=0.02, help="artificial vicosity parameter")
    g.add_argument("--dynamic-viscosity", type=float, default=0.001, help="dynamic viscosity")

    args = argparser.parse_args()

    # geometry inputs (not part of the solver config)
    box_length = args.box_length
    box_height = args.box_height
    fluid_length = args.fluid_length
    fluid_height = args.fluid_height
    h_coef = args.h_coef
    n_boundary_layers = int( np.ceil( 2 * h_coef * args.dp / args.dp ) )

    # compile-time parameters (config.h — C++ template type selection, not dumpable to JSON)
    compile_time_config = {
        "bc_type": args.bc_type,
        "diffusive_term": args.diffusive_term,
        "viscous_term": args.viscous_term,
    }

    # setup dict uses config-compatible keys so it can be dumped directly to config.jsonc
    dambreak_setup = {
        "case-name": "damBreak2D_WCSPH-DBC",
        "output-directory": "results",
        "fluid-particles": "sources/dambreak_fluid.vtk",
        "boundary-particles": "sources/dambreak_boundary.vtk",
        "final-time": 2.0,
        "snapshot-period": 0.05,
        "verbose-intensity": "with-snapshot",
        "numberOfParticles": 0,
        "numberOfAllocatedParticles": 0,
        "numberOfBoundaryParticles": 0,
        "numberOfAllocatedBoundaryParticles": 0,
        "searchRadius": 2 * h_coef * args.dp,
        "domainOrigin-x": 0,
        "domainOrigin-y": 0,
        "domainSize-x": 0,
        "domainSize-y": 0,
        "dp": args.dp,
        "h": h_coef * args.dp,
        "mass": args.density * ( args.dp * args.dp ),
        "rho0": args.density,
        "delta": 0.1,
        "alpha": args.alpha,
        "dynamicViscosity": args.dynamic_viscosity,
        "speedOfSound": args.speed_of_sound,
        "eps": 0.001,
        "external-force-x": 0,
        "external-force-y": -9.81,
        "initial-time-step": args.cfl * ( h_coef * args.dp ) / args.speed_of_sound,
        "CFL": args.cfl,
    }

    # check initial settings
    if fluid_length > box_length or fluid_height > box_length:
        sys.stderr.write( "Size of the fluid block must be smaller than the size of the box." )

    # create necessary folders
    resultsPath = r'./results'
    if not os.path.exists( resultsPath ):
        os.makedirs( resultsPath )

    sourcesPath = r'./sources'
    if not os.path.exists( sourcesPath ):
        os.makedirs( sourcesPath )

    # generate particles (updates numberOfParticles, numberOfBoundaryParticles, domainOrigin in setup)
    generate_dam_break_fluid_particles( dambreak_setup, fluid_length, fluid_height )
    domain_end_x, domain_end_y = generate_dam_break_boundary_particles(
        dambreak_setup, box_length, box_height, n_boundary_layers )

    # compute domain size (updates domainOrigin, domainSize in setup)
    compute_domain_size( dambreak_setup, domain_end_x, domain_end_y )

    print( "Complete example setup:" )
    pprint( dambreak_setup )

    # dump setup directly to config.jsonc
    save_params_to_json( dambreak_setup, 'sources/config.jsonc' )

    # write config header (compile-time C++ types)
    write_config_header( compile_time_config )
