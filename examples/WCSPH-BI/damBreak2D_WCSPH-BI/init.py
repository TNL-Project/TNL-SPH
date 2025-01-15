#! /usr/bin/env python3

import numpy as np
import sys
sys.path.append('../../../src/tools')
import saveParticlesVTK
import math

def compute_hydrostatic_density( ry, fluid_height, rho0, speed_of_sound ):
    hydrostaticPressure = rho0 * 9.81 * ( fluid_height - ry )
    #hydrostaticDensity = ( ( hydrostaticPressure / ( speed_of_sound** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0
    hydrostaticDensity = rho0 + hydrostaticPressure / speed_of_sound**2
    return hydrostaticDensity

def compute_hydrostatic_density_with_profile( rx, ry, fluid_lenght, fluid_height, rho0, speed_of_sound ):
    hydrostaticPressure = rho0 * 9.81 * ( fluid_height - ry ) * math.cos(0.5 * math.pi * ( rx ) / fluid_lenght )
    #hydrostaticDensity = ( ( hydrostaticPressure / ( speed_of_sound** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0
    hydrostaticDensity = rho0 + hydrostaticPressure / speed_of_sound**2
    return hydrostaticDensity

def generate_dam_break_fluid_particles( setup ):
    fluid_rx = []; fluid_ry = []
    fluid_density = []
    dp = setup[ 'dp' ]
    fluid_length = setup[ 'fluid_length' ];
    fluid_height = setup[ 'fluid_height' ]
    fluid_lenght_n = round( fluid_length / dp )
    fluid_height_n = round( fluid_height / dp )
    rho0 = setup[ 'density' ]
    speed_of_sound = setup[ 'speed_of_sound' ]

    for x in range( fluid_lenght_n ):
        for y in range( fluid_height_n ):
            fluid_rx.append( dp * ( x + 1 ) - dp / 2 ) #Ensure that the interface is at 0, shift the block
            fluid_ry.append( dp * ( y + 1 ) - dp / 2 ) #Ensure that the interface is at 0, shift the block
            hydrostatic_density = compute_hydrostatic_density_with_profile(
                                    fluid_rx[ -1 ],
                                    fluid_ry[ -1 ],
                                    fluid_length,
                                    fluid_height,
                                    rho0,
                                    speed_of_sound )
            fluid_density.append( hydrostatic_density )

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
    setup[ "fluid_n" ] = fluid_n

def generate_dam_break_boundary_particles_light( setup ):
    box_rx = []; box_ry = []
    normal_x = []; normal_y = []
    dp = setup[ 'dp' ]
    box_length_n = round( setup[ 'box_length' ] / dp )
    box_height_n = round( setup[ 'box_height' ] / dp )
    rho0 = setup[ 'density' ]

    # left wall
    for y in range( box_height_n ):
        box_rx.append( 0. )
        box_ry.append( ( y + 1 ) * dp - dp / 2)
        normal_x.append( 1. )
        normal_y.append( 0. )

    # bottom wall
    for x in range( box_length_n ):
        box_rx.append( ( x + 1 ) * dp - dp / 2 )
        box_ry.append( 0. )
        normal_x.append( 0. )
        normal_y.append( 1. )

    x_last = box_rx[ -1 ] + dp / 2 # due to discretisation, we need to save last value of bottom wall

    # right wall
    for y in range( box_height_n ):
        box_rx.append( x_last )
        box_ry.append( ( y + 1 ) * dp - dp / 2 )
        normal_x.append( -1. )
        normal_y.append( 0. )

    y_last = box_ry[ -1 ] + dp / 2 # due to discretisation, we need to save last value of bottom wall

    # top wall
    for x in range( box_length_n ):
        box_rx.append( ( x + 1 ) * dp - dp / 2 )
        box_ry.append( y_last )
        normal_x.append( 0. )
        normal_y.append( -1. )

    boundary_n = len( box_rx )
    boundary_r = np.array( ( box_rx, box_ry, np.zeros( boundary_n ) ), dtype=float ).T #!!
    boundary_normal = np.array( ( normal_x, normal_y, np.zeros( boundary_n ) ), dtype=float ).T #!!
    boundary_v = np.zeros( ( boundary_n, 3 ) )
    boundary_rho = rho0 * np.ones( boundary_n )
    boundary_p = np.zeros( boundary_n )
    boundary_elemetnSize = dp * np.ones( boundary_n )
    boundary_ptype = np.ones( boundary_n )
    box_to_write = saveParticlesVTK.create_pointcloud_polydata(
                    boundary_r,
                    boundary_v,
                    boundary_rho,
                    boundary_p,
                    boundary_ptype,
                    normals=boundary_normal,
                    elementSize=boundary_elemetnSize )
    saveParticlesVTK.save_polydata( box_to_write, "sources/dambreak_boundary.vtk" )

    setup[ "boundary_n" ] = boundary_n
    setup[ "domain_origin_x" ] = min( box_rx )
    setup[ "domain_origin_y" ]  = min( box_ry )
    setup[ "domain_end_x" ] = max( box_rx )
    setup[ "domain_end_y" ] = max( box_ry )

def compute_domain_size( setup ):
    search_radius = setup[ "search_radius" ]
    # Resize domain by one layer of cells
    eps = 1.005
    eps_sloshing = 1.2
    domain_origin_x = eps * ( setup[ "domain_origin_x" ] - search_radius )
    domain_origin_y = eps * ( setup[ "domain_origin_y" ] - search_radius )
    domain_end_x = eps * ( setup[ "domain_end_x" ] + search_radius )
    domain_end_y = eps_sloshing * ( setup[ "domain_end_y" ] + search_radius )
    domain_size_x = domain_end_x - domain_origin_x
    domain_size_y = domain_end_y - domain_origin_y

    extra_parameters = {
        "domain_origin_x" : domain_origin_x,
        "domain_origin_y" : domain_origin_y,
        "domain_size_x" : domain_size_x,
        "domain_size_y" : domain_size_y,
    }
    setup.update( extra_parameters )

def write_simulation_params( setup ):
    # write parameters to config file
    with open( 'template/config_template.ini', 'r' ) as file :
      config_file = file.read()

    config_file = config_file.replace( 'placeholderSearchRadius', f'{ setup[ "search_radius" ] }' )
    config_file = config_file.replace( 'placeholderDomainOrigin-x', f'{setup[ "domain_origin_x" ]:.5f}' )
    config_file = config_file.replace( 'placeholderDomainOrigin-y', f'{setup[ "domain_origin_y" ]:.5f}' )
    config_file = config_file.replace( 'placeholderDomainSize-x', f'{setup[ "domain_size_x" ]:.5f}' )
    config_file = config_file.replace( 'placeholderDomainSize-y', f'{setup[ "domain_size_y" ]:.5f}' )

    config_file = config_file.replace( 'placeholderInitParticleDistance', f'{ setup[ "dp" ] }' )
    config_file = config_file.replace( 'placeholderSmoothingLength', f'{ setup[ "smoothing_length" ] }' )
    config_file = config_file.replace( 'placeholderMass', f'{ setup[ "particle_mass" ] }' )
    config_file = config_file.replace( 'placeholderBoundaryElementSize', f'{ setup[ "dp" ] }' )
    config_file = config_file.replace( 'placeholderSpeedOfSound', f'{ setup[ "speed_of_sound" ] }' )
    config_file = config_file.replace( 'placeholderDensity', f'{ setup[ "density" ] }' )
    config_file = config_file.replace( 'placeholderTimeStep', f'{ setup[ "time_step" ] }' )
    config_file = config_file.replace( 'placeholderCFL', f'{ setup[ "cfl" ] }' )
    config_file = config_file.replace( 'placeholderFluidParticles', f'{ setup[ "fluid_n" ] }' )
    config_file = config_file.replace( 'placeholderAllocatedFluidParticles', f'{ setup[ "fluid_n" ] }' )
    config_file = config_file.replace( 'placeholderBoundaryParticles', f'{ setup[ "boundary_n" ] }' )
    config_file = config_file.replace( 'placeholderAllocatedBoundaryParticles', f'{ setup[ "boundary_n" ] }' )

    with open( 'sources/config.ini', 'w' ) as file:
      file.write( config_file )

def configure_and_write_measuretool_parameters():
    # write parameters to config file
    with open( 'template/config-measuretool_template.ini', 'r' ) as file :
      config_file = file.read()
    with open( 'sources/config-measuretool.ini', 'w' ) as file:
      file.write( config_file )

if __name__ == "__main__":
    import sys
    import argparse
    import os
    from pprint import pprint

    argparser = argparse.ArgumentParser(description="Heat equation example initial condition generator")
    g = argparser.add_argument_group("resolution parameters")
    g.add_argument("--dp", type=float, default=0.002, help="initial distance between particles")
    g.add_argument("--h-coef", type=float, default=2, help="smoothing length coefitient")
    g = argparser.add_argument_group("domain parameters")
    g.add_argument("--box-length", type=float, default=1.61, help="length of dam break box")
    g.add_argument("--box-height", type=float, default=0.6, help="height of dam break box")
    g.add_argument("--fluid-length", type=float, default=0.6, help="length of fluid block")
    g.add_argument("--fluid-height", type=float, default=0.3, help="height of fluid block")
    g.add_argument("--n-boundary-layers", type=int, default=3, help="number of boundary layers")
    g = argparser.add_argument_group("simulation parameters")
    g.add_argument("--density", type=float, default=997, help="referential density of the fluid")
    g.add_argument("--speed-of-sound", type=float, default=17.155174146594955, help="speed of sound")
    g.add_argument("--cfl", type=float, default=0.10, help="referential density of the fluid")
    #g = argparser.add_argument_group("control parameters")
    #g.add_argument("--example-dir", type=Path, default=1000, help="referential density of the fluid")

    args = argparser.parse_args()

    dambreak_setup = {
        # geometri proportions
        "box_height" : args.box_height,
        "box_length" : args.box_length,
        "fluid_height" : args.fluid_height,
        "fluid_length" : args.fluid_length,
        # general parameteres
        "dp" : args.dp,
        "h_coef" : args.h_coef,
        "density" : args.density,
        "speed_of_sound" : args.speed_of_sound,
        "cfl" : args.cfl,
        "particle_mass" : args.density * ( args.dp * args.dp ),
        "boundary_element_size" : args.dp * args.dp,
        "smoothing_length" : args.h_coef * args.dp,
        "search_radius" : 2 * args.h_coef * args.dp,
        "time_step" : args.cfl * ( args.h_coef * args.dp ) / args.speed_of_sound
    }

    # check initial settings
    if args.fluid_length > args.box_length or args.fluid_height > args.box_length:
        sys.stderr.write( "Size of the fluid block must be smaller than the size of the box." )

    # create necessary folders
    resultsPath = r'./results'
    if not os.path.exists( resultsPath ):
        os.makedirs( resultsPath )

    sourcesPath = r'./sources'
    if not os.path.exists( sourcesPath ):
        os.makedirs( sourcesPath )

    # generate particles
    generate_dam_break_fluid_particles( dambreak_setup )
    generate_dam_break_boundary_particles_light( dambreak_setup )

    # setup parameters
    compute_domain_size( dambreak_setup )

    print( "Complete example setup:" )
    pprint( dambreak_setup )
    # write simulation params
    write_simulation_params( dambreak_setup )
