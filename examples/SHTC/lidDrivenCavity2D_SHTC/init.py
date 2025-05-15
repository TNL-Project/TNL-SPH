#! /usr/bin/env python3

import numpy as np
import sys
sys.path.append('../../../src/tools')
import saveParticlesVTK

def generate_cavity_fluid_particles( setup ):
    fluid_rx = []; fluid_ry = []
    boundary_rx = []; boundary_ry = []
    boundary_vx = []; boundary_vy = []
    boundary_nx = []; boundary_ny = []
    dp = setup[ 'dp' ]
    dpx = ( 4 / 3 )**( 1 / 4 ) * setup[ 'dp' ]
    dpy = ( 3 / 4 )**( 1 / 4 ) * setup[ 'dp' ]
    cavity_length = setup[ 'cavity_length' ]
    cavity_height = setup[ 'cavity_height' ]
    cavity_lenght_n = round( cavity_length / dpx )
    cavity_height_n = round( cavity_height / dpy )
    rho0 = setup[ 'density' ]

    for x in range( -3, cavity_lenght_n + 4 ):
        for y in range( -3, cavity_height_n + 4 ):

            if y % 2 == 0 : rx = dpx * ( x + 0.0 )
            else: rx = dpx * ( x + 0.5 )
            ry = dpy * ( y + 0 )

            if rx < 0 or rx > cavity_length or ry < 0 or ry > cavity_height:
                boundary_rx.append( rx )
                boundary_ry.append( ry )
                if ry > cavity_height:
                    boundary_vx.append( 1 )
                    boundary_vy.append( 0 )
                else:
                    boundary_vx.append( 0 )
                    boundary_vy.append( 0 )

            else:
                fluid_rx.append( rx )
                fluid_ry.append( ry )

    fluid_n = len( fluid_rx )
    fluid_r = np.array( ( fluid_rx, fluid_ry, np.zeros( fluid_n ) ), dtype=float ).T #!!
    fluid_v = np.zeros( ( fluid_n, 3 ) )
    fluid_rho = rho0 * np.ones( fluid_n )
    fluid_p = np.zeros( fluid_n )
    fluid_ptype = np.zeros( fluid_n )
    fluid_to_write = saveParticlesVTK.create_pointcloud_polydata(
                      fluid_r,
                      fluid_v,
                      fluid_rho,
                      fluid_p,
                      fluid_ptype )
    saveParticlesVTK.save_polydata( fluid_to_write, "sources/lidDrivenCavity_fluid.vtk" )

    boundary_n = len( boundary_rx )
    boundary_r = np.array( ( boundary_rx, boundary_ry, np.zeros( boundary_n ) ), dtype=float ).T #!!
    boundary_v = np.array( ( boundary_vx, boundary_vy, np.zeros( boundary_n ) ), dtype=float ).T #
    boundary_rho = rho0 * np.ones( boundary_n )
    boundary_p = np.zeros( boundary_n )
    boundary_normals = np.zeros( ( boundary_n, 3 ) ) #FIXME: Boundary are not generated since we don't need them!
    boundary_ptype = np.ones( boundary_n )
    box_to_write = saveParticlesVTK.create_pointcloud_polydata(
                    boundary_r,
                    boundary_v,
                    boundary_rho,
                    boundary_p,
                    boundary_ptype,
                    normals = boundary_normals )
    saveParticlesVTK.save_polydata( box_to_write, "sources/lidDrivenCavity_boundary.vtk" )

    # compute potential energy
    mass = rho0 * ( dp * dp )
    Epot0 = mass * 9.81 * np.sum( fluid_ry )
    print( f"Initial potential energy of fluid Epot0: {Epot0}" )
    setup[ "fluid_n" ] = fluid_n

    setup[ "boundary_n" ] = boundary_n
    setup[ "domain_origin_x" ] = min( boundary_rx )
    setup[ "domain_origin_y" ]  = min( boundary_ry )
    setup[ "domain_end_x" ] = max( boundary_rx )
    setup[ "domain_end_y" ] = max( boundary_ry )

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
    config_file = config_file.replace( 'placeholderSpeedOfSoundBulk', f'{ setup[ "speed_of_sound_bulk" ] }' )
    config_file = config_file.replace( 'placeholderSpeedOfSoundShear', f'{ setup[ "speed_of_sound_shear" ] }' )
    config_file = config_file.replace( 'placeholderTau', f'{ setup[ "tau" ] }' )
    config_file = config_file.replace( 'placeholderDensity', f'{ setup[ "density" ] }' )
    config_file = config_file.replace( 'placeholderTimeStep', f'{setup[ "time_step" ] }' )
    config_file = config_file.replace( 'placeholderCFL', f'{ setup[ "cfl" ] }' )
    config_file = config_file.replace( 'placeholderDynamicVicosity', f'{ setup[ "dynamic_viscosity" ] }' )
    config_file = config_file.replace( 'placeholderFluidParticles', f'{ setup[ "fluid_n" ] }' )
    config_file = config_file.replace( 'placeholderAllocatedFluidParticles', f'{ setup[ "fluid_n" ] }' )
    config_file = config_file.replace( 'placeholderBoundaryParticles', f'{ setup[ "boundary_n" ] }' )
    config_file = config_file.replace( 'placeholderAllocatedBoundaryParticles', f'{ setup[ "boundary_n" ] }' )

    with open( 'sources/config.ini', 'w' ) as file:
      file.write( config_file )

def configure_and_write_measuretool_parameters( setup ):
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

    argparser = argparse.ArgumentParser(description="Lid driven example initial condition and setup generator")
    g = argparser.add_argument_group("resolution parameters")
    g.add_argument("--dp", type=float, default=0.01, help="initial distance between particles")
    g.add_argument("--h-coef", type=float, default=2.4, help="smoothing length coefitient")
    g = argparser.add_argument_group("domain parameters")
    g.add_argument("--cavity-length", type=float, default=1.0, help="length of the cavity")
    g.add_argument("--cavity-height", type=float, default=1.0, help="height of the cavity")
    g.add_argument("--top-wall-velocity", type=float, default=1.0, help="height of dam break box")
    g.add_argument("--Re", type=float, default=100, help="cavity Reynlods number")
    g = argparser.add_argument_group("simulation parameters")
    g.add_argument("--density", type=float, default=1, help="referential density of the fluid")
    g.add_argument("--speed-of-sound-bulk", type=float, default=20, help="bulk speed of sound")
    g.add_argument("--speed-of-sound-shear", type=float, default=20, help="shear speed of sound")
    g.add_argument("--cfl", type=float, default=0.05, help="referential density of the fluid")

    args = argparser.parse_args()

    dambreak_setup = {
        # geometri proportions
        "cavity_length" : args.cavity_length,
        "cavity_height" : args.cavity_height,
        "top_wall_velocity" : args.top_wall_velocity,
        # general parameteres
        "dp" : args.dp,
        "h_coef" : args.h_coef,
        "n_boundary_layers" : int( np.ceil( 2 * args.h_coef * args.dp  / args.dp ) ),
        "density" : args.density,
        "speed_of_sound_bulk" : args.speed_of_sound_bulk,
        "speed_of_sound_shear" : args.speed_of_sound_shear,
        "cfl" : args.cfl,
        "particle_mass" : args.density * ( args.dp * args.dp ),
        "smoothing_length" : args.h_coef * args.dp,
        #"search_radius" : 2 * args.h_coef * args.dp,
        "search_radius" : args.h_coef * args.dp,
        "time_step" : args.cfl * ( args.h_coef * args.dp ) / args.speed_of_sound_bulk,
        "tau" : 6 * args.top_wall_velocity * args.cavity_length / ( args.Re * args.speed_of_sound_shear**2 ),
        "dynamic_viscosity" : args.top_wall_velocity * args.cavity_length / args.Re
    }

    # create necessary folders
    resultsPath = r'./results'
    if not os.path.exists( resultsPath ):
        os.makedirs( resultsPath )

    sourcesPath = r'./sources'
    if not os.path.exists( sourcesPath ):
        os.makedirs( sourcesPath )

    # generate particles
    generate_cavity_fluid_particles( dambreak_setup )

    # setup parameters
    compute_domain_size( dambreak_setup )

    print( "Complete example setup:" )
    pprint( dambreak_setup )
    # write simulation params
    write_simulation_params( dambreak_setup )
    configure_and_write_measuretool_parameters( dambreak_setup )
