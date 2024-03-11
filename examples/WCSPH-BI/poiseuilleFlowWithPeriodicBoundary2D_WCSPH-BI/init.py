#! /usr/bin/env python3

import numpy as np
import sys
sys.path.append('../../../src/tools')
import saveParticlesVTK
import domainGrid

def generate_channel_fluid_particles( setup ):
    dp = setup[ "dp" ]
    fluid_rx = []; fluid_ry = []
    fluid_vx = []; fluid_vy = []
    fluid_lenght_n = round( ( setup[ "channel_length" ] - dp ) / dp )
    fluid_height_n = round( ( setup[ "channel_height" ] - dp ) / dp )

    for x in range( fluid_lenght_n ):
        for y in range( fluid_height_n ):
            fluid_rx.append( dp * ( x + 1 ) )
            fluid_ry.append( dp * ( y + 1 ) )
            fluid_vx.append( setup[ "periodicityLeft_velocity_x" ] )
            fluid_vy.append( setup[ "periodicityLeft_velocity_y" ] )

    fluid_n = len( fluid_rx )
    fluid_r = np.array( ( fluid_rx, fluid_ry, np.zeros( fluid_n ) ), dtype=float ).T #!!
    fluid_v = np.array( ( fluid_vx, fluid_vy, np.zeros( fluid_n ) ), dtype=float ).T #!!
    fluid_rho = setup[ "density" ] * np.ones( fluid_n )
    fluid_p = np.zeros( fluid_n )
    fluid_ptype = np.zeros( fluid_n )
    fluid_to_write = saveParticlesVTK.create_pointcloud_polydata( fluid_r, fluid_v, fluid_rho, fluid_p, fluid_ptype )
    saveParticlesVTK.save_polydata( fluid_to_write, "sources/openchannel_fluid.vtk" )
    setup[ "fluid_n" ] = fluid_n

def generate_channel_boundary_particles( setup ):
    dp = setup[ "dp" ]
    n_boundary_layers = setup[ "n_boundary_layers" ]

    box_rx = []; box_ry = []
    normal_x = []; normal_y = []
    box_length_n = round( setup[ "channel_length" ] / dp )
    box_height_n = round( setup[ "channel_height" ] / dp )

    # bottom wall
    for layer in range( n_boundary_layers ):
        #FIXME: Im currently no able to do periodicity with boundary set, so the wall is artificialy extended
        #       which is dane by setting n_boundary_layers fixed to 3
        #       for x in range( box_length_n + ( n_boundary_layers - 1 ) * 2 + 1):
        #for x in range( box_length_n + ( 3 - 1 ) * 2 + 1):
        for x in range( box_length_n - 1 ):
            box_rx.append( ( x + 1 ) * dp )
            box_ry.append( 0. - layer * dp )
            normal_x.append( 0. )
            normal_y.append( 1. )

    # top wall
    for layer in range( n_boundary_layers ):
        #FIXME: Im currently no able to do periodicity with boundary set, so the wall is artificialy extended
        #       which is dane by setting n_boundary_layers fixed to 3
        #       for x in range( box_length_n + ( n_boundary_layers - 1 ) * 2 + 1):
        #for x in range( box_length_n + ( 3 - 1 ) * 2 + 1):
        for x in range( box_length_n - 1 ):
            box_rx.append( ( x + 1 ) * dp )
            box_ry.append( ( setup[ "channel_height" ] ) + layer * dp )
            normal_x.append( 0. )
            normal_y.append( -1. )

    boundary_n = len( box_rx )
    boundary_r = np.array( ( box_rx, box_ry, np.zeros( boundary_n ) ), dtype=float ).T #!!
    boundary_normal = np.array( ( normal_x, normal_y, np.zeros( boundary_n ) ), dtype=float ).T #!!
    boundary_v = np.zeros( ( boundary_n, 3 ) )
    boundary_rho = setup[ "density" ] * np.ones( boundary_n )
    boundary_p = np.zeros( boundary_n )
    boundary_ptype = np.ones( boundary_n )
    box_to_write = saveParticlesVTK.create_pointcloud_polydata( boundary_r, boundary_v, boundary_rho, boundary_p, boundary_ptype,
                                                                boundary_normal )
    saveParticlesVTK.save_polydata( box_to_write, "sources/openchannel_boundary.vtk" )

    setup[ "boundary_n" ] = boundary_n
    setup[ "domain_origin_x" ] = min( box_rx )
    setup[ "domain_origin_y" ]  = min( box_ry )
    setup[ "domain_end_x" ] = max( box_rx )
    setup[ "domain_end_y" ] = max( box_ry )

def generate_channel_open_boundary_particles( setup, prefix ):
    dp = setup[ "dp" ]
    patch_rx = []; patch_ry = []
    patch_vx = []; patch_vy = []

    patch_width_n = setup[ prefix + "_layers" ]
    patch_height_n = round( setup[ prefix + "_height" ] / dp  )

    for x in range( patch_width_n ):
        for y in range( patch_height_n ):
            patch_rx.append( setup[ prefix + "_position_x" ] - setup[ prefix + "_orientation_x" ] * dp *  x  )
            patch_ry.append( setup[ prefix + "_position_y" ] + dp * ( y ) )
            patch_vx.append( setup[ prefix + "_velocity_x" ] )
            patch_vy.append( setup[ prefix + "_velocity_y" ] )

    patch_n = len( patch_rx )
    r = np.array( ( patch_rx, patch_ry, np.zeros( patch_n ) ), dtype=float ).T #!!
    v = np.array( ( patch_vx, patch_vy, np.zeros( patch_n ) ), dtype=float ).T #!!
    rho = setup[ "density" ] * np.ones( patch_n, dtype=float )
    p = np.zeros( patch_n )
    ptype = np.ones( patch_n )

    inletToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
    saveParticlesVTK.save_polydata( inletToWrite, "sources/openchannel_" + prefix + ".vtk" )

    setup[ prefix + "_n" ] = patch_n

def compute_domain_size( setup ):
    search_radius = setup[ "search_radius" ]

    # Resize domain by one layer of cells
    # For BI, we use overlap 2 * search_radius just to be sure we have the additional layer of empty cells
    eps = 1.005
    domain_origin_x = eps * ( setup[ "domain_origin_x" ] -  2 * search_radius )
    domain_origin_y = eps * ( setup[ "domain_origin_y" ] - 2 * search_radius )
    domain_end_x = eps * ( setup[ "domain_end_x" ] + 2 * search_radius )
    domain_end_y = eps * ( setup[ "domain_end_y" ] + 2 * search_radius )
    domain_size_x = domain_end_x - domain_origin_x
    domain_size_y = domain_end_y - domain_origin_y

    extra_parameters = {
        "domain_origin_x" : domain_origin_x,
        "domain_origin_y" : domain_origin_y,
        "domain_size_x" : domain_size_x,
        "domain_size_y" : domain_size_y
    }
    setup.update( extra_parameters )

def write_simulation_params( setup ):

    # write parameters to config file
    with open( 'template/config_template.ini', 'r' ) as file :
      config_file = file.read()

    config_file = config_file.replace( 'placeholderSearchRadius', str( round( setup[ "search_radius" ], 7 ) ) )
    config_file = config_file.replace( 'placeholderDomainOrigin-x', str( round( setup[ "domain_origin_x" ], 5 ) ) )
    config_file = config_file.replace( 'placeholderDomainOrigin-y', str( round( setup[ "domain_origin_y" ], 5 ) ) )
    config_file = config_file.replace( 'placeholderDomainSize-x', str( round( setup[ "domain_size_x" ], 5  ) ) )
    config_file = config_file.replace( 'placeholderDomainSize-y', str( round( setup[ "domain_size_y" ], 5  ) ) )

    config_file = config_file.replace( 'placeholderInitParticleDistance', str( setup[ "dp" ] ) )
    config_file = config_file.replace( 'placeholderSmoothingLength', str( round( setup[ "smoothing_lenght" ], 7 ) ) )
    config_file = config_file.replace( 'placeholderMass', str( setup[ "particle_mass" ] ) )
    config_file = config_file.replace( 'placeholderSpeedOfSound', str( setup[ "speed_of_sound" ] ) )
    config_file = config_file.replace( 'placeholderDensity', str( setup[ "density" ] ) )
    config_file = config_file.replace( 'placeholderTimeStep', str( setup[ "time_step" ] ) )
    config_file = config_file.replace( 'placeholderFluidParticles', str( setup[ "fluid_n" ] ) )
    config_file = config_file.replace( 'placeholderAllocatedFluidParticles', str( setup[ "fluid_n" ] ) )
    config_file = config_file.replace( 'placeholderBoundaryParticles', str( setup[ "boundary_n" ] ) )
    config_file = config_file.replace( 'placeholderAllocatedBoundaryParticles', str( setup[ "boundary_n" ] ) )
    config_file = config_file.replace( 'placeholderPeriodicityLeftParticles', str( setup[ "periodicityLeft_n" ] ) )
    config_file = config_file.replace( 'placeholderPeriodicityLeftAllocatedParticles', str( 2 * setup[ "periodicityLeft_n" ] ) )
    config_file = config_file.replace( 'placeholderPeriodicityRightParticles', str( setup[ "periodicityRight_n" ] ) )
    config_file = config_file.replace( 'placeholderPeriodicityRightAllocatedParticles', str( 2 * setup[ "periodicityRight_n" ] ) )

    config_file = config_file.replace( 'placeholderInletOrientation_x', str( setup[ "periodicityLeft_orientation_x" ] ) )
    config_file = config_file.replace( 'placeholderInletOrientation_y', str( setup[ "periodicityLeft_orientation_y" ] ) )
    config_file = config_file.replace( 'placeholderInletVelocity_x', str( setup[ "periodicityLeft_velocity_x" ] ) )
    config_file = config_file.replace( 'placeholderInletVelocity_y', str( setup[ "periodicityLeft_velocity_y" ] ) )
    config_file = config_file.replace( 'placeholderInletPosition1_x', str( setup[ "periodicityLeft_position_x" ]  + setup[ "dp" ] / 2 ) ) #TODO
    config_file = config_file.replace( 'placeholderInletPosition1_y', str( setup[ "periodicityLeft_position_y" ] ) )
    config_file = config_file.replace( 'placeholderInletPosition2_x', str( setup[ "periodicityLeft_position_x" ]  + setup[ "dp" ] / 2 ) ) #TODO
    config_file = config_file.replace( 'placeholderInletPosition2_y', str( setup[ "periodicityLeft_position_y" ] + setup[ "periodicityLeft_height" ] - setup[ "dp" ] / 2 ) ) #TODO
    config_file = config_file.replace( 'placeholderInletDensity', str( setup[ "density" ] ) )
    config_file = config_file.replace( 'placeholderInletWidth_x', str( round( setup[ "periodicityLeft_width" ], 7 ) ) )
    config_file = config_file.replace( 'placeholderInletWidth_y', str( 0. ) )
    config_file = config_file.replace( 'placeholderInletShiftVector_x', str( round( setup[ "periodicityLeft_shift_vector_x" ], 7 ) ) )
    config_file = config_file.replace( 'placeholderInletShiftVector_y', str( 0. ) )

    config_file = config_file.replace( 'placeholderOutletOrientation_x', str( setup[ "periodicityRight_orientation_x" ] ) )
    config_file = config_file.replace( 'placeholderOutletOrientation_y', str( setup[ "periodicityRight_orientation_y" ] ) )
    config_file = config_file.replace( 'placeholderOutletVelocity_x', str( setup[ "periodicityRight_velocity_x" ] ) )
    config_file = config_file.replace( 'placeholderOutletVelocity_y', str( setup[ "periodicityRight_velocity_y" ] ) )
    config_file = config_file.replace( 'placeholderOutletPosition1_x', str( setup[ "periodicityRight_position_x" ] - setup[ "dp" ] / 2 ) ) #TODO
    config_file = config_file.replace( 'placeholderOutletPosition1_y', str( setup[ "periodicityRight_position_y" ] ) )
    config_file = config_file.replace( 'placeholderOutletPosition2_x', str( setup[ "periodicityRight_position_x" ] - setup[ "dp" ] / 2 ) ) #FIXME
    config_file = config_file.replace( 'placeholderOutletPosition2_y', str( setup[ "periodicityRight_position_y" ] + setup[ "periodicityRight_height" ] - setup[ "dp" ] / 2  ) )
    config_file = config_file.replace( 'placeholderOutletDensity', str( setup[ "density" ] ) )
    config_file = config_file.replace( 'placeholderOutletWidth_x', str( round( setup[ "periodicityRight_width" ], 7 ) ) )
    config_file = config_file.replace( 'placeholderOutletWidth_y', str( 0. ) )
    config_file = config_file.replace( 'placeholderOutletShiftVector_x', str( round( setup[ "periodicityRight_shift_vector_x" ], 7 ) ) )
    config_file = config_file.replace( 'placeholderOutletShiftVector_y', str( 0. ) )

    with open( 'sources/config.ini', 'w' ) as file:
      file.write( config_file )

def write_domain_background_grid( setup ):
    #TODO: Rename the DomainGrid function
    search_radius = setup[ "search_radius" ]
    grid_size_x = round( setup[ "domain_size_x" ] / search_radius )
    grid_size_y = round( setup[ "domain_size_y" ] / search_radius )
    grid_origin_x = setup[ "domain_origin_x" ]
    grid_origin_y = setup[ "domain_origin_y" ]
    grid_sectors = np.zeros( grid_size_x * grid_size_y )

    domainGrid.DomainGrid( grid_size_x,
                           grid_size_y,
                           0,
                           grid_origin_x,
                           grid_origin_y,
                           0,
                           grid_sectors,
                           search_radius,
                           "sources/openchannel_grid.vtk" )

if __name__ == "__main__":
    import sys
    import argparse
    import os

    argparser = argparse.ArgumentParser(description="Periodic channel example initial condition generator")
    g = argparser.add_argument_group("resolution parameters")
    g.add_argument("--dp", type=float, default=0.002, help="initial distance between particles")
    g.add_argument("--h-coef", type=float, default=2**0.5, help="smoothing length coefitient")
    g = argparser.add_argument_group("domain parameters")
    g.add_argument("--channel-length", type=float, default=0.1, help="length of fluid block")
    g.add_argument("--channel-height", type=float, default=0.1, help="height of fluid block")
    g.add_argument("--n-boundary-layers", type=int, default=3, help="number of boundary buffer layers")
    g.add_argument("--v-init", type=int, default=0.0, help="initial velocity")
    g = argparser.add_argument_group("simulation parameters")
    g.add_argument("--density", type=float, default=1000, help="referential density of the fluid")
    g.add_argument("--speed-of-sound", type=float, default=34.3, help="speed of sound")
    g.add_argument("--cfl", type=float, default=0.25, help="referential density of the fluid")

    args = argparser.parse_args()

    openchannel_setup = {
        # general parameteres
        "dp" : args.dp,
        "h_coef" : args.h_coef,
        "n_boundary_layers" : 1,
        "density" : args.density,
        "speed_of_sound" : args.speed_of_sound,
        "cfl" : args.cfl,
        "particle_mass" : args.density * ( args.dp * args.dp ),
        "smoothing_lenght" : args.h_coef * args.dp,
        "search_radius" :  2 * args.h_coef * args.dp,
        "time_step" : args.cfl * ( args.h_coef * args.dp ) / args.speed_of_sound,
        # geometric size
        "channel_length" : args.channel_length,
        "channel_height" : args.channel_height,
        # left periodic boundary condition
        "periodicityLeft_position_x" : 0.,
        "periodicityLeft_position_y" : 0. + args.dp,
        "periodicityLeft_orientation_x" : 1.,
        "periodicityLeft_orientation_y" : 0.,
        "periodicityLeft_layers" : args.n_boundary_layers + 1,
        "periodicityLeft_height" : args.channel_height - args.dp,
         #"periodicityLeft_width" : args.n_boundary_layers * args.dp,
        "periodicityLeft_width" : 1.2 * 2 * args.h_coef * args.dp,
        "periodicityLeft_velocity_x" : args.v_init,
        "periodicityLeft_velocity_y" : 0.,
        "periodicityLeft_reference_point_z" : 0. - 1. * args.n_boundary_layers * args.dp, #TODO: Make this more clear
        "periodicityLeft_reference_point_x" : args.dp,
        "periodicityLeft_shift_vector_x" : args.channel_length - 2 * ( args.dp / 2 ),
        "periodicityLeft_shift_vector_y" : 0.,
        # right periodic condition
        "periodicityRight_position_x" : args.channel_length,
        "periodicityRight_position_y" : 0. + args.dp,
        "periodicityRight_orientation_x" : -1.,
        "periodicityRight_orientation_y" : 0.,
        "periodicityRight_layers" : args.n_boundary_layers + 1,
        "periodicityRight_height" : args.channel_height - args.dp,
         #"periodicityRight_width" : args.n_boundary_layers * args.dp,
        "periodicityRight_width" :  1.2 * 2 * args.h_coef * args.dp,
        "periodicityRight_velocity_x" : args.v_init,
        "periodicityRight_velocity_y" : 0.,
        "periodicityRight_reference_point_z" : args.channel_length - ( -1. ) * args.n_boundary_layers * args.dp, #TODO: Make this more clear
        "periodicityRight_reference_point_x" : args.dp,
        "periodicityRight_shift_vector_x" : ( -1 ) * ( args.channel_length - 2 * ( args.dp / 2 ) ),
        "periodicityRight_shift_vector_y" : 0.
    }

    # create necessary folders
    resultsPath = r'./results'
    if not os.path.exists( resultsPath ):
        os.makedirs( resultsPath )

    sourcesPath = r'./sources'
    if not os.path.exists( sourcesPath ):
        os.makedirs( sourcesPath )

    # generate particles
    generate_channel_fluid_particles( openchannel_setup )
    generate_channel_boundary_particles( openchannel_setup )
    generate_channel_open_boundary_particles( openchannel_setup, "periodicityLeft" )
    generate_channel_open_boundary_particles( openchannel_setup, "periodicityRight" )

    # setup parameters
    compute_domain_size( openchannel_setup )

    print( openchannel_setup )
    # write simulation params
    write_simulation_params( openchannel_setup )

    #write linked list background grid
    write_domain_background_grid( openchannel_setup )
