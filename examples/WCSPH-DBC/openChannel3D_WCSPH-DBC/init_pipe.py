#! /usr/bin/env python3

import numpy as np
import sys
sys.path.append('../../../src/tools')
import saveParticlesVTK
import domainGrid

def velocity_paraboloid_profile( r_x, r_y, r_z, setup ):
    v_max = 0.5
    y_center = 0
    z_center = 0

    v_x = v_max - v_max * ( ( np.abs( r_y - y_center ) / setup[ 'inlet_radius' ] )**2 + ( np.abs( r_z - z_center ) / setup[ 'inlet_radius' ] )**2 )
    v_y = 0.
    v_z = 0.
    return v_x, v_y, v_z

def velocity_constant_profile( r_x, r_y, r_z, setup ):
    return 0.5, 0, 0

def generate_channel_fluid_particles( setup ):
    dp = setup[ "dp" ]
    fluid_rx = []; fluid_ry = []; fluid_rz = []
    fluid_vx = []; fluid_vy = []; fluid_vz = []
    fluid_density = []
    fluid_lenght_n = round( ( setup[ "channel_length" ] - dp ) / dp )
    fluid_height_n = round( ( setup[ "channel_height" ] - dp ) / dp )
    fluid_width_n = round( ( setup[ "channel_width" ] - dp ) / dp )

    fluid_x_ref_point = -0.2
    fluid_y_ref_point = -0.5 * setup[ "channel_width" ]
    fluid_z_ref_point = 0

    rho0 = setup[ 'density' ]
    speed_of_sound = setup[ 'speed_of_sound' ]

    for x in range( fluid_lenght_n ):
        for y in range( fluid_width_n ):
            for z in range( fluid_height_n ):
                fluid_rx.append( fluid_x_ref_point + dp * ( x + 1 ) )
                fluid_ry.append( fluid_y_ref_point + dp * ( y + 1 ) )
                fluid_rz.append( fluid_z_ref_point + dp * ( z + 1 ) )
                fluid_vx.append( setup[ "inlet_velocity_x" ] )
                fluid_vy.append( setup[ "inlet_velocity_y" ] )
                fluid_vz.append( setup[ "inlet_velocity_z" ] )

                hydrostaticPressure = rho0 * 9.81 * ( setup[ 'channel_height' ] - fluid_rz[ -1 ] )
                hydrostaticDensity = ( ( hydrostaticPressure / ( speed_of_sound** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0
                fluid_density.append( hydrostaticDensity )

    fluid_n = len( fluid_rx )
    fluid_r = np.array( ( fluid_rx, fluid_ry, fluid_rz ), dtype=float ).T #!!
    fluid_v = np.array( ( fluid_vx, fluid_vy, fluid_vz ), dtype=float ).T #!!
    #fluid_rho = setup[ "density" ] * np.ones( fluid_n )
    fluid_rho = np.array( fluid_density, dtype=float )
    fluid_p = np.zeros( fluid_n )
    fluid_ptype = np.zeros( fluid_n )
    fluid_to_write = saveParticlesVTK.create_pointcloud_polydata( fluid_r, fluid_v, fluid_rho, fluid_p, fluid_ptype )
    saveParticlesVTK.save_polydata( fluid_to_write, "sources/openchannel_fluid.vtk" )
    setup[ "fluid_n" ] = fluid_n

def generate_channel_boundary_particles( setup ):
    dp = setup[ "dp" ]
    n_boundary_layers = setup[ "n_boundary_layers" ]
    box_rx = []; box_ry = []; box_rz = []
    #ghost_rx = []; ghost_ry = []; ghost_rz = []

    inlet_pipe_radius = setup[ "pipe_radius" ]
    inlet_pipe_center_y = setup[ "inlet_center_y" ]
    inlet_pipe_center_z = setup[ "inlet_center_z" ]
    inlet_pipe_length = setup[ "pipe_length" ]

    inlet_pipe_pos_x = 0.8
    # position = position - offset to include additional layers
    inlet_pipe_pos_y = setup[ "inlet_position_y" ]
    inlet_pipe_pos_z = setup[ "inlet_position_z" ]

    inlet_pipe_length_n = round( inlet_pipe_length / dp )
    # size_n = number of points to cover the radius + the last point + offset to include additional layers
    inlet_pipe_size_y_n = round( ( 2 * ( inlet_pipe_radius ) ) / dp ) + 1 + ( n_boundary_layers + 2 ) * 2
    inlet_pipe_size_z_n = round( ( 2 * ( inlet_pipe_radius ) ) / dp ) + 1 + ( n_boundary_layers + 2 ) * 2

    for layer in range( n_boundary_layers ):
            for y in range( inlet_pipe_size_y_n ):
                for z in range( inlet_pipe_size_z_n ):
                    pos_y = inlet_pipe_pos_y + dp * y
                    pos_z = inlet_pipe_pos_z + dp * z

                    radius_xy = np.sqrt( ( pos_y - inlet_pipe_center_y )**2 + ( pos_z - inlet_pipe_center_z )**2 )
                    lowerBound = inlet_pipe_radius + layer * dp - dp / 2
                    upperBound = inlet_pipe_radius + layer * dp + dp / 2
                    if layer == 0:
                        lowerBound = inlet_pipe_radius + layer * dp

                    if radius_xy < upperBound and radius_xy >= lowerBound:
                        for x in range( inlet_pipe_length_n ):
                            box_rx.append( inlet_pipe_pos_x - dp * ( x + 1 ) )
                            box_ry.append( pos_y )
                            box_rz.append( pos_z )

    box_n = len( box_rx )
    box_r = np.array( ( box_rx, box_ry, box_rz ), dtype=float ).T #!!
    #box_ghostNodes = np.array( ( ghost_rx, ghost_ry, np.zeros( boundary_n ) ), dtype=float ).T #!!
    box_v = np.zeros( ( box_n, 3 ) )
    box_rho = setup[ "density" ] * np.ones( box_n )
    box_p = np.zeros( box_n )
    box_ptype = np.ones( box_n )
    box_to_write = saveParticlesVTK.create_pointcloud_polydata( box_r, box_v, box_rho, box_p, box_ptype )
    saveParticlesVTK.save_polydata( box_to_write, "sources/openchannel_boundary.vtk" )

    setup[ "boundary_n" ] = box_n
    setup[ "domain_origin_x" ] = min( box_rx )
    setup[ "domain_origin_y" ]  = min( box_ry )
    setup[ "domain_origin_z" ]  = min( box_rz )
    setup[ "domain_end_x" ] = max( box_rx )
    setup[ "domain_end_y" ] = max( box_ry )
    setup[ "domain_end_z" ] = max( box_rz )

def generate_channel_open_boundary_particles( setup, prefix ):
    dp = setup[ "dp" ]
    patch_rx = []; patch_ry = []; patch_rz = []
    patch_vx = []; patch_vy = []; patch_vz = []
    patch_density = []

    patch_pos_x = setup[ prefix + "_position_x" ]
    patch_pos_y = setup[ prefix + "_position_y" ]
    patch_pos_z = setup[ prefix + "_position_z" ]

    patch_radius = setup[ prefix + "_radius" ]
    patch_center_x = setup[ prefix + "_center_x" ]
    patch_center_y = setup[ prefix + "_center_y" ]
    patch_center_z = setup[ prefix + "_center_z" ]

    patch_width_n = setup[ prefix + "_layers" ]
    patch_size_y_n = round( 2 * patch_radius / dp ) + 2
    patch_size_z_n = round( 2 * patch_radius / dp ) + 2

    patch_orientation_x = setup[ prefix + "_orientation_x" ]

    for x in range( patch_width_n ):
        for y in range( patch_size_y_n ):
            for z in range( patch_size_z_n ):
                pos_x = patch_pos_x + ( -1 ) * patch_orientation_x * dp * x
                pos_y = patch_pos_y + dp * y
                pos_z = patch_pos_z + dp * z

                if np.sqrt( ( pos_y - patch_center_y )**2 + ( pos_z - patch_center_z )**2 ) < patch_radius:
                    patch_rx.append( pos_x )
                    patch_ry.append( pos_y )
                    patch_rz.append( pos_z )

                    v_x, v_y, v_z = velocity_profile( patch_rx[ -1 ], patch_ry[ -1 ], patch_rz[ -1 ], setup )
                    patch_vx.append( v_x )
                    patch_vy.append( v_y )
                    patch_vz.append( v_z )

    patch_n = len( patch_rx )
    r = np.array( ( patch_rx, patch_ry, patch_rz ), dtype=float ).T #!!
    v = np.array( ( patch_vx, patch_vy, patch_vz ), dtype=float ).T #!!
    rho = setup[ "density" ] * np.ones( patch_n, dtype=float )
    #rho = np.array( patch_density, dtype=float )
    p = np.zeros( patch_n )
    ptype = np.ones( patch_n )

    inletToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
    saveParticlesVTK.save_polydata( inletToWrite, "sources/openchannel_" + prefix + ".vtk" )

    setup[ prefix + "_n" ] = patch_n

def compute_domain_size( setup ):
    search_radius = setup[ "search_radius" ]

    # Resize domain by one layer of cells
    eps = 1.005
    domain_origin_x = eps * ( setup[ "domain_origin_x" ] - search_radius )
    domain_origin_y = eps * ( setup[ "domain_origin_y" ] - search_radius )
    domain_origin_z = eps * ( setup[ "domain_origin_z" ] - search_radius )
    domain_end_x = eps * ( setup[ "domain_end_x" ] + search_radius )
    domain_end_y = eps * ( setup[ "domain_end_y" ] + search_radius )
    domain_end_z = eps * ( setup[ "domain_end_z" ] + search_radius )
    domain_size_x = domain_end_x - domain_origin_x
    domain_size_y = domain_end_y - domain_origin_y
    domain_size_z = domain_end_z - domain_origin_z

    extra_parameters = {
        "domain_origin_x" : domain_origin_x,
        "domain_origin_y" : domain_origin_y,
        "domain_origin_z" : domain_origin_z,
        "domain_size_x" : domain_size_x,
        "domain_size_y" : domain_size_y,
        "domain_size_z" : domain_size_z,
    }
    setup.update( extra_parameters )

def write_simulation_params( setup ):

    # write parameters to config file
    with open( 'template/config_template.ini', 'r' ) as file :
      config_file = file.read()

    config_file = config_file.replace( 'placeholderSearchRadius', str( round( setup[ "search_radius" ], 7 ) ) )
    config_file = config_file.replace( 'placeholderDomainOrigin-x', str( round( setup[ "domain_origin_x" ], 5 ) ) )
    config_file = config_file.replace( 'placeholderDomainOrigin-y', str( round( setup[ "domain_origin_y" ], 5 ) ) )
    config_file = config_file.replace( 'placeholderDomainOrigin-z', str( round( setup[ "domain_origin_z" ], 5 ) ) )
    config_file = config_file.replace( 'placeholderDomainSize-x', str( round( setup[ "domain_size_x" ], 5  ) ) )
    config_file = config_file.replace( 'placeholderDomainSize-y', str( round( setup[ "domain_size_y" ], 5  ) ) )
    config_file = config_file.replace( 'placeholderDomainSize-z', str( round( setup[ "domain_size_z" ], 5  ) ) )

    config_file = config_file.replace( 'placeholderInitParticleDistance', str( setup[ "dp" ] ) )
    config_file = config_file.replace( 'placeholderSmoothingLength', str( round( setup[ "smoothing_length" ], 7 ) ) )
    config_file = config_file.replace( 'placeholderMass', str( round( setup[ "particle_mass" ], 7 ) ) )
    config_file = config_file.replace( 'placeholderSpeedOfSound', str( setup[ "speed_of_sound" ] ) )
    config_file = config_file.replace( 'placeholderDensity', str( setup[ "density" ] ) )
    config_file = config_file.replace( 'placeholderTimeStep', str( setup[ "time_step" ] ) )
    config_file = config_file.replace( 'placeholderFluidParticles', str( setup[ "fluid_n" ] ) )
    config_file = config_file.replace( 'placeholderAllocatedFluidParticles', str( 2 * setup[ "fluid_n" ] ) )
    config_file = config_file.replace( 'placeholderBoundaryParticles', str( setup[ "boundary_n" ] ) )
    config_file = config_file.replace( 'placeholderAllocatedBoundaryParticles', str( setup[ "boundary_n" ] ) )
    config_file = config_file.replace( 'placeholderInletParticles', str( setup[ "inlet_n" ] ) )
    config_file = config_file.replace( 'placeholderAllocatedInletParticles', str( setup[ "inlet_n" ] ) ) #TODO 3 *
    config_file = config_file.replace( 'placeholderOutletParticles', str( setup[ "outlet_n" ] ) )
    config_file = config_file.replace( 'placeholderAllocatedOutletParticles', str( 3 * setup[ "outlet_n" ] ) )

    config_file = config_file.replace( 'placeholderInletOrientation_x', str( setup[ "inlet_orientation_x" ] ) )
    config_file = config_file.replace( 'placeholderInletOrientation_y', str( setup[ "inlet_orientation_y" ] ) )
    config_file = config_file.replace( 'placeholderInletOrientation_z', str( setup[ "inlet_orientation_z" ] ) )
    config_file = config_file.replace( 'placeholderInletVelocity_x', str( setup[ "inlet_velocity_x" ] ) )
    config_file = config_file.replace( 'placeholderInletVelocity_y', str( setup[ "inlet_velocity_y" ] ) )
    config_file = config_file.replace( 'placeholderInletVelocity_z', str( setup[ "inlet_velocity_z" ] ) )
    config_file = config_file.replace( 'placeholderInletPosition1_x', str( setup[ "inlet_position_x" ]  + setup[ "dp" ] / 2 ) ) #TODO
    config_file = config_file.replace( 'placeholderInletPosition1_y', str( setup[ "inlet_position_y" ] ) )
    config_file = config_file.replace( 'placeholderInletPosition1_z', str( setup[ "inlet_position_z" ] ) )
    config_file = config_file.replace( 'placeholderInletPosition2_x', str( setup[ "inlet_position_x" ]  + setup[ "dp" ] / 2 ) ) #TODO
    config_file = config_file.replace( 'placeholderInletPosition2_y', str( setup[ "inlet_position_y" ] + setup[ "inlet_size_y" ] ) )# - setup[ "dp" ] / 2 ) )
    config_file = config_file.replace( 'placeholderInletPosition2_z', str( setup[ "inlet_position_z" ] + setup[ "inlet_size_z" ] ) )# - setup[ "dp" ] / 2 ) )
    config_file = config_file.replace( 'placeholderInletDensity', str( setup[ "density" ] ) )
    config_file = config_file.replace( 'placeholderInletWidth_x', str( round( setup[ "inlet_width" ], 7 ) ) )
    config_file = config_file.replace( 'placeholderInletWidth_y', str( 0. ) )
    config_file = config_file.replace( 'placeholderInletWidth_z', str( 0. ) )

    config_file = config_file.replace( 'placeholderOutletOrientation_x', str( setup[ "outlet_orientation_x" ] ) )
    config_file = config_file.replace( 'placeholderOutletOrientation_y', str( setup[ "outlet_orientation_y" ] ) )
    config_file = config_file.replace( 'placeholderOutletOrientation_z', str( setup[ "outlet_orientation_z" ] ) )
    config_file = config_file.replace( 'placeholderOutletVelocity_x', str( setup[ "outlet_velocity_x" ] ) )
    config_file = config_file.replace( 'placeholderOutletVelocity_y', str( setup[ "outlet_velocity_y" ] ) )
    config_file = config_file.replace( 'placeholderOutletVelocity_z', str( setup[ "outlet_velocity_z" ] ) )
    config_file = config_file.replace( 'placeholderOutletPosition1_x', str( setup[ "outlet_position_x" ] - setup[ "dp" ] / 2 ) ) #TODO
    config_file = config_file.replace( 'placeholderOutletPosition1_y', str( setup[ "outlet_position_y" ] ) )
    config_file = config_file.replace( 'placeholderOutletPosition1_z', str( setup[ "outlet_position_z" ] ) )
    config_file = config_file.replace( 'placeholderOutletPosition2_x', str( setup[ "outlet_position_x" ] - setup[ "dp" ] / 2 ) ) #FIXME
    config_file = config_file.replace( 'placeholderOutletPosition2_y', str( setup[ "outlet_position_y" ] + setup[ "outlet_size_y" ] ) )# - setup[ "dp" ] / 2  ) )
    config_file = config_file.replace( 'placeholderOutletPosition2_z', str( setup[ "outlet_position_z" ] + setup[ "outlet_size_z" ] ) )# - setup[ "dp" ] / 2  ) )
    config_file = config_file.replace( 'placeholderOutletDensity', str( setup[ "density" ] ) )
    config_file = config_file.replace( 'placeholderOutletWidth_x', str( round( setup[ "outlet_width" ], 7 ) ) )
    config_file = config_file.replace( 'placeholderOutletWidth_y', str( 0. ) )
    config_file = config_file.replace( 'placeholderOutletWidth_z', str( 0. ) )

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
    from pprint import pprint

    argparser = argparse.ArgumentParser(description="Heat equation example initial condition generator")
    g = argparser.add_argument_group("resolution parameters")
    g.add_argument("--dp", type=float, default=0.005, help="initial distance between particles")
    g.add_argument("--h-coef", type=float, default=3**0.5, help="smoothing length coefitient")
    g = argparser.add_argument_group("domain parameters")
    g.add_argument("--pipe-length", type=float, default=1.5, help="length of pipe")
    g.add_argument("--pipe-radius", type=float, default=0.06, help="radius of pipe")
    g.add_argument("--n-boundary-layers", type=int, default=3, help="number of boundary layers")
    g.add_argument("--v-init", type=int, default=1.5, help="initial velocity")
    g = argparser.add_argument_group("simulation parameters")
    g.add_argument("--density", type=float, default=1000, help="referential density of the fluid")
    g.add_argument("--speed-of-sound", type=float, default=56.0286, help="speed of sound")
    g.add_argument("--cfl", type=float, default=0.2, help="referential density of the fluid")
    #g = argparser.add_argument_group("control parameters")
    #g.add_argument("--example-dir", type=Path, default=1000, help="referential density of the fluid")

    args = argparser.parse_args()

    channel_with_pipe_setup = {
        # general parameteres
        "dp" : args.dp,
        "h_coef" : args.h_coef,
        "n_boundary_layers" : args.n_boundary_layers,
        "density" : args.density,
        "speed_of_sound" : args.speed_of_sound,
        "cfl" : args.cfl,
        "particle_mass" : args.density * ( args.dp * args.dp * args.dp ),
        "smoothing_length" : args.h_coef * args.dp,
        "search_radius" :  2 * args.h_coef * args.dp,
        "time_step" : args.cfl * ( args.h_coef * args.dp ) / args.speed_of_sound,
        # geometric size
        "pipe_length" : args.pipe_length,
        "pipe_radius" : args.pipe_radius,
        # inlet boundary condition
        "inlet_position_x" : -0.6,
        "inlet_position_y" : -0.07,
        "inlet_position_z" : -0.07,
        "inlet_orientation_x" : 1.,
        "inlet_orientation_y" : 0.,
        "inlet_orientation_z" : 0.,
        "inlet_layers" : args.n_boundary_layers + 1,
        "inlet_radius" : args.pipe_radius,
        "inlet_center_x" : -0.6,
        "inlet_center_y" : 0.,
        "inlet_center_z" : 0.,
        "inlet_size_y" : 0.14,
        "inlet_size_z" : 0.14,
        "inlet_width" : ( args.n_boundary_layers + 1 ) * args.dp,
        "inlet_velocity_x" : args.v_init,
        "inlet_velocity_y" : 0.,
        "inlet_velocity_z" : 0.,
        # outlet boundary conditions
        "outlet_position_x" : 0.7,
        "outlet_position_y" : -0.07,
        "outlet_position_z" : -0.07,
        "outlet_orientation_x" : -1.,
        "outlet_orientation_y" : 0.,
        "outlet_orientation_z" : 0.,
        "outlet_layers" : args.n_boundary_layers + 1,
        "outlet_radius" : args.pipe_radius,
        "outlet_center_x" : 0.3,
        "outlet_center_y" : 0.,
        "outlet_center_z" : 0.,
        "outlet_size_y" : 0.14,
        "outlet_size_z" : 0.14,
        "outlet_width" : ( args.n_boundary_layers + 1 ) * args.dp,
        "outlet_velocity_x" : 0.,
        "outlet_velocity_y" : 0.,
        "outlet_velocity_z" : 0.,
    }

    # create necessary folders
    resultsPath = r'./results'
    if not os.path.exists( resultsPath ):
        os.makedirs( resultsPath )

    sourcesPath = r'./sources'
    if not os.path.exists( sourcesPath ):
        os.makedirs( sourcesPath )

    # define and plot inlet velocity profile
    velocity_profile = velocity_constant_profile
    #plot_inlet_velocity_profile( velocity_profile, channel_with_pipe_setup )

    # generate particles
    #generate_channel_fluid_particles( channel_with_pipe_setup )
    channel_with_pipe_setup[ "fluid_n" ] = 0
    generate_channel_boundary_particles( channel_with_pipe_setup )
    generate_channel_open_boundary_particles( channel_with_pipe_setup, "inlet" )
    generate_channel_open_boundary_particles( channel_with_pipe_setup, "outlet" )

    # setup parameters
    compute_domain_size( channel_with_pipe_setup )

    pprint( channel_with_pipe_setup )
    # write simulation params
    write_simulation_params( channel_with_pipe_setup )

    #write linked list background grid
    write_domain_background_grid( channel_with_pipe_setup )
