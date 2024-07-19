#! /usr/bin/env python3

import numpy as np
import sys
sys.path.append('../../../src/tools')
import saveParticlesVTK
import domainGrid

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
    ghost_rx = []; ghost_ry = []; ghost_rz = []
    box_density = []
    box_length_n = round( setup[ "channel_length" ] / dp )
    box_height_n = round( ( 1.5 * setup[ "channel_height" ] ) / dp )
    box_width_n = round( setup[ "channel_width" ] / dp )

    box_x_ref_point = -0.2
    box_y_ref_point = -0.5 * setup[ "channel_width" ]
    box_z_ref_point = 0

    rho0 = setup[ 'density' ]
    speed_of_sound = setup[ 'speed_of_sound' ]

    # bottom wall
    for layer in range( n_boundary_layers ):
        for y in range( box_width_n - 1 ):
            for x in range( box_length_n + ( n_boundary_layers ) * 2 + 1 ):
                box_rx.append( box_x_ref_point + ( x - ( n_boundary_layers ) ) * dp )
                box_ry.append( box_y_ref_point + dp + y * dp )
                box_rz.append( box_z_ref_point + 0. - layer * dp )
                ghost_rx.append( 0. )
                ghost_ry.append( 0. )
                ghost_rz.append( 0. )

                hydrostaticPressure = rho0 * 9.81 * ( setup[ 'channel_height' ] - box_rz[ -1 ] )
                hydrostaticDensity = ( ( hydrostaticPressure / ( speed_of_sound** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0
                box_density.append( hydrostaticDensity )

    # left wall
    for layer in range( n_boundary_layers ):
        for z in range( box_height_n ):
            for x in range( box_length_n + ( n_boundary_layers ) * 2 + 1 ):
                box_rx.append( box_x_ref_point + ( x - ( n_boundary_layers ) ) * dp )
                box_ry.append( box_y_ref_point + 0. - layer * dp )
                box_rz.append( box_z_ref_point + dp + z * dp )
                ghost_rx.append( 0. )
                ghost_ry.append( 0. )
                ghost_rz.append( 0. )

                hydrostaticPressure = rho0 * 9.81 * ( setup[ 'channel_height' ] - box_rz[ -1 ] )
                hydrostaticDensity = ( ( hydrostaticPressure / ( speed_of_sound** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0
                box_density.append( hydrostaticDensity )

    # right wall
    for layer in range( n_boundary_layers ):
        for z in range( box_height_n ):
            for x in range( box_length_n + ( n_boundary_layers ) * 2 + 1 ):
                box_rx.append( box_x_ref_point + ( x - ( n_boundary_layers ) ) * dp )
                box_ry.append( box_y_ref_point + ( setup[ "channel_width" ] ) + layer * dp )
                box_rz.append( box_z_ref_point + dp + z * dp )
                ghost_rx.append( 0. )
                ghost_ry.append( 0. )
                ghost_rz.append( 0. )

                hydrostaticPressure = rho0 * 9.81 * ( setup[ 'channel_height' ] - box_rz[ -1 ] )
                hydrostaticDensity = ( ( hydrostaticPressure / ( speed_of_sound** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0
                box_density.append( hydrostaticDensity )

    # generate corner lines
    def generate90degCorner( y, z, diry, dirz ):
      for layer_y in range( n_boundary_layers ):
        for layer_z in range( n_boundary_layers ):
            for x in range( box_length_n + ( n_boundary_layers ) * 2 + 1 ):
                box_rx.append( box_x_ref_point + ( x - ( n_boundary_layers ) ) * dp )
                box_ry.append( box_y_ref_point + y + layer_y * dp * diry )
                box_rz.append( box_z_ref_point + z + layer_z * dp * dirz )
                ghost_rx.append( 0. )
                ghost_ry.append( 0. )
                ghost_rz.append( 0. )

                hydrostaticPressure = rho0 * 9.81 * ( setup[ 'channel_height' ] - box_rz[ -1 ] )
                hydrostaticDensity = ( ( hydrostaticPressure / ( speed_of_sound** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0
                box_density.append( hydrostaticDensity )

    generate90degCorner( 0, 0., -1, -1 )
    generate90degCorner( setup[ "channel_width" ], 0., +1, -1 )

    boundary_n = len( box_rx )
    boundary_r = np.array( ( box_rx, box_ry, box_rz ), dtype=float ).T #!!
    boundary_ghostNodes = np.array( ( ghost_rx, ghost_ry, np.zeros( boundary_n ) ), dtype=float ).T #!!
    boundary_v = np.zeros( ( boundary_n, 3 ) )
    #boundary_rho = setup[ "density" ] * np.ones( boundary_n )
    boundary_rho = np.array( box_density, dtype=float )
    boundary_p = np.zeros( boundary_n )
    boundary_ptype = np.ones( boundary_n )
    box_to_write = saveParticlesVTK.create_pointcloud_polydata( boundary_r, boundary_v, boundary_rho, boundary_p, boundary_ptype,
                                                                ghostNodes=boundary_ghostNodes )
    saveParticlesVTK.save_polydata( box_to_write, "sources/openchannel_boundary.vtk" )

    setup[ "boundary_n" ] = boundary_n
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

    patch_width_n = setup[ prefix + "_layers" ]
    patch_size_y_n = round( setup[ prefix + "_size_y" ] / dp  )
    patch_size_z_n = round( setup[ prefix + "_size_z" ] / dp  )

    rho0 = setup[ 'density' ]
    speed_of_sound = setup[ 'speed_of_sound' ]

    for x in range( patch_width_n ):
        for y in range( patch_size_y_n ):
            for z in range( patch_size_z_n ):
                patch_rx.append( setup[ prefix + "_position_x" ] - setup[ prefix + "_orientation_x" ] * dp *  x  )
                patch_ry.append( setup[ prefix + "_position_y" ] + dp * ( y ) )
                patch_rz.append( setup[ prefix + "_position_z" ] + dp * ( z ) )
                patch_vx.append( setup[ prefix + "_velocity_x" ] )
                patch_vy.append( setup[ prefix + "_velocity_y" ] )
                patch_vz.append( setup[ prefix + "_velocity_z" ] )

                hydrostaticPressure = rho0 * 9.81 * ( setup[ 'channel_height' ] - patch_rz[ -1 ] )
                hydrostaticDensity = ( ( hydrostaticPressure / ( speed_of_sound** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0
                patch_density.append( hydrostaticDensity )

    patch_n = len( patch_rx )
    r = np.array( ( patch_rx, patch_ry, patch_rz ), dtype=float ).T #!!
    v = np.array( ( patch_vx, patch_vy, patch_vz ), dtype=float ).T #!!
    #rho = setup[ "density" ] * np.ones( patch_n, dtype=float )
    rho = np.array( patch_density, dtype=float )
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
    config_file = config_file.replace( 'placeholderSmoothingLength', str( round( setup[ "smoothing_lenght" ], 7 ) ) )
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

    argparser = argparse.ArgumentParser(description="Heat equation example initial condition generator")
    g = argparser.add_argument_group("resolution parameters")
    g.add_argument("--dp", type=float, default=0.01, help="initial distance between particles")
    g.add_argument("--h-coef", type=float, default=2, help="smoothing length coefitient")
    g = argparser.add_argument_group("domain parameters")
    g.add_argument("--channel-length", type=float, default=2.2, help="length of fluid block")
    g.add_argument("--channel-height", type=float, default=0.2, help="height of fluid block")
    g.add_argument("--channel-width", type=float, default=1.0, help="width of fluid block")
    g.add_argument("--n-boundary-layers", type=int, default=3, help="number of boundary layers")
    g.add_argument("--v-init", type=int, default=1.5, help="initial velocity")
    g = argparser.add_argument_group("simulation parameters")
    g.add_argument("--density", type=float, default=1000, help="referential density of the fluid")
    g.add_argument("--speed-of-sound", type=float, default=60.0, help="speed of sound")
    g.add_argument("--cfl", type=float, default=0.2, help="referential density of the fluid")
    #g = argparser.add_argument_group("control parameters")
    #g.add_argument("--example-dir", type=Path, default=1000, help="referential density of the fluid")

    args = argparser.parse_args()

    openchannel_setup = {
        # general parameteres
        "dp" : args.dp,
        "h_coef" : args.h_coef,
        "n_boundary_layers" : args.n_boundary_layers,
        "density" : args.density,
        "speed_of_sound" : args.speed_of_sound,
        "cfl" : args.cfl,
        "particle_mass" : args.density * ( args.dp * args.dp * args.dp ),
        "smoothing_lenght" : args.h_coef * args.dp,
        "search_radius" :  2 * args.h_coef * args.dp,
        "time_step" : args.cfl * ( args.h_coef * args.dp ) / args.speed_of_sound,
        # geometric size
        "channel_length" : args.channel_length,
        "channel_height" : args.channel_height + args.dp,
        "channel_width" : args.channel_width,
        # inlet boundary condition
        "inlet_position_x" : -0.2,
        "inlet_position_y" : -0.5 * args.channel_width + args.dp,
        "inlet_position_z" : 0. + args.dp,
        "inlet_orientation_x" : 1.,
        "inlet_orientation_y" : 0.,
        "inlet_orientation_z" : 0.,
        "inlet_layers" : args.n_boundary_layers + 1,
        "inlet_size_y" : args.channel_width - args.dp,
        "inlet_size_z" : args.channel_height,
        "inlet_width" : ( args.n_boundary_layers + 1 ) * args.dp,
        "inlet_velocity_x" : args.v_init,
        "inlet_velocity_y" : 0.,
        "inlet_velocity_z" : 0.,
        # outlet boundary condition
        "outlet_position_x" : -0.2 + args.channel_length,
        "outlet_position_y" : -0.5 * args.channel_width + args.dp,
        "outlet_position_z" : 0. + args.dp,
        "outlet_orientation_x" : -1.,
        "outlet_orientation_y" : 0.,
        "outlet_orientation_z" : 0.,
        "outlet_layers" : args.n_boundary_layers + 1,
        "outlet_size_y" : args.channel_width - args.dp,
        "outlet_size_z" : args.channel_height,
        "outlet_width" : ( args.n_boundary_layers + 1 ) * args.dp,
        "outlet_velocity_x" : args.v_init,
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

    # generate particles
    generate_channel_fluid_particles( openchannel_setup )
    generate_channel_boundary_particles( openchannel_setup )
    generate_channel_open_boundary_particles( openchannel_setup, "inlet" )
    generate_channel_open_boundary_particles( openchannel_setup, "outlet" )

    # setup parameters
    compute_domain_size( openchannel_setup )

    print( openchannel_setup )
    # write simulation params
    write_simulation_params( openchannel_setup )

    #write linked list background grid
    write_domain_background_grid( openchannel_setup )
