#! /usr/bin/env python3

import numpy as np
import sys
sys.path.append('../../../src/tools')
import saveParticlesVTK

def generate_dam_break_fluid_particles( setup ):
    dp = setup[ "dp" ]
    fluid_rx = []; fluid_ry = []
    fluid_lenght_n = round( setup[ "fluid_lenght" ] / dp )
    fluid_height_n = round( setup[ "fluid_height" ] / dp )

    for x in range( fluid_lenght_n ):
        for y in range( fluid_height_n ):
            fluid_rx.append( dp * ( x + 1 ) )
            fluid_ry.append( dp * ( y + 1 ) )

    fluid_n = len( fluid_rx )
    fluid_r = np.array( ( fluid_rx, fluid_ry, np.zeros( fluid_n ) ), dtype=float ).T #!!
    fluid_v = np.zeros( ( fluid_n, 3 ) )
    fluid_rho = setup[ "density" ] * np.ones( fluid_n )
    fluid_p = np.zeros( fluid_n )
    fluid_ptype = np.zeros( fluid_n )
    fluid_to_write = saveParticlesVTK.create_pointcloud_polydata( fluid_r, fluid_v, fluid_rho, fluid_p, fluid_ptype )
    saveParticlesVTK.save_polydata( fluid_to_write, "sources/dambreak_fluid.vtk" )
    setup[ "fluid_n" ] = fluid_n


def generate_dam_break_boundary_particles( setup ):
    dp = setup[ "dp" ]
    n_boundary_layers = setup[ "n_boundary_layers" ]

    box_rx = []; box_ry = []
    ghost_rx = []; ghost_ry = []
    box_length_n = round( box_lenght / dp )
    box_height_n = round( box_height / dp )

    # left wall
    for layer in range( n_boundary_layers ):
        for y in range( box_height_n - 1 ):
            box_rx.append( 0. - layer * dp )
            box_ry.append( ( y + 1 ) * dp )
            ghost_rx.append( 0. + dp * ( layer + 1 ) )
            ghost_ry.append( ( y + 1 ) * dp )

    # bottom wall
    for layer in range( n_boundary_layers ):
        for x in range( box_length_n - n_boundary_layers + 1 ):
            box_rx.append( ( x + 1 ) * dp )
            box_ry.append( 0. - layer * dp )
            ghost_rx.append( ( x + 1 ) * dp )
            ghost_ry.append( 0. + dp * ( layer + 1 ) )

    x_last = box_rx[ -1 ] + dp #due to discretisation, we need to save last value of bottom wall

    # right wall
    for layer in range( n_boundary_layers ):
        for y in range( box_height_n - 1 ):
            box_rx.append( x_last + dp * layer )
            box_ry.append( ( y + 1 ) * dp )
            ghost_rx.append( x_last - dp * ( layer + 1 ) )
            ghost_ry.append( ( y + 1 ) * dp )

    # generate the corners
    def generate90degCorner( x, y, dirx, diry ):
      for layer in range( n_boundary_layers ):
        for k in range( n_boundary_layers ):
          box_rx.append( x + k * dp * dirx )
          box_ry.append( y + layer * dp * diry )
          ghost_rx.append( x + ( k + 1 ) * dp * dirx * ( -1 ) )
          ghost_ry.append( y + ( layer + 1 ) * dp * diry * ( -1 ) )

    generate90degCorner( 0, 0., -1, -1 )
    generate90degCorner( x_last, 0., +1, -1 )

    boundary_n = len( box_rx )
    boundary_r = np.array( ( box_rx, box_ry, np.zeros( boundary_n ) ), dtype=float ).T #!!
    boundary_ghostNodes = np.array( ( ghost_rx, ghost_ry, np.zeros( boundary_n ) ), dtype=float ).T #!!
    boundary_v = np.zeros( ( boundary_n, 3 ) )
    boundary_rho = setup[ "density" ] * np.ones( boundary_n )
    boundary_p = np.zeros( boundary_n )
    boundary_ptype = np.ones( boundary_n )
    box_to_write = saveParticlesVTK.create_pointcloud_polydata( boundary_r, boundary_v, boundary_rho, boundary_p, boundary_ptype,
                                                              ghostNodes=boundary_ghostNodes )
    saveParticlesVTK.save_polydata( box_to_write, "sources/dambreak_boundary.vtk" )

    setup[ "boundary_n" ] = boundary_n
    setup[ "domain_origin_x" ] = min( box_rx )
    setup[ "domain_origin_y" ]  = min( box_ry )
    setup[ "domain_end_x" ] = max( box_rx )
    setup[ "domain_end_y" ] = max( box_ry )

def generate_subdomains_data( setup ):
    number_of_subdomains = setup[ "number_of_subdomains" ]
    subdomains_x = setup[ "subdomains_x" ]
    subdomains_y = setup[ "subdomains_y" ]
    search_radius = setup[ "search_radius" ]

    ptcs_per_subdomain = ( setup[ "fluid_n" ] + setup[ "bound_n" ] / ( setup[ "number_of_subdomains" ] )
    grid_splits_x = []
    grid_splits_y = []

    for subdomain_x in range( subdomains_x - 1 ):
        grid_splits_x.append( math.ceil( fluid_rx[ ptcs_per_subdomain * ( subdomain_x + 1 ) ] / search_radius ) )

    for subdomain_y in range( subdomains_y - 1 ):
        grid_splits_y.append( math.ceil( fluid_ry[ ptcs_per_subdomain * ( subdomain_y + 1 ) ] / search_radius ) )

    grid_sizes_x = []
    grid_origins_x = []
    grid_index_origins_x = []

    for subdomain_x in range( subdomains_x - 1 ):
        if subdomain_x == 0:
            grid_sizes_x.append( grid_splits_x[ subdomain_x ] - 0 )
            grid_origins_x.append( setup[ "domain_origin_x" ] )
            grid_index_origins_x.append( 0 )
        elif subdomain_x == subdomains_x - 1:
            grid_sizes_x.append( setup[ "domain_sizes_x" ] - grid_splits_x[ subdomain_x - 1 ] )
            grid_origins_x.append( setup[ "domain_origin_x" ] + grid_splits_x[ subdomain_x - 1 ] * search_radius )
            grid_index_origins_x.append( grid_splits_x[ subdomain_x - 1 ] )
        else:
            grid_sizes_x.append( grid_splits_x[ subdomain_x - 1 ] - grid_index_origins_x[ subdomain_x -1 ] )
            grid_origins_x.append( setup[ "domain_origin_x" ] + grid_splits_x[ subdomain_x - 1 ] * search_radius )
            grid_index_origins_x.append( grid_sizes_x[ subdomain_x - 1 ] )

    grid_sizes_y = []
    grid_origins_y = []
    grid_index_origins_y = []

    for subdomain_y in range( subdomains_y - 1 ):
        if subdomain_y == 0:
            grid_sizes_y.append( grid_splits_y[ subdomain_y ] - 0 )
            grid_origins_y.append( setup[ "domain_origin_y" ] )
            grid_index_origins_y.append( 0 )
        elif subdomain_y == subdomains_y - 1:
            grid_sizes_y.append( setup[ "domain_sizes_y" ] - grid_splits_y[ subdomain_y - 1 ] )
            grid_origins_y.append( setup[ "domain_origin_y" ] + grid_splits_y[ subdomain_y - 1 ] * search_radius )
            grid_index_origins_y.append( grid_splits_y[ subdomain_y - 1 ] )
        else:
            grid_sizes_y.append( grid_splits_y[ subdomain_y - 1 ] - grid_index_origins_y[ subdomain_y -1 ] )
            grid_origins_y.append( setup[ "domain_origin_y" ] + grid_splits_y[ subdomain_y - 1 ] * search_radius )
            grid_index_origins_y.append( grid_sizes_y[ subdomain_y - 1 ] )

def split_to_subdomains( setup, subdomain_x, subdomain_y ):
    sub_fluid_rx = []; sub_fluid_ry = []
    sub_box_rx = []; sub_box_ry = []

    for subdomain_x in range( subdomains_x - 1 ):
        for subdomain_y in range( subdomains_y - 1 ):

    extra_parameters = {
        "domain_origin_x" : domain_origin_x,
        "domain_origin_y" : domain_origin_y,
        "domain_size_x" : domain_size_x,
        "domain_size_y" : domain_size_y
    }
    setup.update( extra_parameters )

    # load the limits of current subdomains
    if subdomain_x == 0:
        lower_limit_x = grid_origins_x[ subdomain ]
        upper_limit_x = grid_origins_x[ subdomain + 1 ]
    elif subdomain == numberOfSubdomains - 1:
        lower_limit_x = grid_origins_x[ subdomain ] - search_radius
        upper_limit_x = setup[ "domain_origin_x" ] + ( setup[ "domain_size_x" ] + 1 ) * search_radius
    else:
        lower_limit_x = grid_origins_x[ subdomain ]
        upper_limit_x = grid_origins_x[ subdomain + 1 ]

    if subdomain_y == 0:
        lower_limit_y = grid_origins_y[ subdomain ]
        upper_limit_y = grid_origins_y[ subdomain + 1 ]
    elif subdomain == numberOfSubdomains - 1:
        lower_limit_y = grid_origins_y[ subdomain ] - search_radius
        upper_limit_y = setup[ "domain_origin_y" ] + ( setup[ "domain_size_y" ] + 1 ) * search_radius
    else:
        lower_limit_y = grid_origins_y[ subdomain ]
        upper_limit_y = grid_origins_y[ subdomain + 1 ]

    # save the part for current subdomain - fluid
    for i in range ( len( fluid_rx ) ):
        if( fluid_rx[ i ] > lower_limit_x ) and ( fluid_rx[ i ] <= upper_limit_x ):
            sub_fluid_rx.append( fluid_rx[ i ] )
            sub_fluid_ry.append( fluid_ry[ i ] )

        if( fluid_rx[ i ] > lower_limit_x and fluid_rx[ i ] <= upper_limit_x + search_radius ):
            sub_fluid_rx.append( fluid_rx[ i ] )
            sub_fluid_ry.append( fluid_ry[ i ] )

    # save the part for current subdomain - boundary
    for i in range ( len( box_rx ) ):
        if( box_rx[ i ] > lower_limit_x ) and ( box_rx[ i ] <= upper_limit_x ):
            sub_box_rx.append( box_rx[ i ] )
            sub_box_ry.append( box_ry[ i ] )

        if( box_rx[ i ] > upper_limit_x and box_rx[ i ] <= upper_limit_x + searchRadius_h ):
            sub_box_rx.append( box_rx[ i ] )
            sub_box_ry.append( box_ry[ i ] )


def compute_and_write_simulation_params( dp,
                                         smoothing_length_coef,
                                         speed_of_sound,
                                         density,
                                         cfl,
                                         fluid_n,
                                         boundary_n,
                                         domain_limits_x,
                                         domain_limits_y ):

    # Compute remaining parameters
    particle_mass = density * ( dp * dp )
    smoothing_lenght =  round( smoothing_length_coef * dp, 7 )
    search_radius = round( smoothing_lenght * 2 , 7 )
    timeStep = round( cfl * ( smoothing_lenght / speed_of_sound ), 8 )

    # Resize domain by one layer of cells
    eps = 1.005
    eps_sloshing = 1.2
    domain_origin_x = eps * ( domain_limits_x[ 0 ] - search_radius )
    domain_origin_y = eps * ( domain_limits_y[ 0 ] - search_radius )
    domain_end_x = eps * ( domain_limits_x[ 1 ] + search_radius )
    domain_end_y = eps_sloshing * ( domain_limits_y[ 1 ] + search_radius ) #increase size in y due to sloshing

    domain_size_x = domain_end_x - domain_origin_x
    domain_size_y = domain_end_y - domain_origin_y

    # write parameters to config file
    with open( 'template/config_template.ini', 'r' ) as file :
      config_file = file.read()

    config_file = config_file.replace( 'placeholderSearchRadius', str( search_radius ) )
    config_file = config_file.replace( 'placeholderDomainOrigin-x', str( domain_origin_x ) )
    config_file = config_file.replace( 'placeholderDomainOrigin-y', str( domain_origin_y ) )
    config_file = config_file.replace( 'placeholderDomainSize-x', str( round( domain_size_x, 5  ) ) )
    config_file = config_file.replace( 'placeholderDomainSize-y', str( round( domain_size_y, 5  ) ) )

    config_file = config_file.replace( 'placeholderInitParticleDistance', str( dp ) )
    config_file = config_file.replace( 'placeholderSmoothingLength', str( smoothing_lenght ) )
    config_file = config_file.replace( 'placeholderMass', str( particle_mass ) )
    config_file = config_file.replace( 'placeholderSpeedOfSound', str( speed_of_sound ) )
    config_file = config_file.replace( 'placeholderDensity', str( density ))
    config_file = config_file.replace( 'placeholderTimeStep', str( timeStep ) )
    config_file = config_file.replace( 'placeholderFluidParticles', str( fluid_n ) )
    config_file = config_file.replace( 'placeholderAllocatedFluidParticles', str( fluid_n ) )
    config_file = config_file.replace( 'placeholderBoundaryParticles', str( boundary_n ) )
    config_file = config_file.replace( 'placeholderAllocatedBoundaryParticles', str( boundary_n ) )

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

    argparser = argparse.ArgumentParser(description="Heat equation example initial condition generator")
    g = argparser.add_argument_group("distribution parameters")
    g.add_argument("--subdomains-x", type=int, default=1, help="number of subdomains in x direction")
    g.add_argument("--subdomains-y", type=int, default=1, help="number of subdomains in y direction")
    g = argparser.add_argument_group("resolution parameters")
    g.add_argument("--dp", type=float, default=0.002, help="initial distance between particles")
    g.add_argument("--h-coef", type=float, default=2**0.5, help="smoothing length coefitient")
    g = argparser.add_argument_group("domain parameters")
    g.add_argument("--box-length", type=float, default=1.61, help="length of dam break box")
    g.add_argument("--box-height", type=float, default=0.8, help="height of dam break box")
    g.add_argument("--fluid-length", type=float, default=0.6, help="length of fluid block")
    g.add_argument("--fluid-height", type=float, default=0.3, help="height of fluid block")
    g.add_argument("--n-boundary-layers", type=int, default=3, help="number of boundary layers")
    g = argparser.add_argument_group("simulation parameters")
    g.add_argument("--density", type=float, default=1000, help="referential density of the fluid")
    g.add_argument("--speed-of-sound", type=float, default=34.3, help="speed of sound")
    g.add_argument("--cfl", type=float, default=0.25, help="referential density of the fluid")
    #g = argparser.add_argument_group("control parameters")
    #g.add_argument("--example-dir", type=Path, default=1000, help="referential density of the fluid")

    args = argparser.parse_args()

    dambreak_setup = {
        # general parameteres
        "dp" : args.dp,
        "h_coef" : args.h_coef,
        "n_boundary_layers" : args.n_boundary_layers,
        "density" : args.density,
        "speed_of_sound" : args.speed_of_sound,
        "cfl" : args.cfl,
        "particle_mass" : args.density * ( args.dp * args.dp ),
        "smoothing_lenght" : args.h_coef * args.dp,
        "search_radius" :  2 * args.h_coef * args.dp,
        "time_step" : args.cfl * ( args.h_coef * args.dp ) / args.speed_of_sound,
        # geometric size
        "box_length" : args.box_length,
        "box_height" : args.box_height,
        "fluid_length" : args.fluid_length,
        "fluid_height" : args.fluid_height,
        # distribution parameters
        "subdomains_x" : args.subdomains_x,
        "subdomains_y" : args.subdomains_y,
        "number_of_subdomains" : args.subdomains_x * args.subdomains_y
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
    fluid_n = generate_dam_break_fluid_particles( args.dp, args.fluid_length, args.fluid_height, args.density )
    boundary_n, domain_limits_x, domain_limits_y = generate_dam_break_boundary_particles( args.dp,
                                                                                          args.box_length,
                                                                                          args.box_height,
                                                                                          args.n_boundary_layers,
                                                                                          args.density )
    # split to subdomains

    # setup parameters
    compute_and_write_simulation_params( args.dp,
                                         args.h_coef,
                                         args.speed_of_sound,
                                         args.density,
                                         args.cfl,
                                         fluid_n,
                                         boundary_n,
                                         domain_limits_x,
                                         domain_limits_y )

    # setup measuretool
    configure_and_write_measuretool_parameters()
