#! /usr/bin/env python3

import math
import numpy as np
import sys
sys.path.append('../../../src/tools')
import saveParticlesVTK

def generate_dam_break_fluid_particles( setup ):
    dp = setup[ "dp" ]
    fluid_rx = []; fluid_ry = []
    fluid_length_n = round( setup[ "fluid_length" ] / dp )
    fluid_height_n = round( setup[ "fluid_height" ] / dp )

    for x in range( fluid_length_n ):
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

    return fluid_rx, fluid_ry

def generate_dam_break_boundary_particles( setup ):
    dp = setup[ "dp" ]
    n_boundary_layers = setup[ "n_boundary_layers" ]

    box_rx = []; box_ry = []
    ghost_rx = []; ghost_ry = []
    box_length_n = round( setup[ "box_length" ] / dp )
    box_height_n = round( setup[ "box_height" ] / dp )

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

    return box_rx, box_ry

def generate_subdomains_data( setup, fluid_rx, fluid_ry, box_rx, box_ry ):
    number_of_subdomains = setup[ "number_of_subdomains" ]
    subdomains_x = setup[ "subdomains_x" ]
    subdomains_y = setup[ "subdomains_y" ]
    search_radius = setup[ "search_radius" ]

    ptcs_per_subdomain = ( int )( math.ceil( setup[ "fluid_n" ] + setup[ "boundary_n" ] ) / ( setup[ "number_of_subdomains" ] ) )
    grid_splits_x = []
    grid_splits_y = []

    for subdomain_x in range( subdomains_x - 1 ):
        grid_splits_x.append( math.ceil( fluid_rx[ ptcs_per_subdomain * ( subdomain_x + 1 ) ] / search_radius ) )
    for subdomain_y in range( subdomains_y - 1 ):
        grid_splits_y.append( math.ceil( fluid_ry[ ptcs_per_subdomain * ( subdomain_y + 1 ) ] / search_radius ) )

    domain_sizes_x = []
    grid_sizes_x = []
    grid_origins_x = []
    grid_index_origins_x = []

    if subdomains_x == 1:
        #grid_sizes_y.append( setup[ "grid_size_y" ] ) # there is not variable grid size, only domain size
        grid_origins_x.append( setup[ "domain_origin_x" ] )
        grid_index_origins_x.append( 0 )
        domain_sizes_x( setup[ "domain_size_x" ] )
    else:
        #for subdomain_x in range( subdomains_x - 1 ):
        for subdomain_x in range( subdomains_x ):
            if subdomain_x == 0:
                grid_sizes_x.append( grid_splits_x[ subdomain_x ] - 0 )
                grid_origins_x.append( setup[ "domain_origin_x" ] )
                grid_index_origins_x.append( 0 )
                domain_sizes_x.append( grid_splits_x[ subdomain_x ] * search_radius )
            elif subdomain_x == subdomains_x - 1:
                #grid_sizes_x.append( setup[ "domain_size_x" ] - grid_splits_x[ subdomain_x - 1 ] ) //TODO:
                grid_sizes_x.append( math.ceil( setup[ "domain_size_x" ] / search_radius ) - grid_splits_x[ subdomain_x - 1 ] )
                #grid_origins_x.append( setup[ "domain_origin_x" ] + grid_splits_x[ subdomain_x - 1 ] * search_radius )
                #grid_index_origins_x.append( grid_splits_x[ subdomain_x - 1 ] * search_radius )
                grid_origins_x.append( setup[ "domain_origin_x" ] + search_radius * np.sum( grid_sizes_x[ 0 : subdomain_x ] ) ) #TODO: Use this, added domain origin
                grid_index_origins_x.append( np.sum( grid_sizes_x[ 0 : subdomain_x ] ) ) #TODO: Use this.
                domain_sizes_x.append( setup[ "domain_size_x" ]  - grid_splits_x[ subdomain_x - 1 ] * search_radius )
            else:
                grid_sizes_x.append( grid_splits_x[ subdomain_x - 1 ] - grid_index_origins_x[ subdomain_x -1 ] )
                grid_origins_x.append( setup[ "domain_origin_x" ] + grid_splits_x[ subdomain_x - 1 ] * search_radius )
                grid_index_origins_x.append( grid_sizes_x[ subdomain_x - 1 ] )
                #grid_index_origins_x.append( grid_splits_x[ subdomain_x - 1 ] )
                domain_sizes_x.append( grid_splits_x[ subdomain_x - 1 ] * search_radius )

    domain_sizes_y = []
    grid_sizes_y = []
    grid_origins_y = []
    grid_index_origins_y = []

    if subdomains_y == 1:
        #grid_sizes_y.append( setup[ "grid_size_y" ] ) # there is not variable grid size, only domain size
        grid_origins_y.append( setup[ "domain_origin_y" ] )
        grid_index_origins_y.append( 0 )
        domain_sizes_y.append( setup[ "domain_size_y" ] )
    else:
        #for subdomain_y in range( subdomains_y - 1 ):
        for subdomain_y in range( subdomains_y ):
            if subdomain_y == 0:
                grid_sizes_y.append( grid_splits_y[ subdomain_y ] - 0 )
                grid_origins_y.append( setup[ "domain_origin_y" ] )
                grid_index_origins_y.append( 0 )
                domain_sizes_y( grid_splits_y[ subdomain_y ] * search_radius )
            elif subdomain_y == subdomains_y - 1:
                grid_sizes_y.append( setup[ "domain_sizesy" ] - grid_splits_y[ subdomain_y - 1 ] )
                grid_origins_y.append( setup[ "domain_origin_y" ] + grid_splits_y[ subdomain_y - 1 ] * search_radius )
                grid_index_origins_y.append( grid_splits_y[ subdomain_y - 1 ] )
                domain_sizes_y( grid_splits_y[ subdomain_y - 1 ] * search_radius )
            else:
                grid_sizes_y.append( grid_splits_y[ subdomain_y - 1 ] - grid_index_origins_y[ subdomain_y -1 ] )
                grid_origins_y.append( setup[ "domain_origin_y" ] + grid_splits_y[ subdomain_y - 1 ] * search_radius )
                grid_index_origins_y.append( grid_sizes_y[ subdomain_y - 1 ] )
                domain_sizes_y( grid_splits_y[ subdomain_y - 1 ] * search_radius )

    extra_parameters = {
        "grid_splits_x" : grid_splits_x,
        "grid_splits_y" : grid_splits_y,
        "domain_sizes_x" : domain_sizes_x,
        "domain_sizes_y" : domain_sizes_y,
        "grid_sizes_x" : grid_sizes_x,
        "grid_sizes_y" : grid_sizes_y,
        "grid_origins_x" : grid_origins_x,
        "grid_origins_y" : grid_origins_y,
        "grid_index_origins_x" : grid_index_origins_x,
        "grid_index_origins_y" : grid_index_origins_y,
    }
    setup.update( extra_parameters )
    print( "Domain decomposition informations: " )
    pprint( extra_parameters )

def split_to_subdomains( setup, fluid_rx, fluid_ry, box_rx, box_ry ):
    search_radius = setup[ "search_radius" ]
    subdomains_x = setup[ "subdomains_x" ]
    subdomains_y = setup[ "subdomains_y" ]
    grid_origins_x = setup[ "grid_origins_x" ]
    grid_origins_y = setup[ "grid_origins_y" ]

    #for subdomain_x in range( subdomains_x - 1 ):
    #    for subdomain_y in range( subdomains_y - 1 ):
    for subdomain_x in range( subdomains_x ):
        for subdomain_y in range( subdomains_y ):
            sub_fluid_rx = []; sub_fluid_ry = []
            sub_box_rx = []; sub_box_ry = []
            key_prefix = f"subdomain-x-{subdomain_x}-y-{subdomain_y}-"
            print( f"\nProcessing subdomain x: {subdomain_x}, y: {subdomain_y}:" )
            #extra_parameters = {
            #    "domain_origin_x" : domain_origin_x,
            #    "domain_origin_y" : domain_origin_y,
            #    "domain_size_x" : domain_size_x,
            #    "domain_size_y" : domain_size_y,
            #}
            #setup.update( extra_parameters )

            # load the limits of current subdomains
            if subdomains_x == 1:
                eps = 1.005
                lower_limit_x = eps * setup[ "domain_origin_x" ]
                upper_limit_x = eps * setup[ "domain_size_x" ]
            else:
                if subdomain_x == 0:
                    lower_limit_x = grid_origins_x[ subdomain_x ]
                    upper_limit_x = grid_origins_x[ subdomain_x + 1 ]
                elif subdomain_x == subdomains_x - 1:
                    lower_limit_x = grid_origins_x[ subdomain_x ] #- search_radius
                    upper_limit_x =  setup[ "domain_size_x" ]
                else:
                    lower_limit_x = grid_origins_x[ subdomain_x ]
                    upper_limit_x = grid_origins_x[ subdomain_x + 1 ]

            if subdomains_y == 1:
                eps = 1.005
                lower_limit_y = setup[ "domain_origin_y" ]
                upper_limit_y = setup[ "domain_size_y" ]
            else:
                if subdomain_y == 0:
                    lower_limit_y = grid_origins_y[ subdomain_y ]
                    upper_limit_y = grid_origins_y[ subdomain_y + 1 ]
                elif subdomain == subdomains_y - 1:
                    lower_limit_y = grid_origins_y[ subdomain_y ] - search_radius
                    upper_limit_y = setup[ "domain_origin_y" ] + ( setup[ "domain_size_y" ] + 1 ) * search_radius
                else:
                    lower_limit_y = grid_origins_y[ subdomain_y ]
                    upper_limit_y = grid_origins_y[ subdomain_y + 1 ]

            # save the part for current subdomain - fluid
            for i in range ( len( fluid_rx ) ):
                if( fluid_rx[ i ] > lower_limit_x ) and ( fluid_rx[ i ] <= upper_limit_x ) and ( fluid_ry[ i ] > lower_limit_y ) and ( fluid_ry[ i ] <= upper_limit_y ):
                    sub_fluid_rx.append( fluid_rx[ i ] )
                    sub_fluid_ry.append( fluid_ry[ i ] )

                #if( fluid_rx[ i ] > lower_limit_x and fluid_rx[ i ] <= upper_limit_x + search_radius ):
                #    sub_fluid_rx.append( fluid_rx[ i ] )
                #    sub_fluid_ry.append( fluid_ry[ i ] )

            # save the part for current subdomain - boundary
            for i in range ( len( box_rx ) ):
                if( box_rx[ i ] > lower_limit_x ) and ( box_rx[ i ] <= upper_limit_x ) and ( box_ry[ i ] > lower_limit_y ) and ( box_ry[ i ] <= upper_limit_y ):
                    sub_box_rx.append( box_rx[ i ] )
                    sub_box_ry.append( box_ry[ i ] )

                #if( box_rx[ i ] > upper_limit_x and box_rx[ i ] <= upper_limit_x + search_radius ):
                #    sub_box_rx.append( box_rx[ i ] )
                #    sub_box_ry.append( box_ry[ i ] )

            print( f"Save data for subdomain x: {subdomain_x}, y: {subdomain_y}." )
            print( f"lower_limit_x: {lower_limit_x:.2f}, upper_limit_x: {upper_limit_x:.2f}, lower_limit_y: {lower_limit_y:.2f}, upper_limit_y: {upper_limit_y:.2f}" )

            # export subdomain particles
            sub_fluid_n = len( sub_fluid_rx )
            sub_fluid_r = np.array( ( sub_fluid_rx, sub_fluid_ry, np.zeros( sub_fluid_n ) ), dtype=float ).T #!!
            sub_fluid_v = np.zeros( ( sub_fluid_n, 3 ) )
            sub_fluid_rho = setup[ "density" ] * np.ones( sub_fluid_n )
            sub_fluid_p = np.zeros( sub_fluid_n )
            sub_fluid_ptype = np.zeros( sub_fluid_n )
            sub_fluid_to_write = saveParticlesVTK.create_pointcloud_polydata(
                    sub_fluid_r, sub_fluid_v, sub_fluid_rho, sub_fluid_p, sub_fluid_ptype )
            saveParticlesVTK.save_polydata( sub_fluid_to_write, f"sources/dambreak_{key_prefix}fluid.vtk" )

            sub_box_n = len( sub_box_rx )
            sub_box_r = np.array( ( sub_box_rx, sub_box_ry, np.zeros( sub_box_n ) ), dtype=float ).T #!!
            #sub_box_ghostNodes = np.array( ( ghost_rx, ghost_ry, np.zeros( sub_box_n ) ), dtype=float ).T #!!
            sub_box_v = np.zeros( ( sub_box_n, 3 ) )
            sub_box_rho = setup[ "density" ] * np.ones( sub_box_n )
            sub_box_p = np.zeros( sub_box_n )
            sub_box_ptype = np.ones( sub_box_n )
            sub_box_to_write = saveParticlesVTK.create_pointcloud_polydata(
                    sub_box_r, sub_box_v, sub_box_rho, sub_box_p, sub_box_ptype ) # ghostNodes=boundary_ghostNodes )
            saveParticlesVTK.save_polydata( sub_box_to_write, f"sources/dambreak_{key_prefix}boundary.vtk" )

            #TODO: There I probably need to compute subdomain limits, resp. besides sub_flud_n and sub_box_n save also
            #      other domain parameters.

            setup[ f"{key_prefix}fluid_n" ] = sub_fluid_n
            setup[ f"{key_prefix}box_n" ] = sub_box_n

    print( "\n" )

def compute_domain_size( setup ):
    search_radius = setup[ "search_radius" ]

    # Resize domain by one layer of cells
    eps = 1.005
    eps_sloshing = 1.2
    domain_origin_x = eps * ( setup[ "domain_origin_x" ] - search_radius )
    domain_origin_y = eps * ( setup[ "domain_origin_y" ] - search_radius )
    domain_end_x = eps * ( setup[ "domain_end_x" ] + search_radius )
    domain_end_y = eps_sloshing * ( setup[ "domain_end_y" ] + search_radius ) #increase size in y due to sloshing
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
    config_file = config_file.replace( 'placeholderDomainOrigin-x', str( round( setup[ "domain_origin_x" ], 7 ) ) )
    config_file = config_file.replace( 'placeholderDomainOrigin-y', str( round( setup[ "domain_origin_y" ], 7 ) ) )
    config_file = config_file.replace( 'placeholderDomainSize-x', str( round( setup[ "domain_size_x" ], 7  ) ) )
    config_file = config_file.replace( 'placeholderDomainSize-y', str( round( setup[ "domain_size_y" ], 7  ) ) )

    config_file = config_file.replace( 'placeholderInitParticleDistance', str( setup[ "dp" ] ) )
    config_file = config_file.replace( 'placeholderSmoothingLength', str( round( setup[ "smoothing_length" ], 7 ) ) )
    config_file = config_file.replace( 'placeholderMass', str( setup[ "particle_mass" ] ) )
    config_file = config_file.replace( 'placeholderSpeedOfSound', str( setup[ "speed_of_sound" ] ) )
    config_file = config_file.replace( 'placeholderDensity', str( setup[ "density" ] ) )
    config_file = config_file.replace( 'placeholderTimeStep', str( setup[ "time_step" ] ) )
    config_file = config_file.replace( 'placeholderFluidParticles', str( setup[ "fluid_n" ] ) )
    config_file = config_file.replace( 'placeholderAllocatedFluidParticles', str( 2 * setup[ "fluid_n" ] ) )
    config_file = config_file.replace( 'placeholderBoundaryParticles', str( setup[ "boundary_n" ] ) )
    config_file = config_file.replace( 'placeholderAllocatedBoundaryParticles', str( setup[ "boundary_n" ] ) )

    config_file = config_file.replace( 'placeholderNumberOfSubdomains', str( setup[ "subdomains_x" ] * setup[ "subdomains_y" ]  ) )
    config_file = config_file.replace( 'placeholderSubdomains-x', str( setup[ "subdomains_x" ] ) )
    config_file = config_file.replace( 'placeholderSubdomains-y', str( setup[ "subdomains_y" ] ) )

    with open( 'sources/config.ini', 'w' ) as file:
      file.write( config_file )

def write_distributed_domain_params( setup ):
    subdomains_x = setup[ "subdomains_x" ]
    subdomains_y = setup[ "subdomains_y" ]

    # write paramerters to new created config file related to decomposition
    with open( 'sources/config-distributed-domain.ini', "w") as file:
        #file.write( f'# Distributed subdomains global informations\n' )
        #file.write( f'number-of-subdomains = { setup[ f"number_of_subdomains" ]}\n' )
        #file.write( f'subdomains-x = { subdomains_x }\n' )
        #file.write( f'subdomains-y = { subdomains_y }\n' )
        #file.write( f'domainOrigin-x = { setup[ "domain_origin_x" ]:.5}\n' )
        #file.write( f'domainOrigin-y = { setup[ "domain_origin_y" ]:.5}\n' )
        #file.write( f'domainSize-x = { setup[ "domain_size_x" ]:.5}\n' )
        #file.write( f'domainSize-y = { setup[ "domain_size_y" ]:.5}\n' )

        #file.write( f'\n' )
        file.write( f'# Subdomains informations\n' );
        for subdomain_x in range( subdomains_x ):
            for subdomain_y in range( subdomains_y ):
                key_prefix = f"subdomain-x-{subdomain_x}-y-{subdomain_y}-"
                #file.write( f'subdomain-x = { subdomain_x }\n' )
                #file.write( f'subdomain-y = { subdomain_y }\n' )
                file.write( f"{key_prefix}fluid-particles = sources/dambreak_subdomain-x-{subdomain_x}-y-{subdomain_y}-fluid.vtk\n" )
                file.write( f"{key_prefix}boundary-particles = sources/dambreak_subdomain-x-{subdomain_x}-y-{subdomain_y}-boundary.vtk\n" )
                file.write( f'{key_prefix}fluid_n = { setup[ f"{key_prefix}fluid_n" ] }\n' )
                file.write( f'{key_prefix}boundary_n = { setup[ f"{key_prefix}box_n" ] }\n' )
                file.write( f'{key_prefix}fluid_n_allocated = { 2*setup[ f"{key_prefix}fluid_n" ] }\n' )
                file.write( f'{key_prefix}boundary_n_allocated = { 2*setup[ f"{key_prefix}box_n" ] }\n' )
                subdomain_grid_origin_x = setup[ f"grid_origins_x" ][ subdomain_x ]
                subdomain_grid_origin_y = setup[ f"grid_origins_y" ][ subdomain_y ]
                #file.write( f"{key_prefix}origin-x = { subdomain_grid_origin_x:.7f}\n" )
                #file.write( f"{key_prefix}origin-y = { subdomain_grid_origin_y:.7f}\n" )
                file.write( f"{key_prefix}origin-x = { subdomain_grid_origin_x:.7f}\n" )
                file.write( f"{key_prefix}origin-y = { subdomain_grid_origin_y:.7f}\n" )
                #subdomain_grid_size_x = setup[ f"grid_sizes_x" ][ subdomain_x ]
                #subdomain_grid_size_y = setup[ f"grid_sizes_y" ][ subdomain_y ]
                #file.write( f"{key_prefix}girdSize-x = { subdomain_grid_size_x:.2f}\n" )
                #file.write( f"{key_prefix}girdSize-y = { subdomain_grid_size_y:.2f}\n" )
                subdomain_size_x = setup[ f"domain_sizes_x" ][ subdomain_x ]
                subdomain_size_y = setup[ f"domain_sizes_y" ][ subdomain_y ]
                file.write( f"{key_prefix}size-x = { subdomain_size_x:.7f}\n" )
                file.write( f"{key_prefix}size-y = { subdomain_size_y:.7f}\n" )
                file.write( f'\n' )

def configure_and_write_measuretool_parameters():
    # write parameters to config file
    with open( 'template/config-measuretool_template.ini', 'r' ) as file :
      config_file = file.read()
    with open( 'sources/config-measuretool.ini', 'w' ) as file:
      file.write( config_file )

def write_domain_background_grid( setup ):
    import domainGrid
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
                           "sources/dambreak_grid.vtk" )

if __name__ == "__main__":
    import sys
    import argparse
    import os
    from pprint import pprint

    argparser = argparse.ArgumentParser(description="Heat equation example initial condition generator")
    g = argparser.add_argument_group("distribution parameters")
    g.add_argument("--subdomains-x", type=int, default=2, help="number of subdomains in x direction")
    g.add_argument("--subdomains-y", type=int, default=1, help="number of subdomains in y direction")
    g = argparser.add_argument_group("resolution parameters")
    g.add_argument("--dp", type=float, default=0.002, help="initial distance between particles")
    g.add_argument("--h-coef", type=float, default=2, help="smoothing length coefitient")
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
        "smoothing_length" : args.h_coef * args.dp,
        "search_radius" :  round( 2 * args.h_coef * args.dp, 7 ),
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
    fluid_rx, fluid_ry = generate_dam_break_fluid_particles( dambreak_setup )
    box_rx, box_ry = generate_dam_break_boundary_particles( dambreak_setup )

    # setup parameters
    compute_domain_size( dambreak_setup )

    # split to subdomains
    generate_subdomains_data( dambreak_setup, fluid_rx, fluid_ry, box_rx, box_ry )
    split_to_subdomains( dambreak_setup, fluid_rx, fluid_ry, box_rx, box_ry )

    print("Complete example setup:")
    pprint( dambreak_setup )
    # setup parameters
    write_simulation_params( dambreak_setup )

    # write distributed domain information
    write_distributed_domain_params( dambreak_setup )

    # setup measuretool
    #configure_and_write_measuretool_parameters()

    #write linked list background grid
    write_domain_background_grid( dambreak_setup )
