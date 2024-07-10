#! /usr/bin/env python3

import math
import numpy as np
import sys
import subprocess
sys.path.append('../../../src/tools')
import saveParticlesVTK
import vtk
from vtk.numpy_interface import dataset_adapter as dsa

def generate_geometry_with_dualsphysics_gencase( dp ):
    subprocess.check_call( [ './generateGeometryWithDualSPHysicsGenCase.sh', str( dp ) ], cwd='./template/generateGeometryWithDualSPHysicsGenCase/' )

def process_dam_break_fluid_particles( setup ):
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName( f'./sources/genCaseGeometries/dambreak_fluid_dp{setup[ "dp" ]}.vtk' )
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.Update()

    polydata = reader.GetOutput()
    np_points_fluid = dsa.WrapDataObject( polydata ).Points

    fluid_n = len( np_points_fluid )
    fluid_r = np.array( np_points_fluid, dtype=float ) #!!
    fluid_v = np.array( dsa.WrapDataObject( polydata ).PointData[ 'Vel' ], dtype=float )
    fluid_rho = np.array( dsa.WrapDataObject( polydata ).PointData[ 'Rhop' ] )
    fluid_p = np.zeros( fluid_n )
    fluid_ptype = np.zeros( fluid_n )

    fluidToWrite = saveParticlesVTK.create_pointcloud_polydata( fluid_r, fluid_v, fluid_rho, fluid_p, fluid_ptype )
    saveParticlesVTK.save_polydata( fluidToWrite, "sources/dambreak_fluid.vtk" )
    setup[ "fluid_n" ] = fluid_n

    return fluid_r[ :, 0 ], fluid_r[ :, 1 ], fluid_r[ :, 2 ]

def process_dam_break_boundary_particles( setup ):
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName( f'./sources/genCaseGeometries/dambreak_bound_dp{setup[ "dp" ]}.vtk' )
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.Update()

    polydata = reader.GetOutput()
    np_points_box = dsa.WrapDataObject( polydata ).Points

    box_n = len( np_points_box )
    box_r = np.array( np_points_box, dtype=float ) #!!
    box_v = np.array( dsa.WrapDataObject( polydata ).PointData[ 'Vel' ], dtype=float )
    box_rho = np.array( dsa.WrapDataObject( polydata ).PointData[ 'Rhop' ] )
    box_p = np.zeros( box_n )
    box_ptype = np.zeros( box_n )

    boundToWrite = saveParticlesVTK.create_pointcloud_polydata( box_r, box_v, box_rho, box_p, box_ptype )
    saveParticlesVTK.save_polydata( boundToWrite, "sources/dambreak_boundary.vtk" )

    setup[ "boundary_n" ] = box_n
    setup[ "domain_origin_x" ] = min( np_points_box[ :, 0 ] )
    setup[ "domain_origin_y" ] = min( np_points_box[ :, 1 ] )
    setup[ "domain_origin_z" ] = min( np_points_box[ :, 2 ] )
    setup[ "domain_end_x" ] = max( np_points_box[ :, 0 ] )
    setup[ "domain_end_y" ] = max( np_points_box[ :, 1 ] )
    setup[ "domain_end_z" ] = max( np_points_box[ :, 2 ] )

    return box_r[ :, 0 ], box_r[ :, 1 ], box_r[ :, 2 ]

def compute_domain_size( setup ):
    search_radius = setup[ "search_radius" ]

    # Resize domain by one layer of cells
    eps = 1.005
    eps_sloshing = 1.5
    domain_origin_x = eps * ( setup[ "domain_origin_x" ] - search_radius )
    domain_origin_y = eps * ( setup[ "domain_origin_y" ] - search_radius )
    domain_origin_z = eps * ( setup[ "domain_origin_z" ] - search_radius )
    domain_end_x = eps * ( setup[ "domain_end_x" ] + search_radius )
    domain_end_y = eps * ( setup[ "domain_end_y" ] + search_radius )
    domain_end_z = eps_sloshing * ( setup[ "domain_end_z" ] + search_radius ) #increase size in z due to sloshing
    domain_size_x = domain_end_x - domain_origin_x
    domain_size_y = domain_end_y - domain_origin_y
    domain_size_z = domain_end_z - domain_origin_z

    extra_parameters = {
        "domain_origin_x" : domain_origin_x,
        "domain_origin_y" : domain_origin_y,
        "domain_origin_z" : domain_origin_z,
        "domain_size_x" : domain_size_x,
        "domain_size_y" : domain_size_y,
        "domain_size_z" : domain_size_z
    }
    setup.update( extra_parameters )

def generate_subdomains_data( setup, fluid_rx, fluid_ry, box_rx, box_ry ):
    number_of_subdomains = setup[ "number_of_subdomains" ]
    subdomains_x = setup[ "subdomains_x" ]
    subdomains_y = setup[ "subdomains_y" ]
    search_radius = setup[ "search_radius" ]
    overlap_width = setup[ "overlap_width" ]

    # referential origin (its not loaded by solver, added just for orientation)
    setup[ "referential_origin_x" ] = setup[ "domain_origin_x" ] - overlap_width * search_radius
    setup[ "referential_origin_y" ] = setup[ "domain_origin_y" ] - overlap_width * search_radius
    setup[ "referential_origin_z" ] = setup[ "domain_origin_z" ] - overlap_width * search_radius

    # initial split considering both, fluid and boundary
    #ptcs_per_subdomain = ( int )( math.ceil( setup[ "fluid_n" ] + setup[ "boundary_n" ] ) / ( setup[ "number_of_subdomains" ] ) )
    # initial split considering only fluid
    ptcs_per_subdomain = ( int )( math.ceil( setup[ "fluid_n" ] ) / ( setup[ "number_of_subdomains" ] ) )

    grid_splits_x = [ ]
    grid_splits_y = [ ]

    for subdomain_x in range( subdomains_x - 1 ):
        grid_splits_x.append( math.ceil( fluid_rx[ ptcs_per_subdomain * ( subdomain_x + 1 ) ] / search_radius ) )
    for subdomain_y in range( subdomains_y - 1 ):
        grid_splits_y.append( math.ceil( fluid_ry[ ptcs_per_subdomain * ( subdomain_y + 1 ) ] / search_radius ) )

    domain_sizes_x = []
    grid_sizes_x = []
    grid_origins_x = []
    grid_index_origins_x = []
    grid_end_x = []
    grid_index_end_x = []

    if subdomains_x == 1:
        #grid_sizes_y.append( setup[ "grid_size_y" ] ) # there is not variable grid size, only domain size
        grid_origins_x.append( setup[ "domain_origin_x" ] )
        grid_index_origins_x.append( 0 + overlap_width )
        domain_sizes_x( setup[ "domain_size_x" ] )
        grid_sizes_x.append( ( int )( np.ceil( setup[ "domain_size_x" ] / search_radius ) ) )
    else:
        #for subdomain_x in range( subdomains_x - 1 ):
        for subdomain_x in range( subdomains_x ):
            if subdomain_x == 0:
                grid_sizes_x.append( grid_splits_x[ subdomain_x ] - 0 )
                grid_origins_x.append( setup[ "domain_origin_x" ] )
                grid_index_origins_x.append( 0 + overlap_width )
                domain_sizes_x.append( grid_splits_x[ subdomain_x ] * search_radius )
                grid_end_x.append( setup[ "domain_origin_x" ] + (  grid_splits_x[ subdomain_x ] * search_radius ) )
                grid_index_end_x.append( grid_splits_x[ subdomain_x ] )
            elif subdomain_x == subdomains_x - 1:
                #grid_sizes_x.append( setup[ "domain_size_x" ] - grid_splits_x[ subdomain_x - 1 ] ) //TODO:
                grid_sizes_x.append( math.ceil( setup[ "domain_size_x" ] / search_radius ) - grid_splits_x[ subdomain_x - 1 ] )
                #grid_origins_x.append( setup[ "domain_origin_x" ] + grid_splits_x[ subdomain_x - 1 ] * search_radius )
                #grid_index_origins_x.append( grid_splits_x[ subdomain_x - 1 ] * search_radius )
                #FIXME: grid_origins_x.append( setup[ "domain_origin_x" ] + search_radius * np.sum( grid_sizes_x[ 0 : subdomain_x ] ) ) #TODO: Use this, added domain origin
                grid_origins_x.append( setup[ "domain_origin_x" ] + (  grid_splits_x[ subdomain_x - 1 ] * search_radius ) )
                print( f"""Adding grid index origin for submodule x: {subdomain_x}, located at: { setup[ 'domain_origin_x' ] + search_radius * np.sum( grid_sizes_x[ 0 : subdomain_x ] ) },
                      taking into grid_sizes_x[ 0 : subdomain_x ]: {np.sum( grid_sizes_x[ 0 : subdomain_x ] )}""" )
                grid_index_origins_x.append( np.sum( grid_sizes_x[ 0 : subdomain_x ] ) + overlap_width ) #TODO: Use this.
                domain_sizes_x.append( setup[ "domain_size_x" ]  - grid_splits_x[ subdomain_x - 1 ] * search_radius )
                grid_end_x.append( setup[ "domain_size_x" ] )
                grid_index_end_x.append( np.sum( grid_sizes_x[ 0 : subdomain_x ] ) + math.ceil( setup[ "domain_size_x" ] / search_radius ) - grid_splits_x[ subdomain_x - 1 ] )
            else:
                ##grid_sizes_x.append( grid_splits_x[ subdomain_x - 1 ] - grid_index_origins_x[ subdomain_x -1 ] )
                grid_sizes_x.append( grid_splits_x[ subdomain_x ] - grid_splits_x[ subdomain_x - 1 ] )
                grid_origins_x.append( setup[ "domain_origin_x" ] + grid_splits_x[ subdomain_x - 1 ] * search_radius )
                print( f"""Adding grid index origin for submodule x: {subdomain_x}, located at: {setup[ 'domain_origin_x' ] + grid_splits_x[ subdomain_x - 1 ] * search_radius},
                      taking into accout grid_split: {grid_splits_x[ subdomain_x - 1 ]}""" )
                #grid_index_origins_x.append( grid_sizes_x[ subdomain_x - 1 ] )
                grid_index_origins_x.append( np.sum( grid_sizes_x[ 0 : subdomain_x ] ) + overlap_width )
                #domain_sizes_x.append( grid_splits_x[ subdomain_x - 1 ] * search_radius )
                domain_sizes_x.append( ( grid_splits_x[ subdomain_x ] - grid_splits_x[ subdomain_x - 1 ] ) * search_radius )
                grid_end_x.append( setup[ "domain_origin_x" ] + (  grid_splits_x[ subdomain_x ] * search_radius ) )
                grid_index_end_x.append( grid_splits_x[ subdomain_x ] )

    domain_sizes_y = []
    grid_sizes_y = []
    grid_origins_y = []
    grid_index_origins_y = []

    if subdomains_y == 1:
        #grid_sizes_y.append( setup[ "grid_size_y" ] ) # there is not variable grid size, only domain size
        grid_origins_y.append( setup[ "domain_origin_y" ] )
        grid_index_origins_y.append( 0 + overlap_width )
        domain_sizes_y.append( setup[ "domain_size_y" ] )
        grid_sizes_y.append( ( int )( np.ceil( setup[ "domain_size_y" ] / search_radius ) ) )
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
        "gird_ends_x" : grid_end_x,
        "gird_index_end_x" : grid_index_end_x,
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
            sub_fluid_rx = []; sub_fluid_ry = []; sub_fluid_rz = []
            sub_box_rx = []; sub_box_ry = []; sub_box_rz = []
            key_prefix = f"subdomain-x-{subdomain_x}-y-{subdomain_y}-"
            print( f"\nProcessing subdomain x: {subdomain_x}, y: {subdomain_y}:" )

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
                    sub_fluid_rz.append( fluid_rz[ i ] )

                #if( fluid_rx[ i ] > lower_limit_x and fluid_rx[ i ] <= upper_limit_x + search_radius ):
                #    sub_fluid_rx.append( fluid_rx[ i ] )
                #    sub_fluid_ry.append( fluid_ry[ i ] )

            # save the part for current subdomain - boundary
            for i in range ( len( box_rx ) ):
                if( box_rx[ i ] > lower_limit_x ) and ( box_rx[ i ] <= upper_limit_x ) and ( box_ry[ i ] > lower_limit_y ) and ( box_ry[ i ] <= upper_limit_y ):
                    sub_box_rx.append( box_rx[ i ] )
                    sub_box_ry.append( box_ry[ i ] )
                    sub_box_rz.append( box_rz[ i ] )

                #if( box_rx[ i ] > upper_limit_x and box_rx[ i ] <= upper_limit_x + search_radius ):
                #    sub_box_rx.append( box_rx[ i ] )
                #    sub_box_ry.append( box_ry[ i ] )

            print( f"Save data for subdomain x: {subdomain_x}, y: {subdomain_y}." )
            print( f"lower_limit_x: {lower_limit_x:.2f}, upper_limit_x: {upper_limit_x:.2f}, lower_limit_y: {lower_limit_y:.2f}, upper_limit_y: {upper_limit_y:.2f}" )

            # export subdomain particles
            sub_fluid_n = len( sub_fluid_rx )
            sub_fluid_r = np.array( ( sub_fluid_rx, sub_fluid_ry, sub_fluid_rz ), dtype=float ).T #!!
            sub_fluid_v = np.zeros( ( sub_fluid_n, 3 ) )
            sub_fluid_rho = setup[ "density" ] * np.ones( sub_fluid_n )
            sub_fluid_p = np.zeros( sub_fluid_n )
            sub_fluid_ptype = np.zeros( sub_fluid_n )
            sub_fluid_to_write = saveParticlesVTK.create_pointcloud_polydata(
                    sub_fluid_r, sub_fluid_v, sub_fluid_rho, sub_fluid_p, sub_fluid_ptype )
            saveParticlesVTK.save_polydata( sub_fluid_to_write, f"sources/dambreak_{key_prefix}fluid.vtk" )

            sub_box_n = len( sub_box_rx )
            sub_box_r = np.array( ( sub_box_rx, sub_box_ry, sub_box_rz ), dtype=float ).T #!!
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
    config_file = config_file.replace( 'placeholderTimeStep', str( round( setup[ "time_step" ], 8 ) ) )
    config_file = config_file.replace( 'placeholderFluidParticles', str( setup[ "fluid_n" ] ) )
    config_file = config_file.replace( 'placeholderAllocatedFluidParticles', str( setup[ "fluid_n" ] ) )
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
                file.write( f'{key_prefix}boundary_n_allocated = { 3*setup[ f"{key_prefix}box_n" ] }\n' )
                subdomain_grid_origin_x = setup[ f"grid_origins_x" ][ subdomain_x ]
                subdomain_grid_origin_y = setup[ f"grid_origins_y" ][ subdomain_y ]
                file.write( f"{key_prefix}origin-x = { subdomain_grid_origin_x:.7f}\n" )
                file.write( f"{key_prefix}origin-y = { subdomain_grid_origin_y:.7f}\n" )
                file.write( f"{key_prefix}origin-z = { setup[ 'domain_origin_z' ]:.7f}\n" ) #2D decomposition
                subdomain_grid_origin_glob_coords_x = setup[ f"grid_index_origins_x" ][ subdomain_x ]
                subdomain_grid_origin_glob_coords_y = setup[ f"grid_index_origins_y" ][ subdomain_y ]
                file.write( f"{key_prefix}origin-global-coords-x = { subdomain_grid_origin_glob_coords_x }\n" )
                file.write( f"{key_prefix}origin-global-coords-y = { subdomain_grid_origin_glob_coords_y }\n" )
                file.write( f"{key_prefix}origin-global-coords-z = { 0 + setup[ 'overlap_width' ]}\n" ) #2D decomposition
                subdomain_size_x = setup[ f"domain_sizes_x" ][ subdomain_x ]
                subdomain_size_y = setup[ f"domain_sizes_y" ][ subdomain_y ]
                file.write( f"{key_prefix}size-x = { subdomain_size_x:.7f}\n" )
                file.write( f"{key_prefix}size-y = { subdomain_size_y:.7f}\n" )
                file.write( f"{key_prefix}size-z = { setup[ 'domain_size_z' ]:.7f}\n" ) #2D decomposition
                subdomain_grid_dims_x = setup[ f"grid_sizes_x" ][ subdomain_x ]
                subdomain_grid_dims_y = setup[ f"grid_sizes_y" ][ subdomain_y ]
                file.write( f"{key_prefix}grid-dimensions-x = { subdomain_grid_dims_x }\n" )
                file.write( f"{key_prefix}grid-dimensions-y = { subdomain_grid_dims_y }\n" )
                file.write( f"{key_prefix}grid-dimensions-z = { ( int )( np.ceil( setup[ 'domain_size_z'] / setup[ 'search_radius' ] ) ) }\n" ) #2D decomposition
                file.write( f'\n' )

def write_domain_background_grid( setup ):
    import domainGrid
    #TODO: Rename the DomainGrid function
    search_radius = setup[ "search_radius" ]
    grid_size_x = round( setup[ "domain_size_x" ] / search_radius )
    grid_size_y = round( setup[ "domain_size_y" ] / search_radius )
    grid_size_z = round( setup[ "domain_size_z" ] / search_radius )
    grid_origin_x = setup[ "domain_origin_x" ]
    grid_origin_y = setup[ "domain_origin_y" ]
    grid_origin_z = setup[ "domain_origin_z" ]
    grid_sectors = np.zeros( grid_size_x * grid_size_y * grid_size_z )

    domainGrid.DomainGrid3D( grid_size_x,
                           grid_size_y,
                           grid_size_z,
                           grid_origin_x,
                           grid_origin_y,
                           grid_origin_z,
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
    g.add_argument("--subdomains-x", type=int, default=3, help="number of subdomains in x direction")
    g.add_argument("--subdomains-y", type=int, default=1, help="number of subdomains in y direction")
    g.add_argument("--overlap-width", type=int, default=1, help="width of domain overlap in cells")
    g = argparser.add_argument_group("resolution parameters")
    g.add_argument("--dp", type=float, default=0.02, help="initial distance between particles")
    g.add_argument("--h-coef", type=float, default=2, help="smoothing length coefficient")
    g = argparser.add_argument_group("simulation parameters")
    g.add_argument("--density", type=float, default=1000, help="referential density of the fluid")
    g.add_argument("--speed-of-sound", type=float, default=45.17, help="speed of sound")
    g.add_argument("--cfl", type=float, default=0.15, help="referential density of the fluid")
    g = argparser.add_argument_group("control initialization")
    g.add_argument( '--generate-geometry', default=True, action=argparse.BooleanOptionalAction, help="generate new geometry with gencase" )

    args = argparser.parse_args()

    dambreak_setup = {
        # general parameteres
        "dp" : args.dp,
        "h_coef" : args.h_coef,
        "density" : args.density,
        "speed_of_sound" : args.speed_of_sound,
        "cfl" : args.cfl,
        "particle_mass" : args.density * ( args.dp * args.dp * args.dp ),
        "smoothing_length" : args.h_coef * args.dp,
        "search_radius" :  2 * args.h_coef * args.dp,
        "time_step" : args.cfl * ( args.h_coef * args.dp ) / args.speed_of_sound,
        # distribution parameters
        "subdomains_x" : args.subdomains_x,
        "subdomains_y" : args.subdomains_y,
        "number_of_subdomains" : args.subdomains_x * args.subdomains_y,
        "overlap_width" : args.overlap_width
    }

    # create necessary folders
    resultsPath = r'./results'
    if not os.path.exists( resultsPath ):
        os.makedirs( resultsPath )

    sourcesPath = r'./sources'
    if not os.path.exists( sourcesPath ):
        os.makedirs( sourcesPath )

    # generate particles using DualSPHysics genCase
    if args.generate_geometry:
        generate_geometry_with_dualsphysics_gencase( dambreak_setup[ "dp" ] )

    # generate particles
    fluid_rx, fluid_ry, fluid_rz = process_dam_break_fluid_particles( dambreak_setup )
    box_rx, box_ry, box_rz = process_dam_break_boundary_particles( dambreak_setup )

    # setup parameters
    compute_domain_size( dambreak_setup )

    # split to subdomains
    generate_subdomains_data( dambreak_setup, fluid_rx, fluid_ry, box_rx, box_ry )
    split_to_subdomains( dambreak_setup, fluid_rx, fluid_ry, box_rx, box_ry )

    # write simulation params
    pprint( dambreak_setup )
    write_simulation_params( dambreak_setup )

    # write distributed domain information
    write_distributed_domain_params( dambreak_setup )

    #write linked list background grid
    write_domain_background_grid( dambreak_setup )
