#! /usr/bin/env python3

import numpy as np
import math
import sys
sys.path.append('../../../src/tools')
import saveParticlesVTK

def compute_hydrostatic_density( ry, fluid_height, rho0, speed_of_sound ):
    hydrostaticPressure = rho0 * 9.81 * ( fluid_height - ry )
    hydrostaticDensity = ( ( hydrostaticPressure / ( speed_of_sound** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0
    #hydrostaticDensity = rho0 + hydrostaticPressure / speed_of_sound**2
    return hydrostaticDensity

def generate_dam_break_fluid_particles( setup ):
    number_of_subdomains = setup[ 'number_of_subdomains' ]
    for subdomain in range( 0, number_of_subdomains ):
        subdomain_key = f"subdomain-{subdomain}"

        fluid_rx = []; fluid_ry = []
        fluid_density = []
        dp = setup[ 'dp' ] * setup[ f"{subdomain_key}-factor" ]
        fluid_length = setup[ 'fluid_length' ]
        fluid_height = setup[ 'fluid_height' ]
        fluid_lenght_n = round( fluid_length / dp )
        fluid_height_n = round( fluid_height / dp )
        rho0 = setup[ 'density' ]
        speed_of_sound = setup[ 'speed_of_sound' ]

        #x_min = setup[ subdomain_key + 'x_min' ]
        #x_max = setup[ subdomain_key + 'x_max' ]
        #y_min = setup[ subdomain_key + 'y_min' ]
        #y_max = setup[ subdomain_key + 'y_max' ]

        search_radius_L0 = setup[ "search_radius" ]
        grid_origin_x = setup[ "domain_origin_x" ] +  search_radius_L0 * setup[ f'{subdomain_key}-origin_glob_coords_x' ]
        grid_size_x = setup[ f'{subdomain_key}-dimensions_x' ] *  setup[ f'{subdomain_key}-search_radius' ]
        grid_origin_y = setup[ "domain_origin_y" ] +  search_radius_L0 * setup[ f'{subdomain_key}-origin_glob_coords_y' ]
        grid_size_y = setup[ f'{subdomain_key}-dimensions_y' ] *  setup[ f'{subdomain_key}-search_radius' ]
        x_min = grid_origin_x
        x_max = grid_origin_x + grid_size_x
        y_min = grid_origin_y
        y_max = grid_origin_y + grid_size_y
        print( f"Generate fluid particles - subdomain limits: \nx: [ {x_min}, {x_max} ], y: [ {y_min}, {y_max} ]" )

        for x in range( fluid_lenght_n ):
            for y in range( fluid_height_n ):
                rx = dp * ( x + 1 )
                ry = dp * ( y + 1 )
                if ( rx > x_min and rx <= x_max ) and ( ry > y_min and ry <= y_max ):
                    fluid_rx.append( rx )
                    fluid_ry.append( ry )
                    fluid_density.append( compute_hydrostatic_density( fluid_ry[ -1 ], fluid_height, rho0, speed_of_sound ) )

        fluid_n = len( fluid_rx )
        fluid_r = np.array( ( fluid_rx, fluid_ry, np.zeros( fluid_n ) ), dtype=float ).T #!!
        fluid_v = np.zeros( ( fluid_n, 3 ) )
        fluid_rho = np.array( fluid_density, dtype=float )
        fluid_p = np.zeros( fluid_n )
        fluid_ptype = np.zeros( fluid_n )
        fluid_to_write = saveParticlesVTK.create_pointcloud_polydata( fluid_r, fluid_v, fluid_rho, fluid_p, fluid_ptype )
        saveParticlesVTK.save_polydata( fluid_to_write, f"sources/{subdomain_key}-dambreak_fluid.vtk" )

        # compute potential energy
        mass = rho0 * ( dp * dp )
        Epot0 = mass * 9.81 * np.sum( fluid_ry )
        print( f"Initial potential energy of fluid Epot0: {Epot0}" )
        setup[ f"{subdomain_key}-fluid_n" ] = fluid_n

def generate_dam_break_boundary_particles( setup ):
    number_of_subdomains = setup[ 'number_of_subdomains' ]
    for subdomain in range( 0, number_of_subdomains ):
        subdomain_key = f"subdomain-{subdomain}"

        box_rx = []; box_ry = []
        ghost_rx = []; ghost_ry = []
        normal_x = []; normal_y = []
        dp = setup[ 'dp' ] * setup[ f"{subdomain_key}-factor" ]
        n_boundary_layers = setup[ 'n_boundary_layers' ]
        box_length_n = round( setup[ 'box_length' ] / dp )
        box_height_n = round( setup[ 'box_height' ] / dp )
        rho0 = setup[ 'density' ]

        search_radius_L0 = setup[ "search_radius" ]
        grid_origin_x = setup[ "domain_origin_x" ] +  search_radius_L0 * setup[ f'{subdomain_key}-origin_glob_coords_x' ]
        grid_size_x = setup[ f'{subdomain_key}-dimensions_x' ] *  setup[ f'{subdomain_key}-search_radius' ]
        grid_origin_y = setup[ "domain_origin_y" ] +  search_radius_L0 * setup[ f'{subdomain_key}-origin_glob_coords_y' ]
        grid_size_y = setup[ f'{subdomain_key}-dimensions_y' ] *  setup[ f'{subdomain_key}-search_radius' ]
        x_min = grid_origin_x
        x_max = grid_origin_x + grid_size_x
        y_min = grid_origin_y
        y_max = grid_origin_y + grid_size_y

        # left wall
        for layer in range( n_boundary_layers ):
            for y in range( box_height_n - 1 ):
                rx = 0. - layer * dp
                ry = ( y + 1 ) * dp
                if ( rx > x_min and rx <= x_max ) and ( ry > y_min and ry <= y_max ):
                    box_rx.append( 0. - layer * dp )
                    box_ry.append( ( y + 1 ) * dp )
                    ghost_rx.append( 0. + dp * ( layer + 1 ) )
                    ghost_ry.append( ( y + 1 ) * dp )
                    normal_x.append( 1. )
                    normal_y.append( 0. )

        # bottom wall
        for layer in range( n_boundary_layers ):
            for x in range( box_length_n - n_boundary_layers + 1 ):
                rx = ( x + 1 ) * dp
                ry = 0. - layer * dp
                if ( rx > x_min and rx <= x_max ) and ( ry > y_min and ry <= y_max ):
                    box_rx.append( rx )
                    box_ry.append( ry )
                    ghost_rx.append( ( x + 1 ) * dp )
                    ghost_ry.append( 0. + dp * ( layer + 1 ) )
                    normal_x.append( 0. )
                    normal_y.append( 1. )

        x_last = box_rx[ -1 ] + dp #due to discretisation, we need to save last value of bottom wall

        # right wall
        for layer in range( n_boundary_layers ):
            for y in range( box_height_n - 1 ):
                rx = x_last + dp * layer
                ry = ( y + 1 ) * dp
                if ( rx > x_min and rx <= x_max ) and ( ry > y_min and ry <= y_max ):
                    box_rx.append( rx )
                    box_ry.append( ry )
                    ghost_rx.append( x_last - dp * ( layer + 1 ) )
                    ghost_ry.append( ( y + 1 ) * dp )
                    normal_x.append( -1. )
                    normal_y.append( 0. )

        # generate the corners
        def generate90degCorner( x, y, dirx, diry ):
                for layer in range( n_boundary_layers ):
                    for k in range( n_boundary_layers ):
                        rx = x + k * dp * dirx
                        ry = y + layer * dp * diry
                        if ( rx > x_min and rx <= x_max ) and ( ry > y_min and ry <= y_max ):
                            box_rx.append( rx )
                            box_ry.append( ry )
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
        saveParticlesVTK.save_polydata( box_to_write, f"sources/{subdomain_key}-dambreak_boundary.vtk" )

        setup[ f"{subdomain_key}-boundary_n" ] = boundary_n
        #setup[ subdomain_key + "origin_x" ] = min( box_rx )
        #setup[ subdomain_key + "origin_y" ]  = min( box_ry )
        #setup[ subdomain_key + "end_x" ] = max( box_rx )
        #setup[ subdomain_key + "end_y" ] = max( box_ry )

#def compute_domain_size( setup ):
#    search_radius = setup[ "search_radius" ]
#    # Resize domain by one layer of cells
#    eps = 1.005
#    eps_sloshing = 1.2
#    domain_origin_x = eps * ( setup[ "domain_origin_x" ] - search_radius )
#    domain_origin_y = eps * ( setup[ "domain_origin_y" ] - search_radius )
#    domain_end_x = eps * ( setup[ "domain_end_x" ] + search_radius )
#    domain_end_y = eps_sloshing * ( setup[ "domain_end_y" ] + search_radius )
#    domain_size_x = domain_end_x - domain_origin_x
#    domain_size_y = domain_end_y - domain_origin_y
#
#    extra_parameters = {
#        "domain_origin_x" : domain_origin_x,
#        "domain_origin_y" : domain_origin_y,
#        "domain_size_x" : domain_size_x,
#        "domain_size_y" : domain_size_y,
#    }
#    setup.update( extra_parameters )


def compute_subdomain_size( setup ):
    # load coordinates of global domain
    search_radius_L0 = setup[ "search_radius" ]
    glob_grid_origin_x = setup[ "domain_origin_x" ]
    glob_grid_origin_y = setup[ "domain_origin_y" ]
    glob_grid_dimension_x = int( math.ceil( setup[ "domain_size_x" ] / search_radius_L0 ) )
    glob_grid_dimension_y = int( math.ceil( setup[ "domain_size_x" ] / search_radius_L0 ) )

    # compute coords of subdomains
    number_of_subdomains = setup[ 'number_of_subdomains' ]
    for subdomain in range( 0, number_of_subdomains ):
        subdomain_key = f"subdomain-{subdomain}"
        subdomain_x_min = setup[ f"{subdomain_key}-x_min" ]
        subdomain_x_max = setup[ f"{subdomain_key}-x_max" ]
        subdomain_y_min = setup[ f"{subdomain_key}-y_min" ]
        subdomain_y_max = setup[ f"{subdomain_key}-y_max" ]

        subdomain_factor = setup[ f"{subdomain_key}-factor" ]
        subdomain_search_radius = subdomain_factor * search_radius_L0

        #: x dimension
        if subdomain_x_min < glob_grid_origin_x:
            subdomain_origin_x = glob_grid_origin_x
            subdomain_grid_origin_glob_coords_x = 0
        else:
            #ASSERT x_min > domain_end
            subdomain_origin_x = subdomain_x_min
            subdomain_grid_origin_glob_coords_x = int( math.floor( ( subdomain_origin_x - glob_grid_origin_x ) / search_radius_L0 ) )

        if subdomain_x_max > glob_grid_origin_x + glob_grid_dimension_x * search_radius_L0:
            subdomain_glob_dimension_x = glob_grid_dimension_x - subdomain_grid_origin_glob_coords_x
            subdomain_dimension_x = subdomain_glob_dimension_x * int( 1 / subdomain_factor )
        else:
            #ASSERT x_max > domain_end
            subdomain_end_x = subdomain_x_max
            subdomain_grid_end_global_cords_x = int( math.floor( ( subdomain_end_x - glob_grid_origin_x ) / search_radius_L0 ) )
            subdomain_glob_dimension_x = subdomain_grid_end_global_cords_x - subdomain_grid_origin_glob_coords_x
            subdomain_dimension_x = subdomain_glob_dimension_x * int( 1 / subdomain_factor )

        #: y dimension
        if subdomain_y_min < glob_grid_dimension_y:
            subdomain_origin_y = glob_grid_origin_y
            subdomain_grid_origin_glob_coords_y = 0
        else:
            #ASSERT y_min > domain_end
            subdomain_origin_y = subdomain_y_min
            subdomain_grid_origin_glob_coords_y = int( math.floor( ( subdomain_origin_y - glob_grid_origin_y ) / search_radius_L0 ) )

        if subdomain_y_max > glob_grid_origin_y + glob_grid_dimension_y * search_radius_L0:
            subdomain_glob_dimension_y = glob_grid_dimension_y - subdomain_grid_origin_glob_coords_y
            subdomain_dimension_y = subdomain_glob_dimension_y * int( 1 / subdomain_factor )
        else:
            #ASSERT y_max > domain_end
            subdomain_end_y = subdomain_y_max
            subdomain_grid_end_global_cords_y = int( math.floor( ( subdomain_end_y - glob_grid_origin_y ) / search_radius_L0 ) )
            subdomain_glob_dimension_y = subdomain_grid_end_global_cords_y - subdomain_grid_origin_glob_coords_y
            subdomain_dimension_y = subdomain_glob_dimension_y * int( 1 / subdomain_factor )

        subdomain_params = {
            f"{subdomain_key}-search_radius" : subdomain_search_radius,
            f"{subdomain_key}-origin_glob_coords_x" : subdomain_grid_origin_glob_coords_x,
            f"{subdomain_key}-origin_glob_coords_y" : subdomain_grid_origin_glob_coords_y,
            f"{subdomain_key}-dimensions_x" : subdomain_dimension_x,
            f"{subdomain_key}-dimensions_y" : subdomain_dimension_y,
        }
        setup.update( subdomain_params )


   #     if( subdomain_grid_origin_x < domain_origin_x ):
   #         domain_origin_x =  setup[ subdomain_key + "origin_x" ]
   #     if( setup[ subdomain_key + "origin_y" ] < domain_origin_y ):
   #         domain_origin_y =  setup[ subdomain_key + "origin_y" ]
   #     if( setup[ subdomain_key + "end_x" ] < domain_end_x ):
   #         domain_end_x =  setup[ subdomain_key + "end_x" ]
   #     if( setup[ subdomain_key + "end_y" ] < domain_end_y ):
   #         domain_end_y =  setup[ subdomain_key + "end_y" ]

   # #: # Get global domains limit
   # #: fluid_n = setup[ 'subdomain-0-fluid_n' ]
   # #: boundary_n = setup[ 'subdomain-0-boundary_n' ]
   # #: number_of_subdomains = setup[ 'number_of_subdomains' ]
   # #: domain_origin_x = setup[ 'subdomain-0-origin_x' ]
   # #: domain_origin_y = setup[ 'subdomain-0-origin_y' ]
   # #: domain_end_x = setup[ 'subdomain-0-end_x' ]
   # #: domain_end_y = setup[ 'subdomain-0-end_y' ]

   # for subdomain in range( 1, number_of_subdomains ):
   #     subdomain_key = f"subdomain-{subdomain}-"
   #     fluid_n += setup[ subdomain_key + 'fluid_n' ]
   #     boundary_n += setup[ subdomain_key + 'boundary_n' ]

   #     if( setup[ subdomain_key + "origin_x" ] < domain_origin_x ):
   #         domain_origin_x =  setup[ subdomain_key + "origin_x" ]
   #     if( setup[ subdomain_key + "origin_y" ] < domain_origin_y ):
   #         domain_origin_y =  setup[ subdomain_key + "origin_y" ]
   #     if( setup[ subdomain_key + "end_x" ] < domain_end_x ):
   #         domain_end_x =  setup[ subdomain_key + "end_x" ]
   #     if( setup[ subdomain_key + "end_y" ] < domain_end_y ):
   #         domain_end_y =  setup[ subdomain_key + "end_y" ]

   # # Compute origin grid coordinates
   # search_radius = setup[ 'search_radius' ]
   # for subdomain in range( 0, number_of_subdomains ):
   #     subdomain_key = f"subdomain-{subdomain}-"
   #     search_radius = setup[ f"search_radius" ]

   #     setup[ f"{subdomain_key}origin-global-coords-x" ] = int( np.floor( ( setup[ subdomain_key + "origin_x" ] - domain_origin_x ) / search_radius ) )
   #     setup[ f"{subdomain_key}origin-global-coords-y" ] = int( np.floor( ( setup[ subdomain_key + "origin_y" ] - domain_origin_y ) / search_radius ) )
   #     setup[ f"{subdomain_key}grid-dimensions-x" ] = int( 1 / ( setup[ f"{subdomain_key}factor" ] ) * ( setup[ subdomain_key + "end_x" ] - setup[ subdomain_key + "origin_x" ] ) / search_radius )
   #     setup[ f"{subdomain_key}grid-dimensions-y" ] = int( 1 / ( setup[ f"{subdomain_key}factor" ] ) * ( setup[ subdomain_key + "end_y" ] - setup[ subdomain_key + "origin_y" ] ) / search_radius )
   #     setup[ f"{subdomain_key}size-x" ] = setup[ subdomain_key + "end_x" ] - setup[ subdomain_key + "origin_x" ]
   #     setup[ f"{subdomain_key}size-y" ] = setup[ subdomain_key + "end_y" ] - setup[ subdomain_key + "origin_y" ]
   #     domain_size_x = domain_end_x - domain_origin_x
   #     domain_size_y = domain_end_y - domain_origin_y

   # # Save global domain limits
   # domain_parameters = {
   #     "fluid_n" : fluid_n,
   #     "boundary_n" : boundary_n,
   #     "domain_origin_x" : domain_origin_x,
   #     "domain_origin_y" : domain_origin_y,
   #     "domain_end_x" : domain_end_x,
   #     "domain_end_y" : domain_end_y,
   # }
   # setup.update( domain_parameters )



def write_distributed_domain_params( setup ):
    # write paramerters to new created config file related to decomposition
    with open( 'sources/config-distributed-domain.ini', "w") as file:

        file.write( f'# Subdomains informations\n' );
        for subdomain in range( 0, setup[ 'number_of_subdomains' ] ):
            key_prefix = f"subdomain-{subdomain}-"
            #file.write( f'subdomain-x = { subdomain_x }\n' )
            #file.write( f'subdomain-y = { subdomain_y }\n' )
            file.write( f"{key_prefix}fluid-particles = sources/subdomain-{subdomain}-dambreak_fluid.vtk\n" )
            file.write( f"{key_prefix}boundary-particles = sources/subdomain-{subdomain}-dambreak_boundary.vtk\n" )
            file.write( f'{key_prefix}fluid_n = { setup[ f"{key_prefix}fluid_n" ] }\n' )
            file.write( f'{key_prefix}boundary_n = { setup[ f"{key_prefix}boundary_n" ] }\n' )
            file.write( f'{key_prefix}fluid_n_allocated = { 2*setup[ f"{key_prefix}fluid_n" ] }\n' )
            file.write( f'{key_prefix}boundary_n_allocated = { 3*setup[ f"{key_prefix}boundary_n" ] }\n' )
            file.write( f'{key_prefix}refinement-factor = { setup[ f"{key_prefix}factor" ] }\n' )
            search_radius = setup[ f"{key_prefix}search_radius" ]
            file.write( f'{key_prefix}refinement-factor = { setup[ f"{key_prefix}factor" ] }\n' )

            subdomain_grid_origin_glob_coords_x = setup[ f"{key_prefix}origin_glob_coords_x" ]
            subdomain_grid_origin_glob_coords_y = setup[ f"{key_prefix}origin_glob_coords_y" ]
            file.write( f"{key_prefix}origin-global-coords-x = { subdomain_grid_origin_glob_coords_x }\n" )
            file.write( f"{key_prefix}origin-global-coords-y = { subdomain_grid_origin_glob_coords_y }\n" )
            subdomain_grid_dims_x = setup[ f"{key_prefix}dimensions_x" ]
            subdomain_grid_dims_y = setup[ f"{key_prefix}dimensions_y" ]
            file.write( f"{key_prefix}grid-dimensions-x = { subdomain_grid_dims_x }\n" )
            file.write( f"{key_prefix}grid-dimensions-y = { subdomain_grid_dims_y }\n" )

            # useless
            search_radius_L0 = setup[ "search_radius" ]
            subdomain_grid_origin_x = setup[ f"domain_origin_x" ] + search_radius_L0 * subdomain_grid_origin_glob_coords_x
            subdomain_grid_origin_y = setup[ f"domain_origin_y" ] + search_radius_L0 * subdomain_grid_origin_glob_coords_y
            file.write( f"{key_prefix}origin-x = { subdomain_grid_origin_x:.7f}\n" )
            file.write( f"{key_prefix}origin-y = { subdomain_grid_origin_y:.7f}\n" )
            subdomain_size_x = search_radius * subdomain_grid_dims_x
            subdomain_size_y = search_radius * subdomain_grid_dims_y
            file.write( f"{key_prefix}size-x = { subdomain_size_x:.7f}\n" )
            file.write( f"{key_prefix}size-y = { subdomain_size_y:.7f}\n" )

            file.write( f'\n' )

def write_simulation_params( setup ):
    # write parameters to config file
    with open( 'template/config_template.ini', 'r' ) as file :
      config_file = file.read()

    config_file = config_file.replace( 'placeholderSearchRadius', f'{setup[ "search_radius" ] }' )
    config_file = config_file.replace( 'placeholderDomainOrigin-x', f'{setup[ "domain_origin_x" ]:.7f}' )
    config_file = config_file.replace( 'placeholderDomainOrigin-y', f'{setup[ "domain_origin_y" ]:.7f}' )
    config_file = config_file.replace( 'placeholderDomainSize-x', f'{setup[ "domain_size_x" ]:.7f}' )
    config_file = config_file.replace( 'placeholderDomainSize-y', f'{setup[ "domain_size_y" ]:.7f}' )

    config_file = config_file.replace( 'placeholderInitParticleDistance', f'{ setup[ "dp" ] }' )
    config_file = config_file.replace( 'placeholderSmoothingLength', f'{ setup[ "smoothing_length" ] }' )
    config_file = config_file.replace( 'placeholderMass', f'{ setup[ "particle_mass" ] }' )
    config_file = config_file.replace( 'placeholderSpeedOfSound', f'{ setup[ "speed_of_sound" ] }' )
    config_file = config_file.replace( 'placeholderDensity', f'{ setup[ "density" ] }' )
    config_file = config_file.replace( 'placeholderTimeStep', f'{ setup[ "time_step" ] }' )
    config_file = config_file.replace( 'placeholderCFL', f'{ setup[ "cfl" ] }' )
    config_file = config_file.replace( 'placeholderAlpha', f'{ setup[ "alpha" ] }' )
    config_file = config_file.replace( 'placeholderDynamicVicosity', f'{ setup[ "dynamic_viscosity" ] }' )
    config_file = config_file.replace( 'placeholderFluidParticles', f'{ setup[ "fluid_n" ] }' )
    config_file = config_file.replace( 'placeholderAllocatedFluidParticles', f'{ setup[ "fluid_n" ] }' )
    config_file = config_file.replace( 'placeholderBoundaryParticles', f'{ setup[ "boundary_n" ] }' )
    config_file = config_file.replace( 'placeholderAllocatedBoundaryParticles', f'{ setup[ "boundary_n" ] }' )

    with open( 'sources/config.ini', 'w' ) as file:
      file.write( config_file )

    # write parameters to config header file
    with open( 'template/config_template.h', 'r' ) as file :
      config_file = file.read()

    config_file = config_file.replace( '#placeholderBoundaryConditionsType',  setup[ "bc_type" ] )
    config_file = config_file.replace( '#placeholderDiffusiveTerm', setup[ "diffusive_term" ] )
    config_file = config_file.replace( '#placeholderViscosTerm', setup[ "viscous_term" ] )

    with open( 'template/config.h', 'w' ) as file:
      file.write( config_file )

def save_grid( setup ):

    import domainGrid

    for subdomain in range( 0, setup[ 'number_of_subdomains' ] ):
        subdomain_key = f"subdomain-{subdomain}"
        search_radius_L0 = setup[ "search_radius" ]
        domainGrid.domainGrid( setup[ f"{subdomain_key}-dimensions_x" ],
                               setup[ f"{subdomain_key}-dimensions_y" ],
                               1,
                               setup[ "domain_origin_x" ] + search_radius_L0 * setup[ f'{subdomain_key}-origin_glob_coords_x' ],
                               setup[ "domain_origin_y" ] + search_radius_L0 * setup[ f'{subdomain_key}-origin_glob_coords_y' ],
                               0,
                               np.zeros(( setup[ f"{subdomain_key}-dimensions_x" ] * setup[ f"{subdomain_key}-dimensions_y" ] )), # deprecated gridSector,
                               setup[ f"{subdomain_key}-search_radius" ],
                               f"sources/{subdomain_key}dambreak_grid.vtk" )

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

    dambreak_setup = {
        # geometri proportions
        "box_height" : args.box_height,
        "box_length" : args.box_length,
        "fluid_height" : args.fluid_height,
        "fluid_length" : args.fluid_length,
        # general parameteres
        "dp" : args.dp,
        "h_coef" : args.h_coef,
        "n_boundary_layers" : int( np.ceil( 2 * args.h_coef * args.dp  / args.dp ) ),
        "density" : args.density,
        "speed_of_sound" : args.speed_of_sound,
        "cfl" : args.cfl,
        "particle_mass" : args.density * ( args.dp * args.dp ),
        "boundary_element_size" : args.dp * args.dp,
        "smoothing_length" : args.h_coef * args.dp,
        "search_radius" : 2 * args.h_coef * args.dp,
        "time_step" : args.cfl * ( args.h_coef * args.dp ) / args.speed_of_sound,
        "alpha" : args.alpha,
        "dynamic_viscosity" : args.dynamic_viscosity,
        # terms and formulations
        "bc_type" : args.bc_type,
        "diffusive_term" : args.diffusive_term,
        "viscous_term" : args.viscous_term,
        ## refinement
        "number_of_subdomains" : 2,
        # master subdomain
        "subdomain-0-factor" : 1,
        "subdomain-0-x_min" : -np.inf,
        "subdomain-0-x_max" : 1.0,
        "subdomain-0-y_min" : -np.inf,
        "subdomain-0-y_max" : np.inf,
        # finer subdomain
        "subdomain-1-factor" : 0.5,
        "subdomain-1-x_min" : 1.0,
        "subdomain-1-x_max" : np.inf,
        "subdomain-1-y_min" : -np.inf,
        "subdomain-1-y_max" : np.inf,
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
    import init_generate_standard as single_resolution
    single_resolution.generate_dam_break_fluid_particles( dambreak_setup )
    single_resolution.generate_dam_break_boundary_particles( dambreak_setup )
    single_resolution.compute_domain_size( dambreak_setup )

    compute_subdomain_size( dambreak_setup )

    print( "Complete example setup:" )
    pprint( dambreak_setup )

    generate_dam_break_fluid_particles( dambreak_setup )
    generate_dam_break_boundary_particles( dambreak_setup )


    # write simulation params

    save_grid( dambreak_setup )

    write_simulation_params( dambreak_setup )
    write_distributed_domain_params( dambreak_setup )
    configure_and_write_measuretool_parameters( dambreak_setup )
