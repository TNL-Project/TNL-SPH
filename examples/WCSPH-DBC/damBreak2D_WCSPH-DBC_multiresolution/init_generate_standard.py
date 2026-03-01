import numpy as np
import sys
sys.path.append('../../../src/tools')
import saveParticlesVTK

def compute_hydrostatic_density( ry, fluid_height, rho0, speed_of_sound ):
    hydrostaticPressure = rho0 * 9.81 * ( fluid_height - ry )
    hydrostaticDensity = ( ( hydrostaticPressure / ( speed_of_sound** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0
    #hydrostaticDensity = rho0 + hydrostaticPressure / speed_of_sound**2
    return hydrostaticDensity


def generate_dam_break_fluid_particles( setup ):
    fluid_rx = []; fluid_ry = []
    fluid_density = []
    dp = setup [ 'dp' ]
    fluid_length = setup[ 'fluid_length' ];
    fluid_height = setup[ 'fluid_height' ]
    fluid_lenght_n = round( fluid_length / dp )
    fluid_height_n = round( fluid_height / dp )
    rho0 = setup[ 'density' ]
    speed_of_sound = setup[ 'speed_of_sound' ]

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
    setup[ "fluid_n" ] = fluid_n

def generate_dam_break_boundary_particles( setup ):
    box_rx = []; box_ry = []
    ghost_rx = []; ghost_ry = []
    normal_x = []; normal_y = []
    dp = setup[ 'dp' ]
    n_boundary_layers = setup[ 'n_boundary_layers' ]
    box_length_n = round( setup[ 'box_length' ] / dp )
    box_height_n = round( setup[ 'box_height' ] / dp )
    rho0 = setup[ 'density' ]

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
