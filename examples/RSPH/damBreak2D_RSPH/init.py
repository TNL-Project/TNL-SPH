#! /usr/bin/env python3

import numpy as np
import sys
sys.path.append('../../../src/tools')
import saveParticlesVTK

def generate_dam_break_fluid_particles( dp, fluid_lenght, fluid_height, density ):
    fluid_rx = []; fluid_ry = []
    fluid_lenght_n = round( fluid_lenght / dp )
    fluid_height_n = round( fluid_height / dp )

    for x in range( fluid_lenght_n ):
        for y in range( fluid_height_n ):
            fluid_rx.append( dp * ( x + 1 ) )
            fluid_ry.append( dp * ( y + 1 ) )

    fluid_n = len( fluid_rx )
    fluid_r = np.array( ( fluid_rx, fluid_ry, np.zeros( fluid_n ) ), dtype=float ).T #!!
    fluid_v = np.zeros( ( fluid_n, 3 ) )
    fluid_rho = density * np.ones( fluid_n )
    fluid_p = np.zeros( fluid_n )
    fluid_ptype = np.zeros( fluid_n )
    fluid_to_write = saveParticlesVTK.create_pointcloud_polydata( fluid_r, fluid_v, fluid_rho, fluid_p, fluid_ptype )
    saveParticlesVTK.save_polydata( fluid_to_write, "sources/dambreak_fluid.vtk" )

    return fluid_n

def generate_dam_break_boundary_particles( dp, box_lenght, box_height, n_boundary_layers, density ):
    box_rx = []; box_ry = []
    normal_x = []; normal_y = []
    box_length_n = round( box_lenght / dp )
    box_height_n = round( box_height / dp )

    # left wall
    for layer in range( n_boundary_layers ):
        for y in range( box_height_n - 1 ):
            box_rx.append( 0. - layer * dp )
            box_ry.append( ( y + 1 ) * dp )
            normal_x.append( 1. )
            normal_y.append( 0. )

    # bottom wall
    for layer in range( n_boundary_layers ):
        for x in range( box_length_n - 2 ):
            box_rx.append( ( x + 1 ) * dp )
            box_ry.append( 0. - layer * dp )
            normal_x.append( 0. )
            normal_y.append( 1. )

    x_last = box_rx[ -1 ] + dp #due to discretisation, we need to save last value of bottom wall

    # right wall
    for layer in range( n_boundary_layers ):
        for y in range( box_height_n - 1 ):
            box_rx.append( x_last + dp * layer )
            box_ry.append( ( y + 1 ) * dp )
            normal_x.append( -1. )
            normal_y.append( 0. )

    # generate the corners
    def generate90degCorner( x, y, dirx, diry ):
      for layer in range( n_boundary_layers ):
        for k in range( n_boundary_layers ):
          box_rx.append( x + k * dp * dirx )
          box_ry.append( y + layer * dp * diry )
          nx = ( x + ( -1 ) * dirx * dp ) - ( x + k * dp * dirx )
          ny = ( y + ( -1 ) * diry * dp ) - ( y + layer * dp * diry )
          n_norm = ( nx**2 + ny**2 )**0.5
          normal_x.append( nx / n_norm )
          normal_y.append( ny / n_norm )

    generate90degCorner( 0, 0., -1, -1 )
    generate90degCorner( x_last, 0., +1, -1 )

    boundary_n = len( box_rx )
    boundary_r = np.array( ( box_rx, box_ry, np.zeros( boundary_n ) ), dtype=float ).T #!!
    boundary_normal = np.array( ( normal_x, normal_y, np.zeros( boundary_n ) ), dtype=float ).T #!!
    boundary_v = np.zeros( ( boundary_n, 3 ) )
    boundary_rho = density * np.ones( boundary_n )
    boundary_p = np.zeros( boundary_n )
    boundary_ptype = np.ones( boundary_n )
    box_to_write = saveParticlesVTK.create_pointcloud_polydata( boundary_r, boundary_v, boundary_rho, boundary_p, boundary_ptype,
                                                                boundary_normal )
    saveParticlesVTK.save_polydata( box_to_write, "sources/dambreak_boundary.vtk" )

    domain_origin_x = min( box_rx )
    domain_origin_y = min( box_ry )
    domain_end_x = max( box_rx )
    domain_end_y = max( box_ry )
    return boundary_n, [ domain_origin_x, domain_end_x ], [ domain_origin_y, domain_end_y ]


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

    # write config with templated configuration
    with open( 'template/config_template.h', 'r' ) as file :
      config_file = file.read()
    with open( 'sources/config.h', 'w' ) as file:
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
    g = argparser.add_argument_group("resolution parameters")
    g.add_argument("--dp", type=float, default=0.002, help="initial distance between particles")
    g.add_argument("--h-coef", type=float, default=2**0.5, help="smoothing length coefitient")
    g = argparser.add_argument_group("domain parameters")
    g.add_argument("--box-length", type=float, default=1.61, help="length of dam break box")
    g.add_argument("--box-height", type=float, default=1.0, help="height of dam break box")
    g.add_argument("--fluid-length", type=float, default=0.6, help="length of fluid block")
    g.add_argument("--fluid-height", type=float, default=0.3, help="height of fluid block")
    g.add_argument("--n-boundary-layers", type=int, default=3, help="number of boundary layers")
    g = argparser.add_argument_group("simulation parameters")
    g.add_argument("--density", type=float, default=1000, help="referential density of the fluid")
    g.add_argument("--speed-of-sound", type=float, default=34.3, help="speed of sound")
    g.add_argument("--cfl", type=float, default=0.08, help="referential density of the fluid")
    #g = argparser.add_argument_group("control parameters")
    #g.add_argument("--example-dir", type=Path, default=1000, help="referential density of the fluid")

    args = argparser.parse_args()

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
