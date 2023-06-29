#---------------------------------------------------------------------------#
#
# case: periodicOpenChannel2D_WCSPH-DBC
# case destription: 2D open channel, consisting of bottom wall and periodic
#                   boundary conditions driven by gravitational force which
#                   models the slope of the channel
#
#---------------------------------------------------------------------------#

### Parameters of the case necessary for case creation:

# Dimensions of the fluid wall, resp channel [m]:
channelWidth = 0.3
waterLevel = 0.3

# Initial particle distance (dp)[m]:
dp = 0.002

# Smoothing length coefitient:
# - smoothing length (h)[m] = smoothing length coef (Coef_h)[-] * initial particle distance (d_p)[m]
# ( h = Coef_h * dp )
smoothingLentghCoef = 2**0.5

# Referential density of the medium (rho0)[kg/m^3]:
rho0 = 1000.

# Number of boundary layers (n_layer)[-]:
numberOfBoundaryLayers = 3

# Numerical speed of sound (c0)[m/s]:
speedOfSound = 34.3

# Initial time step (dtInit)[s].
# - in case, that initial time step is not defined, is computed automatically.
timeStep = 0.00002

# CFL number (CFL)[-]:
CFLnumber = 0.2

write = '.vtk' #.ptcs or .vtk

#---------------------------------------------------------------------------#

wallL_n = round( channelWidth / dp )
channelWidth_n = round( channelWidth / dp )
waterLevel_n = round( waterLevel / dp )

### Generate fluid particles
fluid_rx = []; fluid_ry = []; fluid_rz = []

rhovect = [] #debug
for x in range( channelWidth_n ):
    for z in range( waterLevel_n ):
        fluid_rx.append( dp / 2 + dp * x )
        fluid_ry.append( 0. ) #we use only 2D case
        fluid_rz.append( dp * ( z + 1 ) )
        rhovect.append( x )

### Generate boundary particles
wall_rx = []; wall_ry = []; wall_rz = []

# bottom wall
for layer in range( numberOfBoundaryLayers ):
    for x in range( wallL_n ):
        wall_rx.append( dp / 2  + x  * dp )
        wall_ry.append( 0. ) #we use only 2D case
        wall_rz.append( 0. - layer * dp )

#---------------------------------------------------------------------------#

### Write fluid particles
if write == '.ptcs':
    with open( "periodicOpenChannel_fluid.ptcs", "w" ) as f:
        f.write( str( len( fluid_rx ) ) + "\n" )
        for i in range( len( fluid_rx ) ):
            f.write( str( round( fluid_rx[ i ], 5 ) ) + " " + str( round( fluid_rz[i], 5 ) ) + " " + \
                     str( round( fluid_ry[ i ], 5 ) ) + " " + str( 0. ) + " " + str( 0. ) + " " + str( 0. ) + " " + \
                     str( round( rho0, 5 ) ) + " " + str( round( p0, 5 ) ) + " " + str( 0 ) + "\n" )

    ### Write fluid particles
    with open("periodicOpenChannel_boundary.ptcs", "w") as f:
        f.write( str( len( wall_rx ) ) + "\n" )
        for i in range( len( wall_rx ) ):
            f.write( str( round( wall_rx[ i ], 5 ) ) + " " + str( round( wall_rz[ i ], 5 ) ) + " " + \
                     str( round( wall_ry[ i ], 5 ) ) + " " + str( 0. ) + " " + str( 0. ) + " " + str( 0. ) + " " + \
                     str( round( rho0, 5 ) ) + " " + str( round( p0, 5 ) ) + " " + str( 1 ) + "\n" )
elif write == '.vtk':
    import sys
    sys.path.append('../../tools/')
    import saveParticlesVTK
    import numpy as np

    r = np.array( ( fluid_rx, fluid_rz, fluid_ry ), dtype=float ).T #!!
    v = np.zeros( ( len( fluid_rx ), 3 ) )
    rho = rho0 * np.ones( len( fluid_rx ) )
    #rho = np.array( rhovect, dtype=float ).T
    p = np.zeros( len( fluid_rx ) )
    ptype = np.zeros( len( fluid_rx ) )

    fluidToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
    saveParticlesVTK.save_polydata( fluidToWrite, "periodicOpenChannel_fluid.vtk" )

    r = np.array( ( wall_rx, wall_rz, wall_ry ), dtype=float ).T #!!
    v = np.zeros( ( len( wall_rx ), 3 ) )
    rho = rho0 * np.ones( len( wall_rx ) )
    p = np.zeros( len( wall_rx ) )
    ptype = np.ones( len( wall_rx ) )

    wallToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
    saveParticlesVTK.save_polydata( wallToWrite, "periodicOpenChannel_boundary.vtk" )
else:
    print( "Invalid particle output type." )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
### Compute remaining parameters
particleMass = rho0 * ( dp * dp )
smoothingLentgh =  round( smoothingLentghCoef * dp, 7 )
searchRadius = round( smoothingLentgh * 2 , 7 )
if not timeStep:
    timeStep = round( CFLnumber * ( smoothingLentgh / speedOfSound ), 8 )
coefB = round( speedOfSound * speedOfSound * rho0 / 7 , 1 )
spaceDimension = 2

#Determine grid size
import math
gridXbegin = ( 0 - searchRadius )
gridYbegin = 1.01 * ( ( min( min( fluid_rz ), min( wall_rz ) ) ) - searchRadius )

gridXend = ( ( channelWidth - dp / 2 ) + searchRadius )
gridYend = 1.01 * ( ( max( max( fluid_rz ), max( wall_rz ) ) ) + searchRadius )

gridXsize = math.ceil( ( gridXend - gridXbegin ) / searchRadius )
gridYsize = math.ceil( ( gridYend - gridYbegin ) / searchRadius )

# Read in the file
with open( 'template/SPHCaseConfig_template.h', 'r' ) as file :
  fileSPHConf = file.read()

# Replace the target string
fileSPHConf = fileSPHConf.replace( 'placeholderDimension', str( spaceDimension ) )
fileSPHConf = fileSPHConf.replace( 'placeholderMass', str( particleMass ) )
fileSPHConf = fileSPHConf.replace( 'placeholderSpeedOfSound', str( speedOfSound ) )
fileSPHConf = fileSPHConf.replace( 'placeholderCoefB', str( coefB ) )
fileSPHConf = fileSPHConf.replace( 'placeholderDensity', str( rho0 ))
fileSPHConf = fileSPHConf.replace( 'placeholderInitParticleDistance', str( dp ) )
fileSPHConf = fileSPHConf.replace( 'placeholderSmoothingLength', str( smoothingLentgh ) )
fileSPHConf = fileSPHConf.replace( 'placeholderTimeStep', str( timeStep ) )

# Write the file out again
with open( 'SPHCaseConfig.h', 'w' ) as file:
  file.write( fileSPHConf )

with open( 'template/ParticlesConfig_template.h', 'r' ) as file :
  fileParticleConf = file.read()

# Replace the target string
fileParticleConf = fileParticleConf.replace( 'placeholderDimension', str( spaceDimension ) )
fileParticleConf = fileParticleConf.replace( 'placeholderFluidParticles', str( len( fluid_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedFluidParticles', str( ( int )( len( fluid_rx ) * 1.2 ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderBoundaryParticles', str( len( wall_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedBoundaryParticles', str( ( int )( len( wall_rx ) * 1.2 ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderSearchRadius', str( searchRadius ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXSize', str( gridXsize ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYSize', str( gridYsize ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXBegin', str( round( gridXbegin, 9  ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYBegin', str( round( gridYbegin, 9  ) ) )
# - periodic boundary conditions:
fileParticleConf = fileParticleConf.replace( 'placeholderPeriodicityLength', str( round( channelWidth, 9  ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderFirstPeriodicColumn', str( 1 ) )
fileParticleConf = fileParticleConf.replace( 'placeholderLastPeriodicColumn', str( gridXsize - 1 ) )

# Write the file out again
with open( 'ParticlesConfig.h', 'w' ) as file:
  file.write( fileParticleConf )

# Read and write (with possible edit) simulation control file.
with open( 'template/SimulationControlConfig.h', 'r' ) as file :
  fileSimulationControl = file.read()

with open( 'SimulationControlConfig.h', 'w' ) as file:
  file.write( fileSimulationControl )

# Read and write (with possible edit) measuretool config file.
with open( 'template/MeasuretoolConfig.h', 'r' ) as file :
  fileMeasuretoolConf = file.read()

# Replace the target string
fileMeasuretoolConf = fileMeasuretoolConf.replace( 'placeholderInitParticleDistance', str( dp ) )
fileMeasuretoolConf = fileMeasuretoolConf.replace( 'placeholderSmoothingLength', str( smoothingLentgh ) )

with open( 'MeasuretoolConfig.h', 'w' ) as file:
  file.write( fileMeasuretoolConf )

#---------------------------------------------------------------------------#
### Parameters of the case necessary for case creation:

searchRadius_h = round( smoothingLentgh * 2, 7 )
gridSector = []
for y in range ( gridYsize ):
    for x in range ( gridXsize ):
        gridSector.append( 0 )

from contextlib import redirect_stdout
def DomainGrid( gridXsize, gridYsize, gridZsize, gridXbegin, gridYbegin, gridZbegin, gridSector, name ):
    with open( name, 'w' ) as f:
        with redirect_stdout(f):
            print( "# vtk DataFile Version 3.0" )
            print( "vtk output" )
            print( "ASCII" )
            print( "DATASET STRUCTURED_POINTS" )
            print( "DIMENSIONS ", gridXsize + 1 , " ", gridYsize + 1, " ", 1 )
            print( "ASPECT_RATIO ", searchRadius_h , " ", searchRadius_h , " ",  searchRadius_h )
            print( "ORIGIN ", gridXbegin , " ", gridYbegin , " ",  0  )
            print( "CELL_DATA ",  gridXsize * gridYsize * 1  )
            print( "SCALARS GridSector int 1 ")
            print( "LOOKUP_TABLE default" )
            for i in range( gridXsize * gridYsize * 1 ):
                print( gridSector[ i ] )

# Write global grid as example
DomainGrid( gridXsize, gridYsize, 1,    # grid size
            gridXbegin, gridYbegin, 0,  # coordinates of rgrid origina
            gridSector,                 # array with index of grid sector
            "periodicOpenChannel_grid.vtk" )     # outputfile name
#---------------------------------------------------------------------------#

import os
resultsPath = r'./results'
if not os.path.exists( resultsPath ):
    os.makedirs( resultsPath )
