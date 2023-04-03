# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#
# damBreak2D_WCSPH-DBC_benchmark
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

boxL = 1.61
boxH = 0.8

fluidL = 0.6
fluidH = 0.3

#dp = 0.002
dp = 0.005
smoothingLentghCoef = 2**0.5*1.3

rho0 = 1000.
p0 = 0.

numberOfBoundaryLayers = 1 #Since we use BI

speedOfSound = 34.3
CFLnumber = 0.2
timeStep = 0.00002 #otherwise is obtained automatically

write = '.vtk' #.ptcs or .vtk

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

boxL_n = round( boxL / dp )
boxH_n = round( boxH / dp )

fluidL_n = round( fluidL / dp )
fluidH_n = round( fluidH / dp )

### Generate fluid particles
fluid_rx = []; fluid_ry = []; fluid_rz = []

for x in range( fluidL_n ):
    for z in range( fluidH_n ):
        fluid_rx.append( dp * ( x + 1 ) )
        fluid_ry.append( 0. ) #we use only 2D case
        fluid_rz.append( dp * ( z + 1 ) )

### Generate boundary particles
box_rx = []; box_ry = []; box_rz = []
#box normals
box_nx = []; box_ny = []; box_nz = []


# left wall
for layer in range( numberOfBoundaryLayers ):
    for z in range( boxH_n - 1 ):
        box_rx.append( 0. - layer * dp )
        box_ry.append( 0. ) #we use only 2D case
        box_rz.append( ( z+1 ) * dp )
        box_nx.append( 1. )
        box_ny.append( 0. )
        box_nz.append( 0. )

# bottom wall
for layer in range( numberOfBoundaryLayers ):
    for x in range( boxL_n + ( numberOfBoundaryLayers - 1 ) * 2 - 2 ):
        box_rx.append( ( x - ( numberOfBoundaryLayers - 1 ) + 1 ) * dp )
        box_ry.append( 0. ) #we use only 2D case
        box_rz.append( 0. - layer * dp )
        box_nx.append( 0. )
        box_ny.append( 0. )
        box_nz.append( 1. )

x_last = box_rx[-1 -(numberOfBoundaryLayers - 1)] #due to discretisation, we need to save last value of bottom wall

# corners
box_rx.append( 0. )
box_ry.append( 0. ) #we use only 2D case
box_rz.append( 0. )
box_nx.append( 0.5**0.5 )
box_ny.append( 0.)
box_nz.append( 0.5**0.5 )

box_rx.append( x_last + dp )
box_ry.append( 0. ) #we use only 2D case
box_rz.append( 0. )
box_nx.append( -0.5**0.5 )
box_ny.append( 0. )
box_nz.append( 0.5**0.5 )

# right wall
for layer in range( numberOfBoundaryLayers ):
    for z in range( boxH_n - 1 ):
        box_rx.append( x_last + dp * layer )
        box_ry.append( 0. ) #we use only 2D case
        box_rz.append( ( z + 1 ) * dp )
        box_nx.append( -1. )
        box_ny.append( 0. )
        box_nz.append( 0. )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

### Write fluid particles
if write == '.ptcs':
    with open( "dambreak_fluid.ptcs", "w" ) as f:
        f.write( str( len( fluid_rx ) ) + "\n" )
        for i in range( len( fluid_rx ) ):
            f.write( str( round( fluid_rx[ i ], 5 ) ) + " " + str( round( fluid_rz[i], 5 ) ) + " " + \
                     str( round( fluid_ry[ i ], 5 ) ) + " " + str( 0. ) + " " + str( 0. ) + " " + str( 0. ) + " " + \
                     str( round( rho0, 5 ) ) + " " + str( round( p0, 5 ) ) + " " + str( 0 ) + "\n" )

    ### Write fluid particles
    with open("dambreak_boundary.ptcs", "w") as f:
        f.write( str( len( box_rx ) ) + "\n" )
        for i in range( len( box_rx ) ):
            f.write( str( round( box_rx[ i ], 5 ) ) + " " + str( round( box_rz[ i ], 5 ) ) + " " + \
                     str( round( box_ry[ i ], 5 ) ) + " " + str( 0. ) + " " + str( 0. ) + " " + str( 0. ) + " " + \
                     str( round( rho0, 5 ) ) + " " + str( round( p0, 5 ) ) + " " + str( 1 ) + "\n" )
elif write == '.vtk':
    import sys
    sys.path.append('../../tools/')
    import saveParticlesVTK
    import numpy as np

    r = np.array( ( fluid_rx, fluid_rz, fluid_ry ), dtype=float ).T #!!
    v = np.zeros( ( len( fluid_rx ), 3 ) )
    rho = rho0 * np.ones( len( fluid_rx ) )
    p = np.zeros( len( fluid_rx ) )
    ptype = np.zeros( len( fluid_rx ) )

    fluidToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
    saveParticlesVTK.save_polydata( fluidToWrite, "dambreak_fluid.vtk" )

    r = np.array( ( box_rx, box_rz, box_ry ), dtype=float ).T #!!
    n = np.array( ( box_nx, box_nz, box_ny ), dtype=float ).T #!!
    v = np.zeros( ( len( box_rx ), 3 ) )
    rho = rho0 * np.ones( len( box_rx ) )
    p = np.zeros( len( box_rx ) )
    ptype = np.ones( len( box_rx ) )

    boxToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype, n )
    saveParticlesVTK.save_polydata( boxToWrite, "dambreak_boundary.vtk" )
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
gridXbegin = 1.005 * ( min( min( fluid_rx ), min( box_rx ) ) ) - searchRadius
gridYbegin = 1.005 * ( min( min( fluid_rz ), min( box_rz ) ) ) - searchRadius

gridXend = 1.005 * ( max( max( fluid_rx ), max( box_rx ) ) ) + searchRadius
gridYend = 1.005 * ( max( max( fluid_rz ), max( box_rz ) ) ) + searchRadius

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
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedFluidParticles', str( len( fluid_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderBoundaryParticles', str( len( box_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedBoundaryParticles', str( len( box_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderSearchRadius', str( searchRadius ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXSize', str( gridXsize ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYSize', str( gridYsize ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXBegin', str( round( gridXbegin, 9  ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYBegin', str( round( gridYbegin, 9  ) ) )

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

with open( 'MeasuretoolConfig.h', 'w' ) as file:
  file.write( fileMeasuretoolConf )
