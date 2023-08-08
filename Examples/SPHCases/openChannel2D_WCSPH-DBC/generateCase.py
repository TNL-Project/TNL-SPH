#---------------------------------------------------------------------------#
#
# case: openChannelWithObstacle (without obstacle atm)
#
#---------------------------------------------------------------------------#
### Parameters of the case necessary for case creation:

## Parameters
dp = 0.002
smoothingLentghCoef = 2**0.5
#smoothingLentghCoef = 2

rho0 = 1000.
p0 = 0.

numberOfBoundaryLayers = 3

speedOfSound = 34.3
CFLnumber = 0.2
timeStep = 0.00002*0.5 #otherwise is obtained automatically

write = '.vtk' #.ptcs or .vtk

boxL = 1.
boxH = 0.3

fluidL = 0.6
fluidH = 0.1
numberOfAllocatedParticles = 50000

## First inlet buffer. ##
inletBufferOrientation_x = 1.
inletBufferOrientation_z = 0.
inletBufferPosition_x = 0.1
inletBufferPosition_z = 0. + dp*1
inletBufferHeight = 0.1
inletBufferLayers = numberOfBoundaryLayers + 1
inletVelocity_x = 1.
inletVelocity_z = 0.

inletBufferWidth = inletBufferLayers * dp # - dp / 2
inletBufferEdge = inletBufferPosition_x + 4 * dp  + dp / 2 #remove, deprecated
inletBufferReferencePoint_x = inletBufferPosition_x - inletBufferOrientation_x * ( inletBufferLayers - 1 ) * dp
inletBufferReferencePoint_z = inletBufferPosition_z - inletBufferOrientation_z * ( inletBufferLayers - 1 ) * dp

## Second inlet buffer. ##
inlet2BufferOrientation_x = -1.
inlet2BufferOrientation_z = 0.
inlet2BufferPosition_x = 0.7 + dp
inlet2BufferPosition_z = 0. + dp*1
inlet2BufferHeight = 0.1
inlet2BufferLayers = numberOfBoundaryLayers + 1
inlet2Velocity_x = 1.5
inlet2Velocity_z = 0.

inlet2BufferWidth = inlet2BufferLayers * dp # - dp / 2
inlet2BufferEdge = inlet2BufferPosition_x  + dp / 2 #remove, deprecated
inlet2BufferReferencePoint_x = inlet2BufferPosition_x - inlet2BufferOrientation_x * ( inlet2BufferLayers - 1 ) * dp
inlet2BufferReferencePoint_z = inlet2BufferPosition_z - inlet2BufferOrientation_z * ( inlet2BufferLayers - 1 ) * dp

#---------------------------------------------------------------------------#

### Create related directories
import os
resultsPath = r'./results'
if not os.path.exists( resultsPath ):
    os.makedirs( resultsPath )

sourcesPath = r'./sources'
if not os.path.exists( sourcesPath ):
    os.makedirs( sourcesPath )

boxL_n = round( boxL / dp )
boxH_n = round( boxH / dp )

fluidL_n = round( fluidL / dp )
fluidH_n = round( fluidH / dp )

inletL_n = inletBufferLayers
inletH_n = round( inletBufferHeight / dp  )

inlet2L_n = inlet2BufferLayers
inlet2H_n = round( inlet2BufferHeight / dp  )

### Generate fluid particles
fluid_rx = []; fluid_ry = []; fluid_rz = []
fluid_density = []

for x in range( fluidL_n ):
    for z in range( fluidH_n ):
        fluid_rx.append( inletBufferPosition_x + dp * ( x + 1 ) )
        fluid_ry.append( 0. ) #we use only 2D case
        fluid_rz.append( dp * ( z + 1 ) )

        hydrostaticPressure = rho0 * 9.81 * ( fluidH - z * dp )
        hydrostaticDensity = ( ( hydrostaticPressure / ( speedOfSound ** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0;
        fluid_density.append( hydrostaticDensity )

### Generate buffer particles
inlet_rx = []; inlet_ry = []; inlet_rz = []
inlet_vx = []; inlet_vy = []; inlet_vz = []
inlet_density = []

for x in range( inletL_n ):
    for z in range( inletH_n ):
        inlet_rx.append( inletBufferPosition_x - inletBufferOrientation_x * dp * ( x ) )
        inlet_ry.append( 0. ) #we use only 2D case
        inlet_rz.append( inletBufferPosition_z + dp * ( z ) )

        inlet_vx.append( inletVelocity_x )
        inlet_vy.append( 0. ) #we use only 2D case
        inlet_vz.append( inletVelocity_z )

        hydrostaticPressure = rho0 * 9.81 * ( fluidH - z * dp )
        hydrostaticDensity = ( ( hydrostaticPressure / ( speedOfSound ** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0;
        inlet_density.append( hydrostaticDensity )

### Generate buffer particles
inlet2_rx = []; inlet2_ry = []; inlet2_rz = []
inlet2_vx = []; inlet2_vy = []; inlet2_vz = []
inlet2_density = []

for x in range( inlet2L_n ):
    for z in range( inlet2H_n ):
        inlet2_rx.append( inlet2BufferPosition_x - inlet2BufferOrientation_x * dp * ( x ) )
        inlet2_ry.append( 0. ) #we use only 2D case
        inlet2_rz.append( inlet2BufferPosition_z + dp * ( z ) )

        inlet2_vx.append( inlet2Velocity_x )
        inlet2_vy.append( 0. ) #we use only 2D case
        inlet2_vz.append( inlet2Velocity_z )

        hydrostaticPressure = rho0 * 9.81 * ( fluidH - z * dp )
        hydrostaticDensity = ( ( hydrostaticPressure / ( speedOfSound ** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0;
        inlet2_density.append( hydrostaticDensity )

### Generate boundary particles
box_rx = []; box_ry = []; box_rz = []
ghost_rx = []; ghost_ry = []; ghost_rz = []
box_density = [];

#:# left wall
#:for layer in range( numberOfBoundaryLayers ):
#:    for z in range( boxH_n - 1 ):
#:        box_rx.append( 0. - layer * dp )
#:        box_ry.append( 0. ) #we use only 2D case
#:        box_rz.append( ( z+1 ) * dp )
#:        box_density.append( rho0 );

# bottom wall
for layer in range( numberOfBoundaryLayers ):
    for x in range( boxL_n + ( numberOfBoundaryLayers - 1 ) * 2 ):
        box_rx.append( ( x - ( numberOfBoundaryLayers - 1 ) ) * dp )
        box_ry.append( 0. ) #we use only 2D case
        box_rz.append( 0. - layer * dp )

        ghost_rx.append( ( x - ( numberOfBoundaryLayers - 1 ) ) * dp )
        ghost_ry.append( 0. + dp * ( layer + 1 ) )
        ghost_rz.append( 0.)

        #hydrostaticPressure = rho0 * 9.81 * ( fluidH - z * dp )
        hydrostaticPressure = rho0 * 9.81 * ( fluidH + layer * dp )
        hydrostaticDensity = ( ( hydrostaticPressure / ( speedOfSound ** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0;
        box_density.append( hydrostaticDensity )

#:x_last = box_rx[-1 -(numberOfBoundaryLayers - 1)] #due to discretisation, we need to save last value of bottom wall
#:
#:# right wall
#:for layer in range( numberOfBoundaryLayers ):
#:    for z in range( boxH_n - 1 ):
#:        box_rx.append( x_last + dp * layer )
#:        box_ry.append( 0. ) #we use only 2D case
#:        box_rz.append( ( z + 1 ) * dp )
#:        box_density.append( rho0 );

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
import sys
sys.path.append('../../tools/')
import saveParticlesVTK
import numpy as np

r = np.array( ( fluid_rx, fluid_rz, fluid_ry ), dtype=float ).T #!!
v = np.zeros( ( len( fluid_rx ), 3 ) )
#rho = rho0 * np.ones( len( fluid_rx ) )
rho = np.array( fluid_density, dtype=float )
p = np.zeros( len( fluid_rx ) )
ptype = np.zeros( len( fluid_rx ) )

fluidToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
saveParticlesVTK.save_polydata( fluidToWrite, "sources/openchannel_fluid.vtk" )

r = np.array( ( box_rx, box_rz, box_ry ), dtype=float ).T #!!
v = np.zeros( ( len( box_rx ), 3 ) )
gn = np.array( ( ghost_rx, ghost_ry, ghost_rz ), dtype=float ).T #!!
#rho = rho0 * np.ones( len( box_rx ) )
rho = np.array( box_density, dtype=float )
p = np.zeros( len( box_rx ) )
ptype = np.ones( len( box_rx ) )

boxToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype, ghostNodes=gn )
saveParticlesVTK.save_polydata( boxToWrite, "sources/openchannel_boundary.vtk" )

r = np.array( ( inlet_rx, inlet_rz, inlet_ry ), dtype=float ).T #!!
v = np.array( ( inlet_vx, inlet_vz, inlet_vy ), dtype=float ).T #!!
#rho = rho0 * np.ones( len( inlet_rx ) )
rho = np.array( inlet_density, dtype=float )
p = np.zeros( len( inlet_rx ) )
ptype = np.ones( len( inlet_rx ) )

inletToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
saveParticlesVTK.save_polydata( inletToWrite, "sources/openchannel_inlet.vtk" )

r = np.array( ( inlet2_rx, inlet2_rz, inlet2_ry ), dtype=float ).T #!!
v = np.array( ( inlet2_vx, inlet2_vz, inlet2_vy ), dtype=float ).T #!!
#rho = rho0 * np.ones( len( inlet2_rx ) )
rho = np.array( inlet2_density, dtype=float )
p = np.zeros( len( inlet2_rx ) )
ptype = np.ones( len( inlet2_rx ) )

inlet2ToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
saveParticlesVTK.save_polydata( inlet2ToWrite, "sources/openchannel_outlet.vtk" )

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

#Determine parameters of inlet1

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
with open( 'sources/SPHCaseConfig.h', 'w' ) as file:
  file.write( fileSPHConf )

with open( 'template/ParticlesConfig_template.h', 'r' ) as file :
  fileParticleConf = file.read()

# Replace the target string
fileParticleConf = fileParticleConf.replace( 'placeholderDimension', str( spaceDimension ) )
#fileParticleConf = fileParticleConf.replace( 'placeholderFluidParticles', str( 0 ) )
fileParticleConf = fileParticleConf.replace( 'placeholderFluidParticles', str( len( fluid_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedFluidParticles', str( numberOfAllocatedParticles ) )
fileParticleConf = fileParticleConf.replace( 'placeholderBoundaryParticles', str( len( box_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedBoundaryParticles', str( len( box_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderInletParticles', str( len( inlet_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedInletParticles', str( len( inlet_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderOutletParticles', str( len( inlet2_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedOutletParticles', str( len( inlet2_rx ) * 3 ) )

fileParticleConf = fileParticleConf.replace( 'placeholderSearchRadius', str( searchRadius ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXSize', str( gridXsize ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYSize', str( gridYsize ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXBegin', str( round( gridXbegin, 9  ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYBegin', str( round( gridYbegin, 9  ) ) )

# Write the file out again
with open( 'sources/ParticlesConfig.h', 'w' ) as file:
  file.write( fileParticleConf )

with open( 'template/OpenBoundaryConfig_template.h', 'r' ) as file :
  fileOBConf = file.read()

#inlet1
fileOBConf = fileOBConf.replace( 'placeholderInletOrientation_x', str( inletBufferOrientation_x ) )
fileOBConf = fileOBConf.replace( 'placeholderInletOrientation_y', str( inletBufferOrientation_z ) )
fileOBConf = fileOBConf.replace( 'placeholderInletVelocity_x', str( inletVelocity_x ) )
fileOBConf = fileOBConf.replace( 'placeholderInletVelocity_y', str( inletVelocity_z ) )
fileOBConf = fileOBConf.replace( 'placeholderInletPosition_x', str( inletBufferPosition_x  + dp/2 ) ) #FIXME
fileOBConf = fileOBConf.replace( 'placeholderInletPosition_y', str( inletBufferPosition_z ) )
fileOBConf = fileOBConf.replace( 'placeholderInletDensity', str( rho0 ) )
fileOBConf = fileOBConf.replace( 'placeholderInletWidth_x', str( round( inletBufferWidth, 7 ) ) )
fileOBConf = fileOBConf.replace( 'placeholderInletWidth_y', str( 0. ) )
fileOBConf = fileOBConf.replace( 'placeholderInletBufferEdge', str( round(  inletBufferEdge, 7 ) ) )

#outlet
fileOBConf = fileOBConf.replace( 'placeholderOutletOrientation_x', str( inlet2BufferOrientation_x ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletOrientation_y', str( inlet2BufferOrientation_z ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletVelocity_x', str( inlet2Velocity_x ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletVelocity_y', str( inlet2Velocity_z ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletPosition_x', str( inlet2BufferPosition_x - dp/2 ) ) #FIXME
fileOBConf = fileOBConf.replace( 'placeholderOutletPosition_y', str( inlet2BufferPosition_z ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletDensity', str( rho0 ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletWidth_x', str( round( inlet2BufferWidth, 7 ) ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletWidth_y', str( 0. ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletBufferEdge', str( round(  inlet2BufferEdge, 7 ) ) )

# Write the file out again
with open( 'sources/OpenBoundaryConfig.h', 'w' ) as file:
  file.write( fileOBConf )

# Read and write (with possible edit) simulation control file.
with open( 'template/SimulationControlConfig.h', 'r' ) as file :
  fileSimulationControl = file.read()

with open( 'sources/SimulationControlConfig.h', 'w' ) as file:
  file.write( fileSimulationControl )

# Read and write (with possible edit) measuretool config file.
with open( 'template/MeasuretoolConfig.h', 'r' ) as file :
  fileMeasuretoolConf = file.read()

with open( 'sources/MeasuretoolConfig.h', 'w' ) as file:
  file.write( fileMeasuretoolConf )
