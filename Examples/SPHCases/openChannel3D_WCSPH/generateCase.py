#---------------------------------------------------------------------------#
#
# case: openChannelWithObstacle (without obstacle atm)
#
#---------------------------------------------------------------------------#
### Parameters of the case necessary for case creation:

smoothingLentghCoef = 2
rho0 = 1000.
speedOfSound = 60.0
CFLnumber = 0.15

# Dimensions of the channel
fluidL = 0.5
fluidH = 0.1

# Obstacle position and size
obstacleCenterX = 0.2
obstacleCenterY = 0.07
obstacleR = 0.01

# Initial particle distance (dp)[m]:
dp = 0.005

#: DtIni=0.0001666666622179967
#: DtMin=8.333333235076177e-06

numberOfAllocatedParticles = 4e6

waterLevel = 0.2
#---------------------------------------------------------------------------#
rho0 = 1000.
p0 = 0.

numberOfBoundaryLayers = 3

#speedOfSound = 34.3
CFLnumber = 0.2
timeStep = 0.00002*0.5 #otherwise is obtained automatically

write = '.vtk' #.ptcs or .vtk

boxL = 1.
boxH = 0.3


## First inlet buffer. ##
inletBufferWidth_start = -0.495
inletBufferWidthY = 0.495 * 2
inletBufferOrientation_x = 1.
inletBufferOrientation_y = 0.
inletBufferOrientation_z = 0.
inletBufferPosition_x = -0.2 - dp
inletBufferPosition_y = 0. + dp
inletBufferPosition_z = 0. + dp*1
inletBufferHeight = waterLevel
inletBufferLayers = numberOfBoundaryLayers + 1
#inletVelocity_x = 1.
inletVelocity_x = 0.5
inletVelocity_y = 0.
inletVelocity_z = 0.

inletBufferWidth = inletBufferLayers * dp # - dp / 2
inletBufferEdge = inletBufferPosition_x + 4 * dp  + dp / 2 #remove, deprecated
inletBufferReferencePoint_x = inletBufferPosition_x - inletBufferOrientation_x * ( inletBufferLayers - 1 ) * dp
inletBufferReferencePoint_y = inletBufferPosition_y - inletBufferOrientation_y * ( inletBufferLayers - 1 ) * dp
inletBufferReferencePoint_z = inletBufferPosition_z - inletBufferOrientation_z * ( inletBufferLayers - 1 ) * dp

## Second inlet buffer. ##
inlet2BufferWidth_start = -0.495
inlet2BufferWidthY = 0.495 * 2
inlet2BufferOrientation_x = -1.
inlet2BufferOrientation_y = 0.
inlet2BufferOrientation_z = 0.
inlet2BufferPosition_x = 2. + dp
inlet2BufferPosition_y = 0. + dp
inlet2BufferPosition_z = 0. + dp*1
inlet2BufferHeight = waterLevel
inlet2BufferLayers = numberOfBoundaryLayers + 1
inlet2Velocity_x = 1.5
inlet2Velocity_y = 0.
inlet2Velocity_z = 0.

inlet2BufferWidth = inlet2BufferLayers * dp # - dp / 2
inlet2BufferEdge = inlet2BufferPosition_x  + dp / 2 #remove, deprecated
inlet2BufferReferencePoint_x = inlet2BufferPosition_x - inlet2BufferOrientation_x * ( inlet2BufferLayers - 1 ) * dp
inlet2BufferReferencePoint_y = inlet2BufferPosition_y - inlet2BufferOrientation_y * ( inlet2BufferLayers - 1 ) * dp
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

import sys
sys.path.append('../../tools')
import saveParticlesVTK
import numpy as np
import vtk

from vtk.numpy_interface import dataset_adapter as dsa

reader = vtk.vtkPolyDataReader()
reader.SetFileName( f'./template/geometry/fluid_input.vtk' )
reader.ReadAllScalarsOn()
reader.ReadAllVectorsOn()
reader.Update()

polydata = reader.GetOutput()
np_points_fluid = dsa.WrapDataObject( polydata ).Points

fluid_r = np.array( np_points_fluid, dtype=float ) #!!
fluid_v = np.array( dsa.WrapDataObject( polydata ).PointData[ 'Vel' ], dtype=float )
fluid_rho = np.array( dsa.WrapDataObject( polydata ).PointData[ 'Rhop' ] )
fluid_p = np.zeros( len( np_points_fluid ) )
fluid_ptype = np.zeros( len( np_points_fluid ) )

fluidToWrite = saveParticlesVTK.create_pointcloud_polydata( fluid_r, fluid_v, fluid_rho, fluid_p, fluid_ptype )
saveParticlesVTK.save_polydata( fluidToWrite, "sources/openchannel_fluid.vtk" )

reader = vtk.vtkPolyDataReader()
reader.SetFileName( f'./template/geometry/boundary_input.vtk' )
reader.ReadAllScalarsOn()
reader.ReadAllVectorsOn()
reader.Update()

polydata = reader.GetOutput()
np_points_box = dsa.WrapDataObject( polydata ).Points

boundary_r = np.array( np_points_box, dtype=float ) #!!
boundary_v = np.array( dsa.WrapDataObject( polydata ).PointData[ 'Vel' ], dtype=float )
boundary_rho = np.array( dsa.WrapDataObject( polydata ).PointData[ 'Rhop' ] )
boundary_p = np.zeros( len( np_points_box ) )
boundary_ptype = np.zeros( len( np_points_box ) )

boundToWrite = saveParticlesVTK.create_pointcloud_polydata( boundary_r, boundary_v, boundary_rho, boundary_p,
                                                            boundary_ptype )
saveParticlesVTK.save_polydata( boundToWrite, "sources/openchannel_boundary.vtk" )
#---------------------------------------------------------------------------#
### Generate buffer particles
inlet_rx = []; inlet_ry = []; inlet_rz = []
inlet_vx = []; inlet_vy = []; inlet_vz = []
inlet_density = []

inletL_n = inletBufferLayers
inletH_n = round( inletBufferHeight / dp )
inletY_n = round( inletBufferWidthY / dp ) + 1
print( "inletL_n: ", inletBufferLayers )
print( "inletH_n: ", inletH_n )
print( "inletY_n: ", inletY_n )


for x in range( inletL_n ):
    for y in range( inletY_n ):
        for z in range( inletH_n ):
            inlet_rx.append( inletBufferPosition_x - inletBufferOrientation_x * dp * ( x ) )
            inlet_ry.append( inletBufferWidth_start + y * dp )
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

for x in range( inletL_n ):
    for y in range( inletY_n ):
        for z in range( inletH_n ):
            inlet2_rx.append( inlet2BufferPosition_x - inlet2BufferOrientation_x * dp * ( x ) )
            inlet2_ry.append( inletBufferWidth_start + y * dp  )
            inlet2_rz.append( inlet2BufferPosition_z + dp * ( z ) )

            inlet2_vx.append( inlet2Velocity_x )
            inlet2_vy.append( 0. ) #we use only 2D case
            inlet2_vz.append( inlet2Velocity_z )

            hydrostaticPressure = rho0 * 9.81 * ( fluidH - z * dp )
            hydrostaticDensity = ( ( hydrostaticPressure / ( speedOfSound ** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0;
            inlet2_density.append( hydrostaticDensity )

r = np.array( ( inlet_rx, inlet_ry, inlet_rz ), dtype=float ).T #!!
v = np.array( ( inlet_vx, inlet_vy, inlet_vz ), dtype=float ).T #!!
#rho = rho0 * np.ones( len( inlet_rx ) )
rho = np.array( inlet_density, dtype=float )
p = np.zeros( len( inlet_rx ) )
ptype = np.ones( len( inlet_rx ) )

inletToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
saveParticlesVTK.save_polydata( inletToWrite, "sources/openchannel_inlet.vtk" )

r = np.array( ( inlet2_rx, inlet2_ry, inlet2_rz ), dtype=float ).T #!!
v = np.array( ( inlet2_vx, inlet2_vy, inlet2_vz ), dtype=float ).T #!!
#rho = rho0 * np.ones( len( inlet2_rx ) )
rho = np.array( inlet2_density, dtype=float )
p = np.zeros( len( inlet2_rx ) )
ptype = np.ones( len( inlet2_rx ) )

inlet2ToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
saveParticlesVTK.save_polydata( inlet2ToWrite, "sources/openchannel_outlet.vtk" )

#---------------------------------------------------------------------------#
### Compute remaining parameters
spaceDimension = 3 #TODO: Move into templates
particleMass = rho0 * ( dp * dp * dp )
smoothingLentgh =  round( smoothingLentghCoef * dp, 7 )
searchRadius = round( smoothingLentgh * 2 , 7 )

if not timeStep:
    timeStep = round( CFLnumber * ( smoothingLentgh / speedOfSound ), 8 )
coefB = round( speedOfSound * speedOfSound * rho0 / 7 , 1 )

#Determine grid size
import math
gridBegin_x = 1.005 * ( min( min( np_points_fluid[ : , 0 ] ), min( np_points_box[ :, 0 ] ) ) - searchRadius )
gridBegin_y = 1.005 * ( min( min( np_points_fluid[ : , 1 ] ), min( np_points_box[ :, 1 ] ) ) - searchRadius )
gridBegin_z = 1.005 * ( min( min( np_points_fluid[ : , 2 ] ), min( np_points_box[ : ,2 ] ) ) - searchRadius )
gridEnd_x = 1.005 * ( max( max( np_points_fluid[ :, 0 ] ), max( np_points_box[ :, 0 ] ) ) + searchRadius )
gridEnd_y = 1.005 * ( max( max( np_points_fluid[ :, 1 ] ), max( np_points_box[ :, 1 ] ) ) + searchRadius )
gridEnd_z = 1.005 * ( max( max( np_points_fluid[ :, 2 ] ), max( np_points_box[ :, 2 ] ) ) + searchRadius )

gridSize_x = math.ceil( ( gridEnd_x - gridBegin_x ) / searchRadius )
gridSize_y = math.ceil( ( gridEnd_y - gridBegin_y ) / searchRadius )
gridSize_z = math.ceil( 1.2 * ( gridEnd_z - gridBegin_z ) / searchRadius )

### Generate configuration files
# Read in the file
with open( 'template/SPHCaseConfig_template.h', 'r' ) as file :
  fileSPHConf = file.read()

fileSPHConf = fileSPHConf.replace( 'placeholderDimension', str( spaceDimension ) )
fileSPHConf = fileSPHConf.replace( 'placeholderMass', str( round( particleMass, 8 ) ) )
fileSPHConf = fileSPHConf.replace( 'placeholderSpeedOfSound', str( speedOfSound ) )
fileSPHConf = fileSPHConf.replace( 'placeholderCoefB', str( coefB ) )
fileSPHConf = fileSPHConf.replace( 'placeholderDensity', str( rho0 ))
fileSPHConf = fileSPHConf.replace( 'placeholderInitParticleDistance', str( dp ) )
fileSPHConf = fileSPHConf.replace( 'placeholderSmoothingLength', str( smoothingLentgh ) )
fileSPHConf = fileSPHConf.replace( 'placeholderTimeStep', str( timeStep ) )

with open( 'sources/SPHCaseConfig.h', 'w' ) as file:
  file.write( fileSPHConf )

# Setup of particle system
with open( 'template/ParticlesConfig_template.h', 'r' ) as file :
  fileParticleConf = file.read()

fileParticleConf = fileParticleConf.replace( 'placeholderDimension', str( spaceDimension ) )
fileParticleConf = fileParticleConf.replace( 'placeholderFluidParticles', str( len( np_points_fluid ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedFluidParticles', str( len( np_points_fluid ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderBoundaryParticles', str( len( np_points_box ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedBoundaryParticles', str( len( np_points_box ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderInletParticles', str( len( inlet_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedInletParticles', str( len( inlet_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderOutletParticles', str( len( inlet2_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedOutletParticles', str( len( inlet2_rx ) * 3 ) )

fileParticleConf = fileParticleConf.replace( 'placeholderSearchRadius', str( searchRadius ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXSize', str( gridSize_x ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYSize', str( gridSize_y ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridZSize', str( gridSize_z ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXBegin', str( round( gridBegin_x, 8  ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYBegin', str( round( gridBegin_y, 8  ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridZBegin', str( round( gridBegin_y, 8  ) ) )

with open( 'sources/ParticlesConfig.h', 'w' ) as file:
  file.write( fileParticleConf )

with open( 'template/OpenBoundaryConfig_template.h', 'r' ) as file :
  fileOBConf = file.read()

#inlet1
fileOBConf = fileOBConf.replace( 'placeholderInletOrientation_x', str( inletBufferOrientation_x ) )
fileOBConf = fileOBConf.replace( 'placeholderInletOrientation_y', str( inletBufferOrientation_y ) )
fileOBConf = fileOBConf.replace( 'placeholderInletOrientation_z', str( inletBufferOrientation_z ) )
fileOBConf = fileOBConf.replace( 'placeholderInletVelocity_x', str( inletVelocity_x ) )
fileOBConf = fileOBConf.replace( 'placeholderInletVelocity_y', str( inletVelocity_y ) )
fileOBConf = fileOBConf.replace( 'placeholderInletVelocity_z', str( inletVelocity_z ) )
fileOBConf = fileOBConf.replace( 'placeholderInletPosition_x', str( inletBufferPosition_x  + dp/2 ) ) #FIXME
fileOBConf = fileOBConf.replace( 'placeholderInletPosition_y', str( inletBufferPosition_y ) )
fileOBConf = fileOBConf.replace( 'placeholderInletPosition_z', str( inletBufferPosition_z ) )
fileOBConf = fileOBConf.replace( 'placeholderInletDensity', str( rho0 ) )
fileOBConf = fileOBConf.replace( 'placeholderInletWidth_x', str( round( inletBufferWidth, 7 ) ) )
fileOBConf = fileOBConf.replace( 'placeholderInletWidth_y', str( 0. ) )
fileOBConf = fileOBConf.replace( 'placeholderInletWidth_z', str( 0. ) )
fileOBConf = fileOBConf.replace( 'placeholderInletBufferEdge', str( round( inletBufferEdge, 7 ) ) )

#outlet
fileOBConf = fileOBConf.replace( 'placeholderOutletOrientation_x', str( inlet2BufferOrientation_x ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletOrientation_y', str( inlet2BufferOrientation_y ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletOrientation_z', str( inlet2BufferOrientation_z ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletVelocity_x', str( inlet2Velocity_x ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletVelocity_y', str( inlet2Velocity_y ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletVelocity_z', str( inlet2Velocity_z ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletPosition_x', str( inlet2BufferPosition_x - dp/2 ) ) #FIXME
fileOBConf = fileOBConf.replace( 'placeholderOutletPosition_y', str( inlet2BufferPosition_y ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletPosition_z', str( inlet2BufferPosition_z ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletDensity', str( rho0 ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletWidth_x', str( round( inlet2BufferWidth, 7 ) ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletWidth_y', str( 0. ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletWidth_z', str( 0. ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletBufferEdge', str( round( inlet2BufferEdge, 7 ) ) )

# Write the file out again
with open( 'sources/OpenBoundaryConfig.h', 'w' ) as file:
  file.write( fileOBConf )

# Setup of simulation control file
with open( 'template/SimulationControlConfig.h', 'r' ) as file :
  fileSimulationControl = file.read()

with open( 'sources/SimulationControlConfig.h', 'w' ) as file:
  file.write( fileSimulationControl )
