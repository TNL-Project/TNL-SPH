#---------------------------------------------------------------------------#
#
# case: openChannelWithObstacle (without obstacle atm)
#
#---------------------------------------------------------------------------#
### Parameters of the case necessary for case creation:

# Dimensions of the channel
fluidL = 0.5
fluidH = 0.1

# Obstacle position and size
obstacleCenterX = 0.2
obstacleCenterY = 0.07
obstacleR = 0.01

# Initial particle distance (dp)[m]:
dp = 0.002

# Smoothing length coefitient:
# smoothing length (h) = smoothing length coef (Coef_h) * initial particle distance (d_p)
# [ h = Coef_h * dp ]
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

numberOfAllocatedParticles = 100000

## Define inlet
inletOrientation_x = 1.
inletOrientation_z = 0.
inletPosition_x = 0.0
inletPosition_z = 0. + dp*1
inletHeight = fluidH
inletLayers = numberOfBoundaryLayers + 1
inletVelocity_x = 0.5
inletVelocity_z = 0.

inletWidth = inletLayers * dp
inletEdge = inletPosition_x + 4 * dp  + dp / 2 #remove, deprecated
inletReferencePoint_x = inletPosition_x - inletOrientation_x * ( inletLayers - 1 ) * dp
inletReferencePoint_z = inletPosition_z - inletOrientation_z * ( inletLayers - 1 ) * dp

## Define outlet
ouletOrientation_x = -1.
ouletOrientation_z = 0.
ouletPosition_x = fluidL + dp
ouletPosition_z = 0. + dp*1
ouletHeight = fluidH
ouletLayers = numberOfBoundaryLayers + 1
outletVelocity_x = 1.5
outletVelocity_z = 0.

ouletWidth = ouletLayers * dp
ouletEdge = ouletPosition_x  + dp / 2 #remove, deprecated
ouletReferencePoint_x = ouletPosition_x - ouletOrientation_x * ( ouletLayers - 1 ) * dp
ouletReferencePoint_z = ouletPosition_z - ouletOrientation_z * ( ouletLayers - 1 ) * dp

#---------------------------------------------------------------------------#

### Create related directories
import os
resultsPath = r'./results'
if not os.path.exists( resultsPath ):
    os.makedirs( resultsPath )

sourcesPath = r'./sources'
if not os.path.exists( sourcesPath ):
    os.makedirs( sourcesPath )

### Generate fluid particles
fluid_rx = []; fluid_ry = []; fluid_rz = []
fluid_density = []

fluidL_n = round( fluidL / dp )
fluidH_n = round( fluidH / dp )

for x in range( fluidL_n ):
    for y in range( fluidH_n ):

        ptcs_x = inletPosition_x + dp * ( x + 1 )
        ptcs_y = dp * ( y + 1 )
        r  = ( ( ptcs_x - obstacleCenterX )**2  + ( ptcs_y - obstacleCenterY )**2 )**0.5

        if ( r  > obstacleR + dp * 0.5 ):
            fluid_rx.append( ptcs_x )
            fluid_ry.append( ptcs_y )
            fluid_rz.append( 0. )

            hydrostaticPressure = rho0 * 9.81 * ( fluidH - y * dp )
            hydrostaticDensity = ( ( hydrostaticPressure / ( speedOfSound ** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0
            fluid_density.append( hydrostaticDensity )

### Generate buffer particles
inlet_rx = []; inlet_ry = []; inlet_rz = []
inlet_vx = []; inlet_vy = []; inlet_vz = []
inlet_density = []

inletL_n = inletLayers
inletH_n = round( inletHeight / dp  )

for x in range( inletL_n ):
    for y in range( inletH_n ):
        inlet_rx.append( inletPosition_x - inletOrientation_x * dp * ( x ) )
        inlet_ry.append( inletPosition_z + dp * ( y ) )
        inlet_rz.append( 0. )

        inlet_vx.append( inletVelocity_x )
        inlet_vy.append( inletVelocity_z )
        inlet_vz.append( 0. )

        hydrostaticPressure = rho0 * 9.81 * ( fluidH - y * dp )
        hydrostaticDensity = ( ( hydrostaticPressure / ( speedOfSound ** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0
        inlet_density.append( hydrostaticDensity )

### Generate buffer particles
outlet_rx = []; outlet_ry = []; outlet_rz = []
outlet_vx = []; outlet_vy = []; outlet_vz = []
outlet_density = []

outletL_n = ouletLayers
outletH_n = round( ouletHeight / dp  )

for x in range( outletL_n ):
    for y in range( outletH_n ):
        outlet_rx.append( ouletPosition_x - ouletOrientation_x * dp * ( x ) )
        outlet_ry.append( ouletPosition_z + dp * ( y ) )
        outlet_rz.append( 0. )

        outlet_vx.append( outletVelocity_x )
        outlet_vy.append( outletVelocity_z )
        outlet_vz.append( 0. )

        hydrostaticPressure = rho0 * 9.81 * ( fluidH - y * dp )
        hydrostaticDensity = ( ( hydrostaticPressure / ( speedOfSound ** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0
        outlet_density.append( hydrostaticDensity )

### Generate boundary particles
box_rx = []; box_ry = []; box_rz = []
ghost_rx = []; ghost_ry = []; ghost_rz = []
box_density = []

boxL_n = round( boxL / dp )
boxH_n = round( boxH / dp )

# bottom wall
for layer in range( numberOfBoundaryLayers ):
    for x in range( boxL_n + ( numberOfBoundaryLayers - 1 ) * 2 ):
        box_rx.append( ( x - ( numberOfBoundaryLayers - 1 ) ) * dp )
        box_ry.append( 0. - layer * dp )
        box_rz.append( 0. )

        ghost_rx.append( ( x - ( numberOfBoundaryLayers - 1 ) ) * dp )
        ghost_ry.append( 0. + dp * ( layer + 1 ) )
        ghost_rz.append( 0.)

        hydrostaticPressure = rho0 * 9.81 * ( fluidH + layer * dp )
        hydrostaticDensity = ( ( hydrostaticPressure / ( speedOfSound ** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0
        box_density.append( hydrostaticDensity )

# obstacle
import numpy as np
for layer in range( numberOfBoundaryLayers ):
    r_layer = obstacleR - layer * dp
    obstacleR_n = int(np.pi * 2 * r_layer / dp)
    angleDP = 2 * np.pi / obstacleR_n

    print("r_layer: ", r_layer)
    print("obstacleR_n: ", obstacleR_n)
    print("angleDP: ", angleDP)

    for phi in range( obstacleR_n ):

        angle = angleDP * phi
        #print( "angle", angle )

        box_rx.append( np.cos( angle ) * r_layer + obstacleCenterX )
        box_ry.append( np.sin( angle ) * r_layer + obstacleCenterY )
        box_rz.append( 0. )

        ghost_rx.append( ( x - ( numberOfBoundaryLayers - 1 ) ) * dp )
        ghost_ry.append( 0. + dp * ( layer + 1 ) )
        ghost_rz.append( 0.)

        hydrostaticPressure = rho0 * 9.81 * ( fluidH + layer * dp )
        hydrostaticDensity = ( ( hydrostaticPressure / ( speedOfSound ** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0
        box_density.append( hydrostaticDensity )

import sys
sys.path.append('../../tools/')
import saveParticlesVTK

fluid_r = np.array( ( fluid_rx, fluid_ry, fluid_rz ), dtype=float ).T #!!
fluid_v = np.zeros( ( len( fluid_rx ), 3 ) )
fluid_rho = np.array( fluid_density, dtype=float )
fluid_p = np.zeros( len( fluid_rx ) )
fluid_ptype = np.zeros( len( fluid_rx ) )

fluidToWrite = saveParticlesVTK.create_pointcloud_polydata( fluid_r, fluid_v, fluid_rho, fluid_p, fluid_ptype )
saveParticlesVTK.save_polydata( fluidToWrite, "sources/openchannel_fluid.vtk" )

boundary_r = np.array( ( box_rx, box_ry, box_rz ), dtype=float ).T #!!
boundary_ghostNodes = np.array( ( ghost_rx, ghost_ry, ghost_rz ), dtype=float ).T #!!
boundary_v = np.zeros( ( len( box_rx ), 3 ) )
boundary_rho = np.array( box_density, dtype=float )
boundary_p = np.zeros( len( box_rx ) )
boundary_ptype = np.ones( len( box_rx ) )

boxToWrite = saveParticlesVTK.create_pointcloud_polydata( boundary_r, boundary_v, boundary_rho, boundary_p, boundary_ptype,
                                                          ghostNodes=boundary_ghostNodes )
saveParticlesVTK.save_polydata( boxToWrite, "sources/openchannel_boundary.vtk" )

inlet_r = np.array( ( inlet_rx, inlet_ry, inlet_rz ), dtype=float ).T #!!
inlet_v = np.array( ( inlet_vx, inlet_vz, inlet_vz ), dtype=float ).T #!!
inlet_rho = np.array( inlet_density, dtype=float )
inlet_p = np.zeros( len( inlet_rx ) )
inlet_ptype = np.ones( len( inlet_rx ) )

inletToWrite = saveParticlesVTK.create_pointcloud_polydata( inlet_r, inlet_v, inlet_rho, inlet_p, inlet_ptype )
saveParticlesVTK.save_polydata( inletToWrite, "sources/openchannel_inlet.vtk" )

outlet_r = np.array( ( outlet_rx, outlet_ry, outlet_rz ), dtype=float ).T #!!
outlet_v = np.array( ( outlet_vx, outlet_vy, outlet_vz ), dtype=float ).T #!!
outlet_rho = np.array( outlet_density, dtype=float )
outlet_p = np.zeros( len( outlet_rx ) )
outlet_ptype = np.ones( len( outlet_rx ) )

outletToWrite = saveParticlesVTK.create_pointcloud_polydata( outlet_r, outlet_v, outlet_rho, outlet_p, outlet_ptype )
saveParticlesVTK.save_polydata( outletToWrite, "sources/openchannel_outlet.vtk" )

### Compute remaining parameters
spaceDimension = 2
particleMass = rho0 * ( dp * dp )
smoothingLentgh =  round( smoothingLentghCoef * dp, 7 )
searchRadius = round( smoothingLentgh * 2 , 7 )

if not timeStep:
    timeStep = round( CFLnumber * ( smoothingLentgh / speedOfSound ), 8 )
coefB = round( speedOfSound * speedOfSound * rho0 / 7 , 1 )

#Determine grid size
from math import ceil
gridBegin_x = 1.005 * ( min( min( fluid_rx ), min( box_rx ) )  - searchRadius )
gridBegin_y = 1.005 * ( min( min( fluid_ry ), min( box_ry ) )  - searchRadius )
gridEnd_x = 1.005 * ( max( max( fluid_rx ), max( box_rx ) )  + searchRadius )
gridEnd_y = 1.2 * ( max( max( fluid_ry ), max( box_ry ) )  + searchRadius )

gridSize_x = ceil( ( gridEnd_x - gridBegin_x ) / searchRadius )
gridSize_y = ceil( ( gridEnd_y - gridBegin_y ) / searchRadius )

### Generate configuration files
# SPH parameters
with open( 'template/SPHCaseConfig_template.h', 'r' ) as file :
  fileSPHConf = file.read()

fileSPHConf = fileSPHConf.replace( 'placeholderDimension', str( spaceDimension ) )
fileSPHConf = fileSPHConf.replace( 'placeholderMass', str( particleMass ) )
fileSPHConf = fileSPHConf.replace( 'placeholderSpeedOfSound', str( speedOfSound ) )
fileSPHConf = fileSPHConf.replace( 'placeholderCoefB', str( coefB ) )
fileSPHConf = fileSPHConf.replace( 'placeholderDensity', str( rho0 ))
fileSPHConf = fileSPHConf.replace( 'placeholderInitParticleDistance', str( dp ) )
fileSPHConf = fileSPHConf.replace( 'placeholderSmoothingLength', str( smoothingLentgh ) )
fileSPHConf = fileSPHConf.replace( 'placeholderTimeStep', str( timeStep ) )

with open( 'sources/SPHCaseConfig.h', 'w' ) as file:
  file.write( fileSPHConf )

# Setup for particle system
with open( 'template/ParticlesConfig_template.h', 'r' ) as file :
  fileParticleConf = file.read()

fileParticleConf = fileParticleConf.replace( 'placeholderDimension', str( spaceDimension ) )
fileParticleConf = fileParticleConf.replace( 'placeholderFluidParticles', str( len( fluid_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedFluidParticles', str( numberOfAllocatedParticles ) )
fileParticleConf = fileParticleConf.replace( 'placeholderBoundaryParticles', str( len( box_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedBoundaryParticles', str( len( box_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderInletParticles', str( len( inlet_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedInletParticles', str( len( inlet_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderOutletParticles', str( len( outlet_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedOutletParticles', str( len( outlet_rx ) * 3 ) )

fileParticleConf = fileParticleConf.replace( 'placeholderSearchRadius', str( searchRadius ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXSize', str( gridSize_x ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYSize', str( gridSize_y ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXBegin', str( round( gridBegin_x, 9  ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYBegin', str( round( gridBegin_y, 9  ) ) )

with open( 'sources/ParticlesConfig.h', 'w' ) as file:
  file.write( fileParticleConf )

# Open boundary parameters
with open( 'template/OpenBoundaryConfig_template.h', 'r' ) as file :
  fileOBConf = file.read()

fileOBConf = fileOBConf.replace( 'placeholderInletOrientation_x', str( inletOrientation_x ) )
fileOBConf = fileOBConf.replace( 'placeholderInletOrientation_y', str( inletOrientation_z ) )
fileOBConf = fileOBConf.replace( 'placeholderInletVelocity_x', str( inletVelocity_x ) )
fileOBConf = fileOBConf.replace( 'placeholderInletVelocity_y', str( inletVelocity_z ) )
fileOBConf = fileOBConf.replace( 'placeholderInletPosition_x', str( inletPosition_x  + dp/2 ) ) #FIXME
fileOBConf = fileOBConf.replace( 'placeholderInletPosition_y', str( inletPosition_z ) )
fileOBConf = fileOBConf.replace( 'placeholderInletDensity', str( rho0 ) )
fileOBConf = fileOBConf.replace( 'placeholderInletWidth_x', str( round( inletWidth, 7 ) ) )
fileOBConf = fileOBConf.replace( 'placeholderInletWidth_y', str( 0. ) )
fileOBConf = fileOBConf.replace( 'placeholderInletBufferEdge', str( round(  inletEdge, 7 ) ) )

fileOBConf = fileOBConf.replace( 'placeholderOutletOrientation_x', str( ouletOrientation_x ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletOrientation_y', str( ouletOrientation_z ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletVelocity_x', str( outletVelocity_x ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletVelocity_y', str( outletVelocity_z ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletPosition_x', str( ouletPosition_x - dp/2 ) ) #FIXME
fileOBConf = fileOBConf.replace( 'placeholderOutletPosition_y', str( ouletPosition_z ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletDensity', str( rho0 ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletWidth_x', str( round( ouletWidth, 7 ) ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletWidth_y', str( 0. ) )
fileOBConf = fileOBConf.replace( 'placeholderOutletBufferEdge', str( round(  ouletEdge, 7 ) ) )

with open( 'sources/OpenBoundaryConfig.h', 'w' ) as file:
  file.write( fileOBConf )

# Setup for simulation control file
with open( 'template/SimulationControlConfig.h', 'r' ) as file :
  fileSimulationControl = file.read()

with open( 'sources/SimulationControlConfig.h', 'w' ) as file:
  file.write( fileSimulationControl )

# Setup the measuretool config
with open( 'template/MeasuretoolConfig.h', 'r' ) as file :
  fileMeasuretoolConf = file.read()

with open( 'sources/MeasuretoolConfig.h', 'w' ) as file:
  file.write( fileMeasuretoolConf )

### Write the domain grid (this is not necessary for the computation)
from domainGrid import DomainGrid

DomainGrid( gridSize_x, gridSize_y, 1,                  # grid size
            gridBegin_x, gridBegin_y, 0,                # coordinates of grid origin
            np.zeros( gridSize_x * gridSize_y ),        # array with index of grid sector
            searchRadius,                               # search radius i. e. cell size
            'sources/openchannel_grid.vtk' )            # outputfile name
