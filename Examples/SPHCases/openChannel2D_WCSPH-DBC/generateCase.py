# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#
# simpleInlet2D_WCSPH-DBC_test
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

## Parameters
dp = 0.002
smoothingLentghCoef = 2**0.5

rho0 = 1000.
p0 = 0.

numberOfBoundaryLayers = 3

speedOfSound = 34.3
CFLnumber = 0.2
timeStep = 0.00002*0.5 #otherwise is obtained automatically

write = '.vtk' #.ptcs or .vtk

boxL = 1.
boxH = 0.8

fluidL = 0.5
fluidH = 0.2

## First inlet buffer. ##
inletBufferOrientation_x = 1.
inletBufferOrientation_z = 0.
inletBufferPosition_x = 0.2
inletBufferPosition_z = 0. + dp*1
inletBufferHeight = 0.1
inletBufferLayers = numberOfBoundaryLayers + 1
inletVelocity_x = 1.
inletVelocity_z = 0.

inletBufferWidth = inletBufferLayers * dp - dp / 2
inletBufferEdge = inletBufferPosition_x + 4 * dp  + dp / 2 #remove, deprecated
inletBufferReferencePoint_x = inletBufferPosition_x - inletBufferOrientation_x * ( inletBufferLayers - 1 ) * dp
inletBufferReferencePoint_z = inletBufferPosition_z - inletBufferOrientation_z * ( inletBufferLayers - 1 ) * dp

## Second inlet buffer. ##
inlet2BufferOrientation_x = -1.
inlet2BufferOrientation_z = 0.
inlet2BufferPosition_x = 0.8
inlet2BufferPosition_z = 0. + dp*1
inlet2BufferHeight = 0.15
inlet2BufferLayers = numberOfBoundaryLayers + 1
inlet2Velocity_x = 1.5
inlet2Velocity_z = 0.

inlet2BufferWidth = inlet2BufferLayers * dp - dp / 2
inlet2BufferEdge = inlet2BufferPosition_x  + dp / 2 #remove, deprecated
inlet2BufferReferencePoint_x = inlet2BufferPosition_x - inlet2BufferOrientation_x * ( inlet2BufferLayers - 1 ) * dp
inlet2BufferReferencePoint_z = inlet2BufferPosition_z - inlet2BufferOrientation_z * ( inlet2BufferLayers - 1 ) * dp

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

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

for x in range( fluidL_n ):
    for z in range( fluidH_n ):
        fluid_rx.append( dp * ( x + 1 ) )
        fluid_ry.append( 0. ) #we use only 2D case
        fluid_rz.append( dp * ( z + 1 ) )

### Generate buffer particles
inlet_rx = []; inlet_ry = []; inlet_rz = []
inlet_vx = []; inlet_vy = []; inlet_vz = []

for x in range( inletL_n ):
    for z in range( inletH_n ):
        inlet_rx.append( inletBufferPosition_x - inletBufferOrientation_x * dp * ( x ) )
        inlet_ry.append( 0. ) #we use only 2D case
        inlet_rz.append( inletBufferPosition_z + dp * ( z ) )

        inlet_vx.append( inletVelocity_x )
        inlet_vy.append( 0. ) #we use only 2D case
        inlet_vz.append( inletVelocity_z )

### Generate buffer particles
inlet2_rx = []; inlet2_ry = []; inlet2_rz = []
inlet2_vx = []; inlet2_vy = []; inlet2_vz = []

for x in range( inlet2L_n ):
    for z in range( inlet2H_n ):
        inlet2_rx.append( inlet2BufferPosition_x - inlet2BufferOrientation_x * dp * ( x ) )
        inlet2_ry.append( 0. ) #we use only 2D case
        inlet2_rz.append( inlet2BufferPosition_z + dp * ( z ) )

        inlet2_vx.append( inlet2Velocity_x )
        inlet2_vy.append( 0. ) #we use only 2D case
        inlet2_vz.append( inlet2Velocity_z )

### Generate boundary particles
box_rx = []; box_ry = []; box_rz = []

# left wall
for layer in range( numberOfBoundaryLayers ):
    for z in range( boxH_n - 1 ):
        box_rx.append( 0. - layer * dp )
        box_ry.append( 0. ) #we use only 2D case
        box_rz.append( ( z+1 ) * dp )

# bottom wall
for layer in range( numberOfBoundaryLayers ):
    for x in range( boxL_n + ( numberOfBoundaryLayers - 1 ) * 2 ):
        box_rx.append( ( x - ( numberOfBoundaryLayers - 1 ) ) * dp )
        box_ry.append( 0. ) #we use only 2D case
        box_rz.append( 0. - layer * dp )

x_last = box_rx[-1 -(numberOfBoundaryLayers - 1)] #due to discretisation, we need to save last value of bottom wall

# right wall
for layer in range( numberOfBoundaryLayers ):
    for z in range( boxH_n - 1 ):
        box_rx.append( x_last + dp * layer )
        box_ry.append( 0. ) #we use only 2D case
        box_rz.append( ( z + 1 ) * dp )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

### Write fluid particles
if write == '.ptcs':
    with open( "openchannel_fluid.ptcs", "w" ) as f:
        f.write( str( len( fluid_rx ) ) + "\n" )
        for i in range( len( fluid_rx ) ):
            f.write( str( round( fluid_rx[ i ], 5 ) ) + " " + str( round( fluid_rz[i], 5 ) ) + " " + \
                     str( round( fluid_ry[ i ], 5 ) ) + " " + str( 0. ) + " " + str( 0. ) + " " + str( 0. ) + " " + \
                     str( round( rho0, 5 ) ) + " " + str( round( p0, 5 ) ) + " " + str( 0 ) + "\n" )
    ### Write boundary particles
    with open("openchannel_boundary.ptcs", "w") as f:
        f.write( str( len( box_rx ) ) + "\n" )
        for i in range( len( box_rx ) ):
            f.write( str( round( box_rx[ i ], 5 ) ) + " " + str( round( box_rz[ i ], 5 ) ) + " " + \
                     str( round( box_ry[ i ], 5 ) ) + " " + str( 0. ) + " " + str( 0. ) + " " + str( 0. ) + " " + \
                     str( round( rho0, 5 ) ) + " " + str( round( p0, 5 ) ) + " " + str( 1 ) + "\n" )
    ### Write fluid particles
    with open("openchannel_inlet.ptcs", "w") as f:
        f.write( str( len( inlet_rx ) ) + "\n" )
        for i in range( len( inlet_rx ) ):
            f.write( str( round( inlet_rx[ i ], 5 ) ) + " " + str( round( inlet_rz[ i ], 5 ) ) + " " + \
                     str( round( inlet_ry[ i ], 5 ) ) + " " + str( round( inlet_vx[ i ] ) ) + " " + str( round( inlet_vy[ i ] ) ) + " " + str( round( inlet_vz[ i ] ) ) + " " + \
                     str( round( rho0, 5 ) ) + " " + str( round( p0, 5 ) ) + " " + str( 10 ) + "\n" )
    ### Write fluid particles
    with open("openchannel_inlet2.ptcs", "w") as f:
        f.write( str( len( inlet2_rx ) ) + "\n" )
        for i in range( len( inlet2_rx ) ):
            f.write( str( round( inlet2_rx[ i ], 5 ) ) + " " + str( round( inlet2_rz[ i ], 5 ) ) + " " + \
                     str( round( inlet2_ry[ i ], 5 ) ) + " " + str( round( inlet2_vx[ i ] ) ) + " " + str( round( inlet2_vy[ i ] ) ) + " " + str( round( inlet2_vz[ i ] ) ) + " " + \
                     str( round( rho0, 5 ) ) + " " + str( round( p0, 5 ) ) + " " + str( 10 ) + "\n" )
elif write == '.vtk':
    import sys
    sys.path.append('../../tools/')
    import saveParticlesVTK
    import numpy as np

    #r = np.array( ( fluid_rx, fluid_rz, fluid_ry ), dtype=float ).T #!!
    #v = np.zeros( ( len( fluid_rx ), 3 ) )
    #rho = rho0 * np.ones( len( fluid_rx ) )
    #p = np.zeros( len( fluid_rx ) )
    #ptype = np.zeros( len( fluid_rx ) )

    #fluidToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
    #saveParticlesVTK.save_polydata( fluidToWrite, "openchannel_fluid.vtk" )

    r = np.array( ( box_rx, box_rz, box_ry ), dtype=float ).T #!!
    v = np.zeros( ( len( box_rx ), 3 ) )
    rho = rho0 * np.ones( len( box_rx ) )
    p = np.zeros( len( box_rx ) )
    ptype = np.ones( len( box_rx ) )

    boxToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
    saveParticlesVTK.save_polydata( boxToWrite, "openchannel_boundary.vtk" )

    r = np.array( ( inlet_rx, inlet_rz, inlet_ry ), dtype=float ).T #!!
    v = np.array( ( inlet_vx, inlet_vz, inlet_vy ), dtype=float ).T #!!
    rho = rho0 * np.ones( len( inlet_rx ) )
    p = np.zeros( len( inlet_rx ) )
    ptype = np.ones( len( inlet_rx ) )

    inletToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
    saveParticlesVTK.save_polydata( inletToWrite, "openchannel_inlet.vtk" )

    r = np.array( ( inlet2_rx, inlet2_rz, inlet2_ry ), dtype=float ).T #!!
    v = np.array( ( inlet2_vx, inlet2_vz, inlet2_vy ), dtype=float ).T #!!
    rho = rho0 * np.ones( len( inlet2_rx ) )
    p = np.zeros( len( inlet2_rx ) )
    ptype = np.ones( len( inlet2_rx ) )

    inlet2ToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
    saveParticlesVTK.save_polydata( inlet2ToWrite, "openchannel_outlet.vtk" )
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
with open( 'SPHCaseConfig.h', 'w' ) as file:
  file.write( fileSPHConf )

with open( 'template/ParticlesConfig_template.h', 'r' ) as file :
  fileParticleConf = file.read()

# Replace the target string
fileParticleConf = fileParticleConf.replace( 'placeholderDimension', str( spaceDimension ) )
fileParticleConf = fileParticleConf.replace( 'placeholderFluidParticles', str( 0 ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedFluidParticles', str( len( fluid_rx ) * 2 ) )
fileParticleConf = fileParticleConf.replace( 'placeholderBoundaryParticles', str( len( box_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedBoundaryParticles', str( len( box_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderInletParticles', str( len( inlet_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedInletParticles', str( len( inlet_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderOutletParticles', str( len( inlet2_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedOutletParticles', str( len( inlet2_rx ) ) )

fileParticleConf = fileParticleConf.replace( 'placeholderSearchRadius', str( searchRadius ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXSize', str( gridXsize ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYSize', str( gridYsize ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXBegin', str( round( gridXbegin, 9  ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYBegin', str( round( gridYbegin, 9  ) ) )

# Write the file out again
with open( 'ParticlesConfig.h', 'w' ) as file:
  file.write( fileParticleConf )

with open( 'template/OpenBoundaryConfig_template.h', 'r' ) as file :
  fileOBConf = file.read()

#inlet1
fileOBConf = fileOBConf.replace( 'placeholderOBP1Orientation_x', str( inletBufferOrientation_x ) )
fileOBConf = fileOBConf.replace( 'placeholderOBP1Orientation_y', str( inletBufferOrientation_z ) )
fileOBConf = fileOBConf.replace( 'placeholderOBP1Velocity_x', str( inletVelocity_x ) )
fileOBConf = fileOBConf.replace( 'placeholderOBP1Velocity_y', str( inletVelocity_z ) )
fileOBConf = fileOBConf.replace( 'placeholderOBP1Position_x', str( inletBufferPosition_x ) )
fileOBConf = fileOBConf.replace( 'placeholderOBP1Position_y', str( inletBufferPosition_z ) )
fileOBConf = fileOBConf.replace( 'placeholderOBP1Density', str( rho0 ) )
fileOBConf = fileOBConf.replace( 'placeholderOBP1Width_x', str( round( inletBufferWidth, 7 ) ) )
fileOBConf = fileOBConf.replace( 'placeholderOBP1Width_y', str( 0. ) )
fileOBConf = fileOBConf.replace( 'placeholderOBP1BufferEdge', str( round(  inletBufferEdge, 7 ) ) )

#outlet
fileOBConf = fileOBConf.replace( 'placeholderOBP2Orientation_x', str( inlet2BufferOrientation_x ) )
fileOBConf = fileOBConf.replace( 'placeholderOBP2Orientation_y', str( inlet2BufferOrientation_z ) )
fileOBConf = fileOBConf.replace( 'placeholderOBP2Velocity_x', str( inlet2Velocity_x ) )
fileOBConf = fileOBConf.replace( 'placeholderOBP2Velocity_y', str( inlet2Velocity_z ) )
fileOBConf = fileOBConf.replace( 'placeholderOBP2Position_x', str( inlet2BufferPosition_x ) )
fileOBConf = fileOBConf.replace( 'placeholderOBP2Position_y', str( inlet2BufferPosition_z ) )
fileOBConf = fileOBConf.replace( 'placeholderOBP2Density', str( rho0 ) )
fileOBConf = fileOBConf.replace( 'placeholderOBP2Width_x', str( round( inlet2BufferWidth, 7 ) ) )
fileOBConf = fileOBConf.replace( 'placeholderOBP2Width_y', str( 0. ) )
fileOBConf = fileOBConf.replace( 'placeholderOBP2BufferEdge', str( round(  inlet2BufferEdge, 7 ) ) )

# Write the file out again
with open( 'OpenBoundaryConfig.h', 'w' ) as file:
  file.write( fileOBConf )

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

import os
resultsPath = r'./results'
if not os.path.exists( resultsPath ):
    os.makedirs( resultsPath )
