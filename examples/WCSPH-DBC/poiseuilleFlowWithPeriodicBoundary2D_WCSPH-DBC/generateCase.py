#---------------------------------------------------------------------------#
#
# case:
#
#---------------------------------------------------------------------------#
### Parameters of the case necessary for case creation:

# Channel width H [m]:
channelH = 0.1

#Channel length L [m]
channelL = 0.1

#Inlet velocity vx_0 [m/s]:
inletVelocity = 0.5

import argparse

parser = argparse.ArgumentParser()
parser.add_argument( "-resolution", default=0.002, type=float )
args = parser.parse_args()

# Initial particle distance (dp):
dp =  args.resolution

# Smoothing length coefitient:
# smoothing length (h) = smoothing length coef (Coef_h) * initial particle distance (d_p)
# [ h = Coef_h * dp ]
smoothingLentghCoef = 2**0.5

# Referential density of the medium (rho0):
rho0 = 1000.

# Numerical speed of sound (c0):
speedOfSound = 34.3

# Number of boundary layers:
numberOfBoundaryLayers = 3

# Initial time step.
# In case, that initial time step is not defined, is computed automatically.
timeStep = 0.00001

# CFL number (CFL):
CFLnumber = 0.2

#Number of allocated particles (necessary for inlet)
# In case, that number of allocated particles is not defined, is computed automatically.
numberOfAllocatedParticles = 10000
#---------------------------------------------------------------------------#

boxL = 0.1
boxH = 0.3

fluidL = 0.1 - dp
fluidH = 0.1 - dp

# Inlet buffer parameters
inletBufferOrientation_x = 1.
inletBufferOrientation_z = 0.
inletBufferPosition_x = 0.0
inletBufferPosition_z = 0. + dp #right bottom corner as referential point
inletBufferHeight = channelH - dp
inletBufferLayers = numberOfBoundaryLayers + 1
inletVelocity_x = inletVelocity
inletVelocity_z = 0.

inletBufferWidth = inletBufferLayers * dp # - dp / 2
inletBufferEdge = inletBufferPosition_x + 4 * dp  + dp / 2 #remove, deprecated
inletBufferReferencePoint_x = inletBufferPosition_x - inletBufferOrientation_x * ( inletBufferLayers - 1 ) * dp
inletBufferReferencePoint_z = inletBufferPosition_z - inletBufferOrientation_z * ( inletBufferLayers - 1 ) * dp

## Second inlet buffer. ##
inlet2BufferOrientation_x = -1.
inlet2BufferOrientation_z = 0.
inlet2BufferPosition_x = channelL
inlet2BufferPosition_z = 0. + dp #left bottom corner as referential point
inlet2BufferHeight = channelH - dp
inlet2BufferLayers = numberOfBoundaryLayers + 1
inlet2Velocity_x = inletVelocity #initialize with inlet velocity
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


### Generate fluid particles
fluid_rx = []; fluid_ry = []; fluid_rz = []
fluid_density = []

fluidL_n = round( fluidL / dp )
fluidH_n = round( fluidH / dp )

for x in range( fluidL_n ):
    for y in range( fluidH_n ):
        fluid_rx.append( inletBufferPosition_x + dp * ( x + 1 ) )
        fluid_ry.append( dp * ( y + 1 ) )
        fluid_rz.append( 0. )

        hydrostaticPressure = rho0 * 9.81 * ( fluidH - y * dp )
        hydrostaticDensity = ( ( hydrostaticPressure / ( speedOfSound ** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0
        fluid_density.append( rho0 )

### Generate buffer particles
inlet_rx = []; inlet_ry = []; inlet_rz = []
inlet_vx = []; inlet_vy = []; inlet_vz = []
inlet_density = []

inletL_n = inletBufferLayers
inletH_n = round( inletBufferHeight / dp  )

for x in range( inletL_n ):
    for y in range( inletH_n ):
        inlet_rx.append( inletBufferPosition_x - inletBufferOrientation_x * dp * ( x ) )
        inlet_ry.append( inletBufferPosition_z + dp * ( y ) )
        inlet_rz.append( 0. )

        inlet_vx.append( inletVelocity_x )
        inlet_vy.append( inletVelocity_z )
        inlet_vz.append( 0. )

        hydrostaticPressure = rho0 * 9.81 * ( fluidH - y * dp )
        hydrostaticDensity = ( ( hydrostaticPressure / ( speedOfSound ** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0
        inlet_density.append( rho0 )

### Generate buffer particles
inlet2_rx = []; inlet2_ry = []; inlet2_rz = []
inlet2_vx = []; inlet2_vy = []; inlet2_vz = []
inlet2_density = []

inlet2L_n = inlet2BufferLayers
inlet2H_n = round( inlet2BufferHeight / dp  )

for x in range( inlet2L_n ):
    for y in range( inlet2H_n ):
        inlet2_rx.append( inlet2BufferPosition_x - inlet2BufferOrientation_x * dp * ( x ) )
        inlet2_ry.append( inlet2BufferPosition_z + dp * ( y ) )
        inlet2_rz.append( 0. )

        inlet2_vx.append( inlet2Velocity_x )
        inlet2_vy.append( inlet2Velocity_z )
        inlet2_vz.append( 0. )

        hydrostaticPressure = rho0 * 9.81 * ( fluidH - y * dp )
        hydrostaticDensity = ( ( hydrostaticPressure / ( speedOfSound ** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0
        inlet2_density.append( rho0 )

### Generate boundary particles
box_rx = []; box_ry = []; box_rz = []
ghost_rx = []; ghost_ry = []; ghost_rz = []
box_density = []

boxL_n = round( boxL / dp )
boxH_n = round( boxH / dp )


#:# left wall
#:for layer in range( numberOfBoundaryLayers ):
#:    for z in range( boxH_n - 1 ):
#:        box_rx.append( 0. - layer * dp )
#:        box_ry.append( 0. ) #we use only 2D case
#:        box_rz.append( ( z+1 ) * dp )
#:        box_density.append( rho0 );

# bottom wall
for layer in range( numberOfBoundaryLayers ):
    for x in range( boxL_n + ( inletBufferLayers - 1 ) * 2 + 1):
        #box_rx.append( ( x - ( numberOfBoundaryLayers - 1 ) ) * dp )
        box_rx.append( ( x - ( inletBufferLayers - 1 ) ) * dp )
        box_ry.append( 0. - layer * dp )
        box_rz.append( 0. )

        ghost_rx.append( ( x - ( numberOfBoundaryLayers - 1 ) ) * dp )
        ghost_ry.append( 0. + dp * ( layer + 1 ) )
        ghost_rz.append( 0.)

        #hydrostaticPressure = rho0 * 9.81 * ( fluidH - z * dp )
        hydrostaticPressure = rho0 * 9.81 * ( fluidH + layer * dp )
        hydrostaticDensity = ( ( hydrostaticPressure / ( speedOfSound ** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0
        box_density.append( rho0 )

# top wall
for layer in range( numberOfBoundaryLayers ):
    for x in range( boxL_n + ( inletBufferLayers - 1 ) * 2 + 1 ):
        #box_rx.append( ( x - ( numberOfBoundaryLayers - 1 ) ) * dp )
        box_rx.append( ( x - ( inletBufferLayers - 1 ) ) * dp )
        box_ry.append( ( fluidH + dp ) + layer * dp )
        box_rz.append( 0. ) #we use only 2D case

        ghost_rx.append( ( x - ( numberOfBoundaryLayers - 1 ) ) * dp )
        ghost_ry.append( 0. + dp * ( layer + 1 ) )
        ghost_rz.append( 0.)

        #hydrostaticPressure = rho0 * 9.81 * ( fluidH - z * dp )
        hydrostaticPressure = rho0 * 9.81 * ( fluidH + layer * dp )
        hydrostaticDensity = ( ( hydrostaticPressure / ( speedOfSound ** 2 * rho0 / 7 ) + 1 )**( 1./7. ) )  * rho0
        box_density.append( rho0 )

import sys
sys.path.append('../../../src/tools/')
import saveParticlesVTK
import numpy as np

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

boxToWrite = saveParticlesVTK.create_pointcloud_polydata( boundary_r, boundary_v, boundary_rho, boundary_p, boundary_ptype, ghostNodes=boundary_ghostNodes )
saveParticlesVTK.save_polydata( boxToWrite, "sources/openchannel_boundary.vtk" )

r = np.array( ( inlet_rx, inlet_ry, inlet_rz ), dtype=float ).T #!!
v = np.array( ( inlet_vx, inlet_vy, inlet_vz ), dtype=float ).T #!!
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
gridYbegin = 1.005 * ( min( min( fluid_ry ), min( box_ry ) ) ) - searchRadius

gridXend = 1.005 * ( max( max( fluid_rx ), max( box_rx ) ) ) + searchRadius
gridYend = 1.005 * ( max( max( fluid_ry ), max( box_ry ) ) ) + searchRadius * 10 #TODO: TEMP FIX

gridXsize = math.ceil( ( gridXend - gridXbegin ) / searchRadius )
gridYsize = math.ceil( ( gridYend - gridYbegin ) / searchRadius )

#Determine parameters of inlet1
leftShiftVector_x = inlet2BufferPosition_x - inletBufferPosition_x - dp
rightShiftVector_x = ( -1 ) * leftShiftVector_x
numberOfParticlesPerCell = 15

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

with open( 'template/PeriodicBoundaryConfig_template.h', 'r' ) as file :
  fileOBConf = file.read()

#inlet1
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityLeftOrientation_x', str( inletBufferOrientation_x ) )
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityLeftOrientation_y', str( inletBufferOrientation_z ) )
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityLeftPosition_x', str( inletBufferPosition_x  + dp/2 ) ) #FIXME
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityLeftPosition_y', str( inletBufferPosition_z ) )
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityLeftWidth_x', str( round( inletBufferWidth, 7 ) ) )
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityLeftWidth_y', str( 0. ) )
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityLeftHeigth_x', str( 0. ) )
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityLeftHeigth_y', str( round( inletBufferHeight, 7 ) ) )
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityLeftShiftVector_x', str( round( leftShiftVector_x, 7 ) ) )
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityLeftShiftVector_y', str( 0. ) )

fileOBConf = fileOBConf.replace( 'placeholderPeriodicityLeftFirstPoint_x', str( round( inlet2BufferPosition_x, 7 ) ) )
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityLeftFirstPoint_y', str( round( inlet2BufferPosition_z, 7 ) ) )
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityLeftSecondPoint_x', str(round( inlet2BufferPosition_x - dp/2 + searchRadius * 1.1 , 7 ) ) )
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityLeftSecondPoint_y', str( round( inlet2BufferPosition_z, 7 ) ) )

#outlet
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityRightOrientation_x', str( inlet2BufferOrientation_x ) )
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityRightOrientation_y', str( inlet2BufferOrientation_z ) )
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityRightVelocity_x', str( inlet2Velocity_x ) )
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityRightVelocity_y', str( inlet2Velocity_z ) )
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityRightPosition_x', str( inlet2BufferPosition_x - dp/2 ) ) #FIXME
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityRightPosition_y', str( inlet2BufferPosition_z ) )
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityRightDensity', str( rho0 ) )
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityRightWidth_x', str( round( inlet2BufferWidth, 7 ) ) )
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityRightWidth_y', str( 0. ) )
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityRightHeigth_x', str( 0. ) )
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityRightHeigth_y', str( round( inlet2BufferHeight, 7 ) ) )
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityRightShiftVector_x', str( round( rightShiftVector_x, 7 ) ) )
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityRightShiftVector_y', str( 0. ) )

fileOBConf = fileOBConf.replace( 'placeholderPeriodicityRightFirstPoint_x', str( round( inlet2BufferPosition_x - dp/2 - searchRadius * 1.1 , 7 ) ) )
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityRightFirstPoint_y', str( round( inlet2BufferPosition_z, 7 ) ) )
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityRightSecondPoint_x', str( inlet2BufferPosition_x - dp/2 ) )
fileOBConf = fileOBConf.replace( 'placeholderPeriodicityRightSecondPoint_y', str( round( inlet2BufferPosition_z, 7 ) ) )

# BOTH
fileOBConf = fileOBConf.replace( 'placeholderNumberOfParticlesPerCell', str( numberOfParticlesPerCell ) )

# Write the file out again
with open( 'sources/PeriodicBoundaryConfig.h', 'w' ) as file:
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

# Save domain grid.
searchRadius_h = round( smoothingLentgh * 2 , 7 )
from contextlib import redirect_stdout

def DomainGrid( gridXsize, gridYsize, gridZsize, gridXbegin, gridYbegin, gridZbegin, gridSector, name ):
    with open( name, 'w' ) as f:
        with redirect_stdout(f):
            print( "# vtk DataFile Version 3.0" )
            print( "vtk output" )
            print( "ASCII" )
            #print( "DATASET STRUCTURED_GRID" )
            print( "DATASET STRUCTURED_POINTS" )
            print( "DIMENSIONS ", gridXsize + 1 , " ", gridYsize + 1, " ", 1 )
            print( "ASPECT_RATIO ", searchRadius_h , " ", searchRadius_h , " ",  searchRadius_h )
            print( "ORIGIN ", gridXbegin , " ", gridYbegin , " ",  0  )
            print( "CELL_DATA ",  gridXsize * gridYsize * 1  )
            print( "SCALARS GridSector int 1 ")
            print( "LOOKUP_TABLE default" )
            for i in range( gridXsize * gridYsize * 1 ):
                print( gridSector[ i ] )

DomainGrid( gridXsize, gridYsize, 1,                    # grid size
            gridXbegin, gridYbegin, 0,                  # coordinates of grid origin
            np.zeros( gridXsize * gridYsize  ),         # array with index of grid sector
            'sources/openchannel_grid.vtk' )            # outputfile name
