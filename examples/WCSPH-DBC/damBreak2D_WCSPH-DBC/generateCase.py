#---------------------------------------------------------------------------#
#
# case: damBreak2D_RSPH
#
#---------------------------------------------------------------------------#
### Parameters of the case necessary for case creation:

# Dimensions of box:
boxL = 1.61
boxH = 0.8

# Dimensions of the fluid box
fluidL = 0.6
fluidH = 0.3

# Initial particle distance (dp):
dp = 0.002

# Smoothing length coefitient:
# smoothing length (h) = smoothing length coef (Coef_h) * initial particle distance (d_p)
# [ h = Coef_h * dp ]
smoothingLentghCoef = 2**0.5

# Referential density of the medium (rho0):
rho0 = 1000.
p0 = 0. #DELETE

# Number of boundary layers:
numberOfBoundaryLayers = 3

# Numerical speed of sound (c0):
speedOfSound = 34.3

# Initial time step.
# In case, that initial time step is not defined, is computed automatically.
timeStep = 0.00002

# CFL number (CFL):
CFLnumber = 0.2
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

fluidL_n = round( fluidL / dp )
fluidH_n = round( fluidH / dp )

for x in range( fluidL_n ):
    for y in range( fluidH_n ):
        fluid_rx.append( dp * ( x + 1 ) )
        fluid_ry.append( dp * ( y + 1 ) )
        fluid_rz.append( 0. )

### Generate boundary particles
box_rx = []; box_ry = []; box_rz = []
ghost_rx = []; ghost_ry = []; ghost_rz = []

boxL_n = round( boxL / dp )
boxH_n = round( boxH / dp )

# left wall
for layer in range( numberOfBoundaryLayers ):
    for y in range( boxH_n - 1 ):
        box_rx.append( 0. - layer * dp )
        box_ry.append( ( y+1 ) * dp )
        box_rz.append( 0. )

        ghost_rx.append( 0. + dp * ( layer + 1 ) )
        ghost_ry.append( ( y + 1 ) * dp )
        ghost_rz.append( 0.)

# bottom wall
for layer in range( numberOfBoundaryLayers ):
    for x in range( boxL_n - numberOfBoundaryLayers + 1 ):
        box_rx.append( ( x + 1 ) * dp )
        box_ry.append( 0. - layer * dp )
        box_rz.append( 0. ) #we use only 2D case

        ghost_rx.append( ( x + 1 ) * dp )
        ghost_ry.append( 0. + dp * ( layer + 1 ) )
        ghost_rz.append( 0.)

x_last = box_rx[ -1 ] + dp #due to discretisation, we need to save last value of bottom wall

# right wall
for layer in range( numberOfBoundaryLayers ):
    for y in range( boxH_n - 1 ):
        box_rx.append( x_last + dp * layer )
        box_ry.append( ( y + 1 ) * dp )
        box_rz.append( 0. ) #we use only 2D case

        ghost_rx.append( x_last - dp * ( layer + 1 ) )
        ghost_ry.append( ( y + 1 ) * dp )
        ghost_rz.append( 0. )

# generate the corners
def generate90degCorner( x, y, dirx, diry ):
  for layer in range( numberOfBoundaryLayers ):
    for k in range( numberOfBoundaryLayers ):
      box_rx.append( x + k * dp * dirx )
      box_ry.append( y + layer * dp * diry )
      box_rz.append( 0. )

      ghost_rx.append( x + ( k + 1 ) * dp * dirx * ( -1 ) )
      ghost_ry.append( y + ( layer + 1 ) * dp * diry * ( -1 ) )
      ghost_rz.append( 0. )

generate90degCorner( 0, 0., -1, -1 )
generate90degCorner( x_last, 0., +1, -1 )

import sys
sys.path.append('../../../src/tools')
import saveParticlesVTK
import numpy as np

fluid_r = np.array( ( fluid_rx, fluid_ry, fluid_rz ), dtype=float ).T #!!
fluid_v = np.zeros( ( len( fluid_rx ), 3 ) )
fluid_rho = rho0 * np.ones( len( fluid_rx ) )
fluid_p = np.zeros( len( fluid_rx ) )
fluid_ptype = np.zeros( len( fluid_rx ) )

fluidToWrite = saveParticlesVTK.create_pointcloud_polydata( fluid_r, fluid_v, fluid_rho, fluid_p, fluid_ptype )
saveParticlesVTK.save_polydata( fluidToWrite, "sources/dambreak_fluid.vtk" )

boundary_r = np.array( ( box_rx, box_ry, box_rz ), dtype=float ).T #!!
boundary_ghostNodes = np.array( ( ghost_rx, ghost_ry, ghost_rz ), dtype=float ).T #!!
boundary_v = np.zeros( ( len( box_rx ), 3 ) )
boundary_rho = rho0 * np.ones( len( box_rx ) )
boundary_p = np.zeros( len( box_rx ) )
boundary_ptype = np.ones( len( box_rx ) )

boxToWrite = saveParticlesVTK.create_pointcloud_polydata( boundary_r, boundary_v, boundary_rho, boundary_p, boundary_ptype,
                                                          ghostNodes=boundary_ghostNodes )
saveParticlesVTK.save_polydata( boxToWrite, "sources/dambreak_boundary.vtk" )

### Compute remaining parameters
spaceDimension = 2 #TODO: Move into templates.
particleMass = rho0 * ( dp * dp )
smoothingLentgh =  round( smoothingLentghCoef * dp, 7 )
searchRadius = round( smoothingLentgh * 2 , 7 )

if not timeStep:
    timeStep = round( CFLnumber * ( smoothingLentgh / speedOfSound ), 8 )
coefB = round( speedOfSound * speedOfSound * rho0 / 7 , 1 )

### Compute remaining domain parameters
from math import ceil
gridBegin_x = 1.005 * ( min( min( fluid_rx ), min( box_rx ) ) - searchRadius )
gridBegin_y = 1.005 * ( min( min( fluid_rz ), min( box_ry ) ) - searchRadius )
gridEnd_x = 1.005 * ( max( max( fluid_rx ), max( box_rx ) ) + searchRadius )
gridEnd_y = 1.2 * ( max( max( fluid_rz ), max( box_ry ) ) + searchRadius )

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

# Setup of particle system
with open( 'template/ParticlesConfig_template.h', 'r' ) as file :
  fileParticleConf = file.read()

fileParticleConf = fileParticleConf.replace( 'placeholderDimension', str( spaceDimension ) )
fileParticleConf = fileParticleConf.replace( 'placeholderFluidParticles', str( len( fluid_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedFluidParticles', str( len( fluid_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderBoundaryParticles', str( len( box_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedBoundaryParticles', str( len( box_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderSearchRadius', str( searchRadius ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXSize', str( gridSize_x ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYSize', str( gridSize_y ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXBegin', str( round( gridBegin_x, 9  ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYBegin', str( round( gridBegin_y, 9  ) ) )

with open( 'sources/ParticlesConfig.h', 'w' ) as file:
  file.write( fileParticleConf )

# Setup of simulation control file
with open( 'template/SimulationControlConfig.h', 'r' ) as file :
  fileSimulationControl = file.read()

with open( 'sources/SimulationControlConfig.h', 'w' ) as file:
  file.write( fileSimulationControl )

# Setup the measuretool config
with open( 'template/MeasuretoolConfig.h', 'r' ) as file :
  fileMeasuretoolConf = file.read()

fileMeasuretoolConf = fileMeasuretoolConf.replace( 'placeholderInitParticleDistance', str( dp ) )
fileMeasuretoolConf = fileMeasuretoolConf.replace( 'placeholderSmoothingLength', str( smoothingLentgh ) )

with open( 'sources/MeasuretoolConfig.h', 'w' ) as file:
  file.write( fileMeasuretoolConf )
