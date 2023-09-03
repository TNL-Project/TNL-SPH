#---------------------------------------------------------------------------#
#
# damBreak3D_WCSPH-DBC_benchmark
#
#---------------------------------------------------------------------------#
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument( "-resolution", default=0.02, type=float )
args = parser.parse_args()
#---------------------------------------------------------------------------#

dp =  args.resolution
print( f'Initial particle distance: {dp}.' )

spaceDimension = 3
smoothingLentghCoef = 2
rho0 = 1000.
speedOfSound = 45.17
CFLnumber = 0.15

# Setup number of subdomains (devices), to split the problem [-]:
numberOfSubdomains = 3

#---------------------------------------------------------------------------#

### Compute remaining parameters
particleMass = rho0 * ( dp * dp * dp )
smoothingLentgh =  round( smoothingLentghCoef * dp, 7 )
searchRadius = round( smoothingLentgh * 2 , 7 )
timeStep = round( CFLnumber * ( smoothingLentgh / speedOfSound ), 8 )
coefB = round( speedOfSound * speedOfSound * rho0 / 7 , 1 )

### Create necessary directories for problem settings and to store the results
import os
resultsPath = r'./results'
if not os.path.exists( resultsPath ):
    os.makedirs( resultsPath )

sourcesPath = r'./sources'
if not os.path.exists( sourcesPath ):
    os.makedirs( sourcesPath )

### Process external geometries - load particles and fields to compute domain dimensions and create subdomains
import sys
sys.path.append('../../tools')

from generateCaseTool import processExternalParticleVTKFile
np_points_fluid = processExternalParticleVTKFile( 'template/geometry/fluid_dp0-02.vtk', "sources/dambreak_fluid.vtk", [] )
np_points_box = processExternalParticleVTKFile( 'template/geometry/boundary_dp0-02.vtk', "sources/dambreak_boundary.vtk", [] )

### Compute domain limits
from generateCaseTool import getDomainLimits
gridBegin, gridEnd, gridSize = getDomainLimits( np_points_fluid, np_points_box, searchRadius, 1., 1., 1.4 )

### Split the particle data into subdomains
from splitToSubdomains import *

subdomains = Subdomains( numberOfSubdomains )
subdomains.divideIntoSubdomains( searchRadius, np_points_fluid, gridSize, gridBegin )
subdomainConfigString, inputFluidFilesString, inputBoundaryFilesString = subdomains.generateSubdomains( np_points_fluid, [], np_points_box, [] )

#---------------------------------------------------------------------------#

# Read in the file
with open( 'template/SPHCaseConfig_template.h', 'r' ) as file :
  fileSPHConf = file.read()

# Replace the target string
fileSPHConf = fileSPHConf.replace( 'placeholderDimension', str( spaceDimension ) )
fileSPHConf = fileSPHConf.replace( 'placeholderMass', str( round( particleMass, 8 ) ) )
fileSPHConf = fileSPHConf.replace( 'placeholderSpeedOfSound', str( speedOfSound ) )
fileSPHConf = fileSPHConf.replace( 'placeholderCoefB', str( coefB ) )
fileSPHConf = fileSPHConf.replace( 'placeholderDensity', str( rho0 ))
fileSPHConf = fileSPHConf.replace( 'placeholderInitParticleDistance', str( dp ) )
fileSPHConf = fileSPHConf.replace( 'placeholderSmoothingLength', str( smoothingLentgh ) )
fileSPHConf = fileSPHConf.replace( 'placeholderTimeStep', str( timeStep ) )

# Write the file out again
with open( 'sources/SPHCaseConfig.h', 'w' ) as file:
  file.write( fileSPHConf )

#---------------------------------------------------------------------------#

# Read and write (with possible edit) simulation control file.
with open( 'template/SimulationControlConfig_template.h', 'r' ) as file :
  fileSimulationControl = file.read()
  fileSimulationControl = fileSimulationControl.replace( '#placeholderNumberOfSubdomains', str( numberOfSubdomains ) )
  fileSimulationControl = fileSimulationControl.replace( '#placeholderInputFluidFiles', inputFluidFilesString )
  fileSimulationControl = fileSimulationControl.replace( '#placeholderInputBoundaryFiles', inputBoundaryFilesString )

with open( 'sources/SimulationControlConfig.h', 'w' ) as file:
  file.write( fileSimulationControl )

#---------------------------------------------------------------------------#

with open( 'template/ParticlesConfig_template.h', 'r' ) as file :
  fileParticleConf = file.read()
  fileParticleConf = fileParticleConf.replace( '#placeholderSubdomainInfo', subdomainConfigString )
  fileParticleConf = fileParticleConf.replace( '#placeholderNumberOfSubdomains', str( numberOfSubdomains ) )

# Write the file out again
with open( 'sources/ParticlesConfig.h', 'w' ) as file:
  file.write( fileParticleConf )

#---------------------------------------------------------------------------#
