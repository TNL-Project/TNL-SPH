#---------------------------------------------------------------------------#
#
# damBreak3D_WCSPH-DBC_benchmark
#
#---------------------------------------------------------------------------#

smoothingLentghCoef = 2
rho0 = 1000.
speedOfSound = 45.17
CFLnumber = 0.15

import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument( "-resolution", default=0.02, type=float )
#python3.9 and above: parser.add_argument( '--generateNewGeometry', default=False, action=argparse.BooleanOptionalAction )
parser.add_argument( '--generateNewGeometry', default=False, action='store_true' )
args = parser.parse_args()

dp =  args.resolution
generateNewGeometry = args.generateNewGeometry

print( f'Initial particle distance: {dp}.' )
print( f'Generate new files with geometry: {generateNewGeometry}.' )



if dp == 0.01:
    CFLnumber=0.15

if dp == 0.005:
    #CFLnumber = 0.095
    CFLnumber = 0.094
##---------------------------------------------------------------------------#

### Create related directories
import os
resultsPath = r'./results'
if not os.path.exists( resultsPath ):
    os.makedirs( resultsPath )

sourcesPath = r'./sources'
if not os.path.exists( sourcesPath ):
    os.makedirs( sourcesPath )

### Generate or load and process the geometry
if generateNewGeometry:
    import subprocess
    subprocess.check_call( [ './generateGeometryWithDualSPHysicsGenCase.sh', str( dp ) ], cwd='./template/generateGeometryWithDualSPHysicsGenCase/' )

import sys
sys.path.append('../../../tools')
import saveParticlesVTK
import numpy as np
import vtk

from vtk.numpy_interface import dataset_adapter as dsa

reader = vtk.vtkPolyDataReader()
reader.SetFileName( f'./template/generatedGeometries/dambreak_fluid_dp{dp}.vtk' )
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
saveParticlesVTK.save_polydata( fluidToWrite, "sources/dambreak_fluid.vtk" )

reader = vtk.vtkPolyDataReader()
reader.SetFileName( f'./template/generatedGeometries/dambreak_bound_dp{dp}.vtk' )
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
saveParticlesVTK.save_polydata( boundToWrite, "sources/dambreak_boundary.vtk" )

### Compute remaining parameters
spaceDimension = 3 #TODO: Move into templates
particleMass = rho0 * ( dp * dp * dp )
smoothingLentgh =  round( smoothingLentghCoef * dp, 7 )
searchRadius = round( smoothingLentgh * 2 , 7 )

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
gridSize_z = math.ceil( 1.5 * ( gridEnd_z - gridBegin_z ) / searchRadius )

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
fileParticleConf = fileParticleConf.replace( 'placeholderSearchRadius', str( searchRadius ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXSize', str( gridSize_x ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYSize', str( gridSize_y ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridZSize', str( gridSize_z ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXBegin', str( round( gridBegin_x, 8  ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYBegin', str( round( gridBegin_y, 8  ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridZBegin', str( round( gridBegin_y, 8  ) ) )

with open( 'sources/ParticlesConfig.h', 'w' ) as file:
  file.write( fileParticleConf )

# Setup of simulation control file
with open( 'template/SimulationControlConfig.h', 'r' ) as file :
  fileSimulationControl = file.read()

with open( 'sources/SimulationControlConfig.h', 'w' ) as file:
  file.write( fileSimulationControl )
