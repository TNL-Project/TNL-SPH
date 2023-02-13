# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#
# damBreak3D_WCSPH-DBC_benchmark
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument( "-resolution" )
args = parser.parse_args()
print("Input resolution: .....", args.resolution)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

boxL = 1.61
boxH = 0.8

fluidL = 0.6
fluidH = 0.3

dp = float( args.resolution )
#dp = 0.02
smoothingLentghCoef = 2

rho0 = 1000.
p0 = 0.

numberOfBoundaryLayers = 3

#speedOfSound = 45.17167357703276
speedOfSound = 45.17
CFLnumber = 0.15

if dp == 0.01:
    CFLnumber=0.15

if dp == 0.005:
    #CFLnumber = 0.095
    CFLnumber = 0.094

#timeStep = 0.00002 #otherwise is obtained automatically

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

import sys
sys.path.append('../../../tools')
import saveParticlesVTK
import numpy as np
import vtk

from vtk.numpy_interface import dataset_adapter as dsa

reader = vtk.vtkPolyDataReader()
reader.SetFileName( './template/damBreak3D_WCSPH-DBC_out/damBreak3D_WCSPH-DBC_Fluid.vtk' )
reader.ReadAllScalarsOn()
reader.ReadAllVectorsOn()
reader.Update()

#It is possible to access PointData, CellData, FieldData, Points (subclasses of vtkPointSet only), Polygons (vtkPolyData only) this way.
polydata = reader.GetOutput()
np_points_fluid = dsa.WrapDataObject( polydata ).Points

r = np.array( np_points_fluid, dtype=float ) #!!
v = np.array( dsa.WrapDataObject( polydata ).PointData[ 'Vel' ], dtype=float )
rho = np.array( dsa.WrapDataObject( polydata ).PointData[ 'Rhop' ] )
p = np.zeros( len( np_points_fluid ) )
ptype = np.zeros( len( np_points_fluid ) )

fluidToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
saveParticlesVTK.save_polydata( fluidToWrite, "dambreak_fluid.vtk" )

reader = vtk.vtkPolyDataReader()
reader.SetFileName( './template/damBreak3D_WCSPH-DBC_out/damBreak3D_WCSPH-DBC_Bound.vtk' )
reader.ReadAllScalarsOn()
reader.ReadAllVectorsOn()
reader.Update()

#It is possible to access PointData, CellData, FieldData, Points (subclasses of vtkPointSet only), Polygons (vtkPolyData only) this way.
polydata = reader.GetOutput()
np_points_box = dsa.WrapDataObject( polydata ).Points
#np_data = dsa.WrapDataObject( polydata ).PointData[ 'Vel.m_average' ]

r = np.array( np_points_box, dtype=float ) #!!
v = np.zeros( ( len( np_points_box ), 3 ) )
rho = rho0 * np.ones( len( np_points_box ) )
p = np.zeros( len( np_points_box ) )
ptype = np.zeros( len( np_points_box ) )

fluidToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
saveParticlesVTK.save_polydata( fluidToWrite, "dambreak_bound.vtk" )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
### Compute remaining parameters
particleMass = rho0 * ( dp * dp * dp )
smoothingLentgh =  round( smoothingLentghCoef * dp, 7 )
searchRadius = round( smoothingLentgh * 2 , 7 )
#if not timeStep:
timeStep = round( CFLnumber * ( smoothingLentgh / speedOfSound ), 8 )
coefB = round( speedOfSound * speedOfSound * rho0 / 7 , 1 )
spaceDimension = 3

#Determine grid size
import math
gridXbegin = 1.005 * ( min( min( np_points_fluid[ : , 0 ] ), min( np_points_box[ :, 0 ] ) ) ) - searchRadius
gridYbegin = 1.005 * ( min( min( np_points_fluid[ : , 1 ] ), min( np_points_box[ :, 1 ] ) ) ) - searchRadius
gridZbegin = 1.005 * ( min( min( np_points_fluid[ : , 2 ] ), min( np_points_box[ : ,2 ] ) ) ) - searchRadius

gridXend = 1.005 * ( max( max( np_points_fluid[ :, 0 ] ), max( np_points_box[ :, 0 ] ) ) ) + searchRadius
gridYend = 1.005 * ( max( max( np_points_fluid[ :, 1 ] ), max( np_points_box[ :, 1 ] ) ) ) + searchRadius
gridZend = 1.005 * ( max( max( np_points_fluid[ :, 2 ] ), max( np_points_box[ :, 2 ] ) ) ) + searchRadius

gridXsize = math.ceil( ( gridXend - gridXbegin ) / searchRadius )
gridYsize = math.ceil( ( gridYend - gridYbegin ) / searchRadius )
gridZsize = math.ceil( ( gridZend - gridZbegin ) / searchRadius ) + 5
if dp == 0.005:
    gridZsize += 30

# Read in the file
with open( 'template/SPHCaseConfig_template.h', 'r' ) as file :
  fileSPHConf = file.read()

# Replace the target string
fileSPHConf = fileSPHConf.replace( 'placeholderDimension', str( spaceDimension ) )
fileSPHConf = fileSPHConf.replace( 'placeholderMass', str( round( particleMass, 8 ) ) )
fileSPHConf = fileSPHConf.replace( 'placeholderSpeedOfSound', str( speedOfSound ) )
fileSPHConf = fileSPHConf.replace( 'placeholderCoefB', str( coefB ) )
fileSPHConf = fileSPHConf.replace( 'placeholderDensity', str( rho0 ))
fileSPHConf = fileSPHConf.replace( 'placeholderSmoothingLength', str( smoothingLentgh ) )
fileSPHConf = fileSPHConf.replace( 'placeholderTimeStep', str( timeStep ) )

# Write the file out again
with open( 'SPHCaseConfig.h', 'w' ) as file:
  file.write( fileSPHConf )

with open( 'template/ParticlesConfig_template.h', 'r' ) as file :
  fileParticleConf = file.read()

# Replace the target string
fileParticleConf = fileParticleConf.replace( 'placeholderDimension', str( spaceDimension ) )
fileParticleConf = fileParticleConf.replace( 'placeholderFluidParticles', str( len( np_points_fluid ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedFluidParticles', str( len( np_points_fluid ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderBoundaryParticles', str( len( np_points_box ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedBoundaryParticles', str( len( np_points_box ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderSearchRadius', str( searchRadius ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXSize', str( gridXsize ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYSize', str( gridYsize ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridZSize', str( gridZsize ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXBegin', str( round( gridXbegin, 8  ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYBegin', str( round( gridYbegin, 8  ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridZBegin', str( round( gridYbegin, 8  ) ) )

# Write the file out again
with open( 'ParticlesConfig.h', 'w' ) as file:
  file.write( fileParticleConf )
