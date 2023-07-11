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

#dp = float( args.resolution )
dp = 0.02
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
#sys.path.append('../../../tools')
sys.path.append('../../tools')
import saveParticlesVTK
import numpy as np
import vtk

from vtk.numpy_interface import dataset_adapter as dsa

reader = vtk.vtkPolyDataReader()
#reader.SetFileName( './template/damBreak3D_WCSPH-DBC_out/damBreak3D_WCSPH-DBC_Fluid.vtk' )
reader.SetFileName( './fluid_new.vtk' )
#reader.SetFileName( './template/damBreak3D_WCSPH-DBC_out/damBreak3D_WCSPH-DBC_Bound.vtk' )
reader.ReadAllScalarsOn()
reader.ReadAllVectorsOn()
reader.Update()

#It is possible to access PointData, CellData, FieldData, Points (subclasses of vtkPointSet only), Polygons (vtkPolyData only) this way.
polydata = reader.GetOutput()
np_points_fluid = dsa.WrapDataObject( polydata ).Points

r = np.array( np_points_fluid, dtype=float ) #!!
v = np.array( dsa.WrapDataObject( polydata ).PointData[ 'Vel' ], dtype=float )
#v = np.array( ( np_points_fluid[ :, 0 ], np_points_fluid[ :, 1 ], np_points_fluid[ :, 2 ], np.zeros(  len( np_points_fluid ) ) ), dtype=float )
rho = np.array( dsa.WrapDataObject( polydata ).PointData[ 'Rhop' ] )
p = np.zeros( len( np_points_fluid ) )
ptype = np.zeros( len( np_points_fluid ) )

fluidToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
saveParticlesVTK.save_polydata( fluidToWrite, "dambreak_fluid.vtk" )

"""
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
#r = np.array( ( np_points_box[ :, 0 ], np_points_box[ :, 1 ], np_points_box[ :, 2 ], np.zeros(  len( np_points_box ) ) ), dtype=float ) #!!
#v = np.zeros( ( len( np_points_box ), 4 ) )
v = np.array( ( np_points_box[ :, 0 ], np_points_box[ :, 1 ], np_points_box[ :, 2 ], np.zeros(  len( np_points_box ) ) ), dtype=float )
rho = rho0 * np.ones( len( np_points_box ) )
p = np.zeros( len( np_points_box ) )
ptype = np.zeros( len( np_points_box ) )

#fluidToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
#saveParticlesVTK.save_polydata( fluidToWrite, "dambreak_bound.vtk" )
"""

reader = vtk.vtkPolyDataReader()
#reader.SetFileName( './template/damBreak3D_WCSPH-DBC_out/damBreak3D_WCSPH-DBC_Fluid.vtk' )
#reader.SetFileName( './template/damBreak3D_WCSPH-DBC_out/damBreak3D_WCSPH-DBC_Bound.vtk' )
reader.SetFileName( './boundary_new.vtk' )
#reader.SetFileName( './template/damBreak3D_WCSPH-DBC_out/boundary_new.vtk' )
reader.ReadAllScalarsOn()
reader.ReadAllVectorsOn()
reader.Update()

#It is possible to access PointData, CellData, FieldData, Points (subclasses of vtkPointSet only), Polygons (vtkPolyData only) this way.
polydata = reader.GetOutput()
np_points_box = dsa.WrapDataObject( polydata ).Points

r = np.array( np_points_box, dtype=float ) #!!
v = np.array( dsa.WrapDataObject( polydata ).PointData[ 'Vel' ], dtype=float )
#v = np.array( ( np_points_bound[ :, 0 ], np_points_fluid[ :, 1 ], np_points_fluid[ :, 2 ], np.zeros(  len( np_points_fluid ) ) ), dtype=float )
rho = np.array( dsa.WrapDataObject( polydata ).PointData[ 'Rhop' ] )
p = np.zeros( len( np_points_box ) )
ptype = np.zeros( len( np_points_box ) )

print("Number of boundary particles to write: ", len(r) )
boundToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
saveParticlesVTK.save_polydata( boundToWrite, "dambreak_bound.vtk" )

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
gridXbegin = 1.005 * ( min( min( np_points_fluid[ : , 0 ] ), min( np_points_box[ :, 0 ] ) ) - searchRadius )
gridYbegin = 1.005 * ( min( min( np_points_fluid[ : , 1 ] ), min( np_points_box[ :, 1 ] ) ) - searchRadius )
gridZbegin = 1.005 * ( min( min( np_points_fluid[ : , 2 ] ), min( np_points_box[ : ,2 ] ) ) - searchRadius )

gridXend = 1.005 * ( max( max( np_points_fluid[ :, 0 ] ), max( np_points_box[ :, 0 ] ) ) + searchRadius )
gridYend = 1.005 * ( max( max( np_points_fluid[ :, 1 ] ), max( np_points_box[ :, 1 ] ) ) + searchRadius )
gridZend = 1.005 * ( max( max( np_points_fluid[ :, 2 ] ), max( np_points_box[ :, 2 ] ) ) + searchRadius )

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
fileSPHConf = fileSPHConf.replace( 'placeholderInitParticleDistance', str( dp ) )
fileSPHConf = fileSPHConf.replace( 'placeholderSmoothingLength', str( smoothingLentgh ) )
fileSPHConf = fileSPHConf.replace( 'placeholderTimeStep', str( timeStep ) )

# Write the file out again
with open( 'SPHCaseConfig.h', 'w' ) as file:
  file.write( fileSPHConf )

#with open( 'template/ParticlesConfig_template.h', 'r' ) as file :
#  fileParticleConf = file.read()
#
## Replace the target string
#fileParticleConf = fileParticleConf.replace( 'placeholderDimension', str( spaceDimension ) )
#fileParticleConf = fileParticleConf.replace( 'placeholderFluidParticles', str( len( np_points_fluid ) ) )
#fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedFluidParticles', str( len( np_points_fluid ) ) )
#fileParticleConf = fileParticleConf.replace( 'placeholderBoundaryParticles', str( len( np_points_box ) ) )
#fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedBoundaryParticles', str( len( np_points_box ) ) )
#fileParticleConf = fileParticleConf.replace( 'placeholderSearchRadius', str( searchRadius ) )
#fileParticleConf = fileParticleConf.replace( 'placeholderGridXSize', str( gridXsize ) )
#fileParticleConf = fileParticleConf.replace( 'placeholderGridYSize', str( gridYsize ) )
#fileParticleConf = fileParticleConf.replace( 'placeholderGridZSize', str( gridZsize ) )
#fileParticleConf = fileParticleConf.replace( 'placeholderGridXBegin', str( round( gridXbegin, 8  ) ) )
#fileParticleConf = fileParticleConf.replace( 'placeholderGridYBegin', str( round( gridYbegin, 8  ) ) )
#fileParticleConf = fileParticleConf.replace( 'placeholderGridZBegin', str( round( gridYbegin, 8  ) ) )
#
## Write the file out again
#with open( 'ParticlesConfig.h', 'w' ) as file:
#  file.write( fileParticleConf )

# Read and write (with possible edit) simulation control file.
with open( 'template/SimulationControlConfig_template.h', 'r' ) as file :
  fileSimulationControl = file.read()

with open( 'SimulationControlConfig.h', 'w' ) as file:
  file.write( fileSimulationControl )

import os
resultsPath = r'./results'
if not os.path.exists( resultsPath ):
    os.makedirs( resultsPath )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# split
#tot: numberOfPtcsTotal = len( fluid_rx ) + len( box_rx )
# split based on fluid:
fluid_rx = np_points_fluid[ :, 0 ]
fluid_ry = np_points_fluid[ :, 1 ]
fluid_rz = np_points_fluid[ :, 2 ]

box_rx = np_points_box[ :, 0 ]
box_ry = np_points_box[ :, 1 ]
box_rz = np_points_box[ :, 2 ]

numberOfPtcsTotal = len( fluid_rx )

middleParticle = ( int )( numberOfPtcsTotal / 2 )
print(" midleParticleNumber: ", middleParticle )

searchRadius_h = round( smoothingLentgh * 2 , 7 )
gridXsplit = math.ceil( fluid_rx[ middleParticle ] / searchRadius_h )
print(" searchRadius: ", gridXsplit )

print( "Grid-global size: [ ", gridXsize , ",", gridYsize, ",", gridZsize, " ]." )
print( "Grid-global begin: [ ", gridXbegin , ",", gridYbegin, ",", gridZsize, " ]." )

gridXsplitBegin = gridXbegin + gridXsplit * searchRadius
print( "Grid-G1 size: [ ", gridXsplit , ",", gridYsize, ",", gridZsize, " ]." )
print( "Grid-G1 begin: [ ", gridXbegin , ",", gridYbegin, ",", gridZbegin, " ]." )

gridXsplitBegin = gridXbegin + gridXsplit * searchRadius
print( "Grid-G2 size: [ ", gridXsize - gridXsplit , ",", gridYsize, ",", gridZsize, " ]." )
print( "Grid-G2 begin: [ ", gridXsplitBegin , ",", gridYbegin, ",", gridZbegin, " ]." )

#=====================================================================================================

fluid_rx_g1 = []
fluid_ry_g1 = []
fluid_rz_g1 = []
ptype_fluid_g1 = []

box_rx_g1 = []
box_ry_g1 = []
box_rz_g1 = []
ptype_box_g1 = []

counter_g1_nptcs = 0
counter_g1_nptcsReal = 0

for i in range ( len( fluid_rx ) ):
    if( fluid_rx[ i ] <= gridXsplitBegin ):
        fluid_rx_g1.append( fluid_rx[ i ] )
        fluid_ry_g1.append( fluid_ry[ i ] )
        fluid_rz_g1.append( fluid_rz[ i ] )
        ptype_fluid_g1.append( 0 )
        counter_g1_nptcs += 1
        counter_g1_nptcsReal += 1
    if( fluid_rx[ i ] > gridXsplitBegin and fluid_rx[ i ] <= gridXsplitBegin + searchRadius_h ):
        fluid_rx_g1.append( fluid_rx[ i ] )
        fluid_ry_g1.append( fluid_ry[ i ] )
        fluid_rz_g1.append( fluid_rz[ i ] )
        ptype_fluid_g1.append( 1 )
        counter_g1_nptcsReal += 1

for i in range ( len( box_rx ) ):
    if( box_rx[ i ] <= gridXsplitBegin ):
        box_rx_g1.append( box_rx[ i ] )
        box_ry_g1.append( box_ry[ i ] )
        box_rz_g1.append( box_rz[ i ] )
        ptype_box_g1.append( 0 )
    if( box_rx[ i ] > gridXsplitBegin and box_rx[ i ] <= gridXsplitBegin + searchRadius_h ):
        box_rx_g1.append( box_rx[ i ] )
        box_ry_g1.append( box_ry[ i ] )
        box_rz_g1.append( box_rz[ i ] )
        ptype_box_g1.append( 1 )

rg1 = np.array( ( fluid_rx_g1, fluid_ry_g1, fluid_rz_g1 ), dtype=float ).T #!!
vg1 = np.zeros( ( len( fluid_rx_g1 ), 3 ) )
rhog1 = rho0 * np.ones( len( fluid_rx_g1 ) )
pg1 = np.zeros( len( fluid_rx_g1 ) )
#ptypeg1 = np.zeros( len( fluid_rx_g1 ) )
ptypeg1 = np.array( ( ptype_fluid_g1 ), dtype=float ).T

fluidToWrite = saveParticlesVTK.create_pointcloud_polydata( rg1, vg1, rhog1, pg1, ptypeg1 )
saveParticlesVTK.save_polydata( fluidToWrite, "dambreak_fluid_g1.vtk" )

rg1 = np.array( ( box_rx_g1, box_ry_g1, box_rz_g1 ), dtype=float ).T #!!
vg1 = np.zeros( ( len( box_rx_g1 ), 3 ) )
rhog1 = rho0 * np.ones( len( box_rx_g1 ) )
pg1 = np.zeros( len( box_rx_g1 ) )
#ptypeg1 = np.ones( len( box_rx_g1 ) )
ptypeg1 = np.array( ( ptype_box_g1 ), dtype=float ).T

boxToWrite = saveParticlesVTK.create_pointcloud_polydata( rg1, vg1, rhog1, pg1, ptypeg1 )
saveParticlesVTK.save_polydata( boxToWrite, "dambreak_boundary_g1.vtk" )

#=====================================================================================================

fluid_rx_g2 = []
fluid_ry_g2 = []
fluid_rz_g2 = []
ptype_fluid_g2 = []
counter_g2_nptcs = 0
counter_g2_nptcsReal = 0

box_rx_g2 = []
box_ry_g2 = []
box_rz_g2 = []
ptype_box_g2 = []

for i in range ( len( fluid_rx ) ):
    if( fluid_rx[ i ] >= gridXsplitBegin ):
        fluid_rx_g2.append( fluid_rx[ i ] )
        fluid_ry_g2.append( fluid_ry[ i ] )
        fluid_rz_g2.append( fluid_rz[ i ] )
        ptype_fluid_g2.append( 2 )
        counter_g2_nptcs += 1
        counter_g2_nptcsReal += 1
    if( fluid_rx[ i ] > gridXsplitBegin - searchRadius_h and fluid_rx[ i ] < gridXsplitBegin ):
        fluid_rx_g2.append( fluid_rx[ i ] )
        fluid_ry_g2.append( fluid_ry[ i ] )
        fluid_rz_g2.append( fluid_rz[ i ] )
        ptype_fluid_g2.append( 3 )
        counter_g2_nptcsReal += 1

for i in range ( len( box_rx ) ):
    if( box_rx[ i ] >= gridXsplitBegin ):
        box_rx_g2.append( box_rx[ i ] )
        box_ry_g2.append( box_ry[ i ] )
        box_rz_g2.append( box_rz[ i ] )
        ptype_box_g2.append( 2 )
    if( box_rx[ i ] > gridXsplitBegin - searchRadius_h and box_rx[ i ] < gridXsplitBegin ):
        box_rx_g2.append( box_rx[ i ] )
        box_ry_g2.append( box_ry[ i ] )
        box_rz_g2.append( box_rz[ i ] )
        ptype_box_g2.append( 3 )

rg2 = np.array( ( fluid_rx_g2, fluid_ry_g2, fluid_rz_g2 ), dtype=float ).T #!!
vg2 = np.zeros( ( len( fluid_rx_g2 ), 3 ) )
rhog2 = rho0 * np.ones( len( fluid_rx_g2 ) )
pg2 = np.zeros( len( fluid_rx_g2 ) )
#ptypeg2 = np.zeros( len( fluid_rx_g2 ) )
ptypeg2 = np.array( ( ptype_fluid_g2 ), dtype=float ).T

fluidToWrite = saveParticlesVTK.create_pointcloud_polydata( rg2, vg2, rhog2, pg2, ptypeg2 )
saveParticlesVTK.save_polydata( fluidToWrite, "dambreak_fluid_g2.vtk" )

rg2 = np.array( ( box_rx_g2, box_ry_g2, box_rz_g2 ), dtype=float ).T #!!
vg2 = np.zeros( ( len( box_rx_g2 ), 3 ) )
rhog2 = rho0 * np.ones( len( box_rx_g2 ) )
pg2 = np.zeros( len( box_rx_g2 ) )
#ptypeg2 = np.ones( len( box_rx_g2 ) )
ptypeg2 = np.array( ( ptype_box_g2 ), dtype=float ).T

boxToWrite = saveParticlesVTK.create_pointcloud_polydata( rg2, vg2, rhog2, pg2, ptypeg2 )
saveParticlesVTK.save_polydata( boxToWrite, "dambreak_boundary_g2.vtk" )

#---------------------------------------------------------------------------------------------------

#====================================================================================================

gridCoordinates_x = []
gridCoordinates_y = []
gridCoordinates_z = []

gridSector = []

for z in range ( gridZsize ):
    for y in range ( gridYsize ):
        for x in range ( gridXsize ):
            gridCoordinates_x.append( x * searchRadius_h )
            gridCoordinates_y.append( y * searchRadius_h )
            #gridCoordinates_z.append( 0. )
            gridCoordinates_z.append( z * searchRadius_h )

            if x < gridXsplit:
                gridSector.append( 0 )
            elif x == gridXsplit:
                gridSector.append( 1 )
            elif x == gridXsplit + 1:
                gridSector.append( 2 )
            elif x > gridXsplit + 1:
                gridSector.append( 3 )
            else:
                printf(" Invalid grid coordinates! ")



from contextlib import redirect_stdout

def DomainGrid( gridXsize, gridYsize, gridZsize, gridXbegin, gridYbegin, gridZbegin, gridSector, name ):
    #with open( "distributedGrid.vtk", 'w' ) as f:
    with open( name, 'w' ) as f:
        with redirect_stdout(f):
            print( "# vtk DataFile Version 3.0" )
            print( "vtk output" )
            print( "ASCII" )
            #print( "DATASET STRUCTURED_GRID" )
            print( "DATASET STRUCTURED_POINTS" )
            print( "DIMENSIONS ", gridXsize + 1 , " ", gridYsize + 1, " ", gridZsize + 1 )
            #print( "POINTS ", gridXsize * gridYsize * 1, " float" )
            #for i in range( gridXsize * gridYsize * 1 ):
            #    print( "", gridCoordinates_x[ i ], " ", gridCoordinates_y[ i ], " ", gridCoordinates_z[ i ] )
            #print( "ASPECT_RATIO 1 1 1" )
            print( "ASPECT_RATIO ", searchRadius_h , " ", searchRadius_h , " ",  searchRadius_h )
            print( "ORIGIN ", gridXbegin , " ", gridYbegin , " ",  gridZbegin  )
            ##print( "POINT_DATA ",  gridXsize * gridYsize * 1  )
            print( "CELL_DATA ",  gridXsize * gridYsize * gridZsize  )
            #print( "FIELD FieldData 1" )
            #print( "GridSector 1",  gridXsize * gridYsize * 1 , " int" )
            print( "SCALARS GridSector int 1 ")
            print( "LOOKUP_TABLE default" )
            for i in range( gridXsize * gridYsize * gridZsize ):
                print( gridSector[ i ] )
            #print( "VECTORS GridIndex int 3 ")
            #print( "LOOKUP_TABLE default" )
            #for i in range( gridXsize * gridYsize * 1 ):
            #    print( gridSector[ i ], " ", gridSector[ i ], " ", gridSector[ i ] )
    print("Done.")

print( f'Grid sector size: {len(gridSector)}, grid size: { gridXsize * gridYsize * gridZsize }' )
# Write global grid as example
DomainGrid( gridXsize, gridYsize, gridZsize,    # grid size
            gridXbegin, gridYbegin, gridZbegin,  # coordinates of rgrid origina
            gridSector,                 # array with index of grid sector
            "distributedGrid.vtk" )     # outputfile name

gridSector_g1 = []

for z in range ( gridZsize ):
    for y in range ( gridYsize ):
        for x in range ( gridXsplit + 1 ):

            if x < gridXsplit:
                gridSector_g1.append( 0 )
            elif x == gridXsplit:
                gridSector_g1.append( 1 )
            elif x == gridXsplit + 1:
                gridSector_g1.append( 2 )
            elif x > gridXsplit + 1:
                gridSector_g1.append( 3 )
            else:
                printf(" Invalid grid coordinates! ")

# Write local grid G1
DomainGrid( gridXsplit + 1, gridYsize, gridZsize,   # grid size
            gridXbegin, gridYbegin, gridZbegin,      # coordinates of rgrid origina
            gridSector_g1,                  # array with index of grid sector
            "distributedGrid_g1.vtk" )      # outputfile name

gridSector_g2 = []

grid2Size = ( gridXsize - gridXsplit + 1 ) * gridYsize * 1

for z in range ( gridZsize ):
    for y in range ( gridYsize ):
        for x in range ( gridXsplit - 1 , gridXsize ):

            if x < gridXsplit:
                gridSector_g2.append( 0 )
            elif x >= gridXsplit:
                gridSector_g2.append( 1 )
            else:
                printf(" Invalid grid coordinates! ")

print( grid2Size  )
print( len( gridSector_g2 ) )

# Write local grid G2
DomainGrid( gridXsize - gridXsplit + 1, gridYsize, gridZsize,
            gridXbegin + ( gridXsplit - 1 ) * searchRadius_h, gridYbegin, gridZbegin,
            gridSector_g2,
            "distributedGrid_g2.vtk" )


#====================================================================================================

with open( 'template/ParticlesConfig_template.h', 'r' ) as file :
  fileParticleConf = file.read()

# Replace the target string
#fileParticleConf = fileParticleConf.replace( 'placeholderSubdomain', 'G0' )
fileParticleConf = fileParticleConf.replace( 'placeholderDimension', str( spaceDimension ) )

fileParticleConf = fileParticleConf.replace( 'placeholderFluidParticlesG0', str( len( fluid_rx_g1 ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedFluidParticlesG0', str( len( fluid_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderBoundaryParticlesG0', str( len( box_rx_g1 ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedBoundaryParticlesG0', str( len( box_rx ) ) )

fileParticleConf = fileParticleConf.replace( 'placeholderFluidParticlesG1', str( len( fluid_rx_g2 ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedFluidParticlesG1', str( len( fluid_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderBoundaryParticlesG1', str( len( box_rx_g2 ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedBoundaryParticlesG1', str( len( box_rx ) ) )


fileParticleConf = fileParticleConf.replace( 'placeholderSearchRadius', str( searchRadius ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXSize', str( gridXsize ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYSize', str( gridYsize ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridZSize', str( gridZsize ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXBegin', str( round( gridXbegin, 9  ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYBegin', str( round( gridYbegin, 9  ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridZBegin', str( round( gridZbegin, 9  ) ) )

# Distributed subdomain info G0
fileParticleConf = fileParticleConf.replace( 'placeholderParticleIdxStartG0', str( 0 ) )
fileParticleConf = fileParticleConf.replace( 'placeholderParticleIdxRealStartG0', str( 0 ) )
fileParticleConf = fileParticleConf.replace( 'placeholderParticleIdxEndG0', str( counter_g1_nptcs - 1  ) )
fileParticleConf = fileParticleConf.replace( 'placeholderParticleIdxRealEndG0', str( counter_g1_nptcsReal - 1) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridIdxOverlapStartG0', str( 0 ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridIdxStartG0', str( 0 ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridIdxOverlapEndG0', str( gridXsplit ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridIdxEndG0', str( gridXsplit - 1 ) )

# Distributed subdomain info G1
fileParticleConf = fileParticleConf.replace( 'placeholderParticleIdxStartG1', str( 0 ) )
fileParticleConf = fileParticleConf.replace( 'placeholderParticleIdxRealStartG1', str( 0 ) )
fileParticleConf = fileParticleConf.replace( 'placeholderParticleIdxEndG1', str( counter_g1_nptcs - 1  ) )
fileParticleConf = fileParticleConf.replace( 'placeholderParticleIdxRealEndG1', str( counter_g1_nptcsReal - 1) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridIdxOverlapStartG1', str( gridXsplit - 1 ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridIdxStartG1', str( gridXsplit ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridIdxOverlapEndG1', str( gridXsize ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridIdxEndG1', str( gridXsize ) )

# Write the file out again
with open( 'ParticlesConfig.h', 'w' ) as file:
  file.write( fileParticleConf )

#----------------------------------------------------------------------------------------------------
