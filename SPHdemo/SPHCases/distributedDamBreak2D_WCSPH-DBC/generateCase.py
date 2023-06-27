# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#
# damBreak2D_WCSPH-DBC_benchmark
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

boxL = 1.61
boxH = 0.8

fluidL = 0.6
fluidH = 0.3

dp = 0.002
smoothingLentghCoef = 2**0.5

rho0 = 1000.
p0 = 0.

numberOfBoundaryLayers = 3

speedOfSound = 34.3
CFLnumber = 0.2
timeStep = 0.00002 #otherwise is obtained automatically

numberOfProcessors = 2


write = '.vtk' #.ptcs or .vtk

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

boxL_n = round( boxL / dp )
boxH_n = round( boxH / dp )

fluidL_n = round( fluidL / dp )
fluidH_n = round( fluidH / dp )

### Generate fluid particles
fluid_rx = []; fluid_ry = []; fluid_rz = []

for x in range( fluidL_n ):
    for z in range( fluidH_n ):
        fluid_rx.append( dp * ( x + 1 ) )
        fluid_ry.append( 0. ) #we use only 2D case
        fluid_rz.append( dp * ( z + 1 ) )

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
    for x in range( boxL_n + ( numberOfBoundaryLayers - 1 ) * 2 + 1 ):
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
    with open( "dambreak_fluid.ptcs", "w" ) as f:
        f.write( str( len( fluid_rx ) ) + "\n" )
        for i in range( len( fluid_rx ) ):
            f.write( str( round( fluid_rx[ i ], 5 ) ) + " " + str( round( fluid_rz[i], 5 ) ) + " " + \
                     str( round( fluid_ry[ i ], 5 ) ) + " " + str( 0. ) + " " + str( 0. ) + " " + str( 0. ) + " " + \
                     str( round( rho0, 5 ) ) + " " + str( round( p0, 5 ) ) + " " + str( 0 ) + "\n" )

    ### Write fluid particles
    with open("dambreak_boundary.ptcs", "w") as f:
        f.write( str( len( box_rx ) ) + "\n" )
        for i in range( len( box_rx ) ):
            f.write( str( round( box_rx[ i ], 5 ) ) + " " + str( round( box_rz[ i ], 5 ) ) + " " + \
                     str( round( box_ry[ i ], 5 ) ) + " " + str( 0. ) + " " + str( 0. ) + " " + str( 0. ) + " " + \
                     str( round( rho0, 5 ) ) + " " + str( round( p0, 5 ) ) + " " + str( 1 ) + "\n" )
elif write == '.vtk':
    import sys
    sys.path.append('../../tools/')
    import saveParticlesVTK
    import numpy as np

    r = np.array( ( fluid_rx, fluid_rz, fluid_ry ), dtype=float ).T #!!
    v = np.zeros( ( len( fluid_rx ), 3 ) )
    rho = rho0 * np.ones( len( fluid_rx ) )
    p = np.zeros( len( fluid_rx ) )
    ptype = np.zeros( len( fluid_rx ) )

    fluidToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
    saveParticlesVTK.save_polydata( fluidToWrite, "dambreak_fluid.vtk" )

    r = np.array( ( box_rx, box_rz, box_ry ), dtype=float ).T #!!
    v = np.zeros( ( len( box_rx ), 3 ) )
    rho = rho0 * np.ones( len( box_rx ) )
    p = np.zeros( len( box_rx ) )
    ptype = np.ones( len( box_rx ) )

    boxToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
    saveParticlesVTK.save_polydata( boxToWrite, "dambreak_boundary.vtk" )
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
fileParticleConf = fileParticleConf.replace( '_placeholderSubdomain', '' )
fileParticleConf = fileParticleConf.replace( 'placeholderDimension', str( spaceDimension ) )
fileParticleConf = fileParticleConf.replace( 'placeholderFluidParticles', str( len( fluid_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedFluidParticles', str( len( fluid_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderBoundaryParticles', str( len( box_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedBoundaryParticles', str( len( box_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderSearchRadius', str( searchRadius ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXSize', str( gridXsize ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYSize', str( gridYsize ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXBegin', str( round( gridXbegin, 9  ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYBegin', str( round( gridYbegin, 9  ) ) )

# Write the file out again
with open( 'ParticlesConfig.h', 'w' ) as file:
  file.write( fileParticleConf )

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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# split
#tot: numberOfPtcsTotal = len( fluid_rx ) + len( box_rx )
# split based on fluid:
numberOfPtcsTotal = len( fluid_rx )

middleParticle = ( int )( numberOfPtcsTotal / 2 )
print(" midleParticleNumber: ", middleParticle )

searchRadius_h = round( smoothingLentgh * 2 , 7 )
gridXsplit = math.ceil( fluid_rx[ middleParticle ] / searchRadius_h )
print(" searchRadius: ", gridXsplit )

print( "Grid-global size: [ ", gridXsize , ",", gridYsize, " ]." )
print( "Grid-global begin: [ ", gridXbegin , ",", gridYbegin, " ]." )

gridXsplitBegin = gridXbegin + gridXsplit * searchRadius
print( "Grid-G1 size: [ ", gridXsplit , ",", gridYsize, " ]." )
print( "Grid-G1 begin: [ ", gridXbegin , ",", gridYbegin, " ]." )

gridXsplitBegin = gridXbegin + gridXsplit * searchRadius
print( "Grid-G2 size: [ ", gridXsize - gridXsplit , ",", gridYsize, " ]." )
print( "Grid-G2 size: [ ", gridXsplitBegin , ",", gridYsize, " ]." )

#nbscell# = [int(cidx[i]) for i in (neighbors(x, y)) if cellFullList_ptcsidx_first[int(cidx[i])] != 0]

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

rg1 = np.array( ( fluid_rx_g1, fluid_rz_g1, fluid_ry_g1 ), dtype=float ).T #!!
vg1 = np.zeros( ( len( fluid_rx_g1 ), 3 ) )
rhog1 = rho0 * np.ones( len( fluid_rx_g1 ) )
pg1 = np.zeros( len( fluid_rx_g1 ) )
#ptypeg1 = np.zeros( len( fluid_rx_g1 ) )
ptypeg1 = np.array( ( ptype_fluid_g1 ), dtype=float ).T

fluidToWrite = saveParticlesVTK.create_pointcloud_polydata( rg1, vg1, rhog1, pg1, ptypeg1 )
saveParticlesVTK.save_polydata( fluidToWrite, "dambreak_fluid_g1.vtk" )

rg1 = np.array( ( box_rx_g1, box_rz_g1, box_ry_g1 ), dtype=float ).T #!!
vg1 = np.zeros( ( len( box_rx_g1 ), 3 ) )
rhog1 = rho0 * np.ones( len( box_rx_g1 ) )
pg1 = np.zeros( len( box_rx_g1 ) )
#ptypeg1 = np.ones( len( box_rx_g1 ) )
ptypeg1 = np.array( ( ptype_box_g1 ), dtype=float ).T

boxToWrite = saveParticlesVTK.create_pointcloud_polydata( rg1, vg1, rhog1, pg1, ptypeg1 )
saveParticlesVTK.save_polydata( boxToWrite, "dambreak_boundary_g1.vtk" )

#----------------------------------------------------------------------------------------------------

with open( 'template/ParticlesConfig_template.h', 'r' ) as file :
  fileParticleConf = file.read()

# Replace the target string
fileParticleConf = fileParticleConf.replace( 'placeholderSubdomain', 'g1' )
fileParticleConf = fileParticleConf.replace( 'placeholderDimension', str( spaceDimension ) )
fileParticleConf = fileParticleConf.replace( 'placeholderFluidParticles', str( len( fluid_rx_g1 ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedFluidParticles', str( len( fluid_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderBoundaryParticles', str( len( box_rx_g1 ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedBoundaryParticles', str( len( box_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderSearchRadius', str( searchRadius ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXSize', str( gridXsize ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYSize', str( gridYsize ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXBegin', str( round( gridXbegin, 9  ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYBegin', str( round( gridYbegin, 9  ) ) )

# Write the file out again
with open( 'ParticlesConfig_g1.h', 'w' ) as file:
  file.write( fileParticleConf )

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

rg2 = np.array( ( fluid_rx_g2, fluid_rz_g2, fluid_ry_g2 ), dtype=float ).T #!!
vg2 = np.zeros( ( len( fluid_rx_g2 ), 3 ) )
rhog2 = rho0 * np.ones( len( fluid_rx_g2 ) )
pg2 = np.zeros( len( fluid_rx_g2 ) )
#ptypeg2 = np.zeros( len( fluid_rx_g2 ) )
ptypeg2 = np.array( ( ptype_fluid_g2 ), dtype=float ).T

fluidToWrite = saveParticlesVTK.create_pointcloud_polydata( rg2, vg2, rhog2, pg2, ptypeg2 )
saveParticlesVTK.save_polydata( fluidToWrite, "dambreak_fluid_g2.vtk" )

rg2 = np.array( ( box_rx_g2, box_rz_g2, box_ry_g2 ), dtype=float ).T #!!
vg2 = np.zeros( ( len( box_rx_g2 ), 3 ) )
rhog2 = rho0 * np.ones( len( box_rx_g2 ) )
pg2 = np.zeros( len( box_rx_g2 ) )
#ptypeg2 = np.ones( len( box_rx_g2 ) )
ptypeg2 = np.array( ( ptype_box_g2 ), dtype=float ).T

boxToWrite = saveParticlesVTK.create_pointcloud_polydata( rg2, vg2, rhog2, pg2, ptypeg2 )
saveParticlesVTK.save_polydata( boxToWrite, "dambreak_boundary_g2.vtk" )

#---------------------------------------------------------------------------------------------------

with open( 'template/ParticlesConfig_template.h', 'r' ) as file :
  fileParticleConf = file.read()

# Replace the target string
fileParticleConf = fileParticleConf.replace( 'placeholderSubdomain', 'g2' )
fileParticleConf = fileParticleConf.replace( 'placeholderDimension', str( spaceDimension ) )
fileParticleConf = fileParticleConf.replace( 'placeholderFluidParticles', str( len( fluid_rx_g2 ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedFluidParticles', str( len( fluid_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderBoundaryParticles', str( len( box_rx_g2 ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderAllocatedBoundaryParticles', str( len( box_rx ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderSearchRadius', str( searchRadius ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXSize', str( gridXsize ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYSize', str( gridYsize ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridXBegin', str( round( gridXbegin, 9  ) ) )
fileParticleConf = fileParticleConf.replace( 'placeholderGridYBegin', str( round( gridYbegin, 9  ) ) )

# Write the file out again
with open( 'ParticlesConfig_g2.h', 'w' ) as file:
  file.write( fileParticleConf )

#====================================================================================================


with open( 'template/DistributedSimulationConfig_template.h', 'r' ) as file :
  fileSubdomains = file.read()

fileSubdomains = fileSubdomains.replace( 'placeholderParticleIdxStartG1', str( 0 ) )
fileSubdomains = fileSubdomains.replace( 'placeholderParticleIdxRealStartG1', str( 0 ) )
fileSubdomains = fileSubdomains.replace( 'placeholderParticleIdxEndG1', str( counter_g1_nptcs - 1  ) )
fileSubdomains = fileSubdomains.replace( 'placeholderParticleIdxRealEndG1', str( counter_g1_nptcsReal - 1) )
fileSubdomains = fileSubdomains.replace( 'placeholderGridIdxOverlapStartG1', str( 0 ) )
fileSubdomains = fileSubdomains.replace( 'placeholderGridIdxStartG1', str( 0 ) )
fileSubdomains = fileSubdomains.replace( 'placeholderGridIdxOverlapEndG1', str( gridXsplit ) )
fileSubdomains = fileSubdomains.replace( 'placeholderGridIdxEndG1', str( gridXsplit - 1 ) )

fileSubdomains = fileSubdomains.replace( 'placeholderParticleIdxStartG2', str( 0 ) )
fileSubdomains = fileSubdomains.replace( 'placeholderParticleIdxRealStartG2', str( 0 ) )
fileSubdomains = fileSubdomains.replace( 'placeholderParticleIdxEndG2', str( counter_g1_nptcs - 1  ) )
fileSubdomains = fileSubdomains.replace( 'placeholderParticleIdxRealEndG2', str( counter_g1_nptcsReal - 1) )
fileSubdomains = fileSubdomains.replace( 'placeholderGridIdxOverlapStartG2', str( gridXsplit - 1 ) )
fileSubdomains = fileSubdomains.replace( 'placeholderGridIdxStartG2', str( gridXsplit ) )
fileSubdomains = fileSubdomains.replace( 'placeholderGridIdxOverlapEndG2', str( gridXsize ) )
fileSubdomains = fileSubdomains.replace( 'placeholderGridIdxEndG2', str( gridXsize ) )

# Write the file out again
with open( 'DistributedSimulationConfig.h', 'w' ) as file:
  file.write( fileSubdomains )

#====================================================================================================

gridCoordinates_x = []
gridCoordinates_y = []
gridCoordinates_z = []

gridSector = []

for y in range ( gridYsize ):
    for x in range ( gridXsize ):
        gridCoordinates_x.append( x * searchRadius_h )
        gridCoordinates_y.append( y * searchRadius_h )
        gridCoordinates_z.append( 0. )

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
            print( "DIMENSIONS ", gridXsize + 1 , " ", gridYsize + 1, " ", 1 )
            #print( "POINTS ", gridXsize * gridYsize * 1, " float" )
            #for i in range( gridXsize * gridYsize * 1 ):
            #    print( "", gridCoordinates_x[ i ], " ", gridCoordinates_y[ i ], " ", gridCoordinates_z[ i ] )
            #print( "ASPECT_RATIO 1 1 1" )
            print( "ASPECT_RATIO ", searchRadius_h , " ", searchRadius_h , " ",  searchRadius_h )
            print( "ORIGIN ", gridXbegin , " ", gridYbegin , " ",  0  )
            ##print( "POINT_DATA ",  gridXsize * gridYsize * 1  )
            print( "CELL_DATA ",  gridXsize * gridYsize * 1  )
            #print( "FIELD FieldData 1" )
            #print( "GridSector 1",  gridXsize * gridYsize * 1 , " int" )
            print( "SCALARS GridSector int 1 ")
            print( "LOOKUP_TABLE default" )
            for i in range( gridXsize * gridYsize * 1 ):
                print( gridSector[ i ] )
            #print( "VECTORS GridIndex int 3 ")
            #print( "LOOKUP_TABLE default" )
            #for i in range( gridXsize * gridYsize * 1 ):
            #    print( gridSector[ i ], " ", gridSector[ i ], " ", gridSector[ i ] )
    print("Done.")

# Write global grid as example
DomainGrid( gridXsize, gridYsize, 1,    # grid size
            gridXbegin, gridYbegin, 0,  # coordinates of rgrid origina
            gridSector,                 # array with index of grid sector
            "distributedGrid.vtk" )     # outputfile name

gridSector_g1 = []

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
DomainGrid( gridXsplit + 1, gridYsize, 1,   # grid size
            gridXbegin, gridYbegin, 0,      # coordinates of rgrid origina
            gridSector_g1,                  # array with index of grid sector
            "distributedGrid_g1.vtk" )      # outputfile name

gridSector_g2 = []

grid2Size = ( gridXsize - gridXsplit + 1 ) * gridYsize * 1

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
DomainGrid( gridXsize - gridXsplit + 1, gridYsize, 1,
            gridXbegin + ( gridXsplit - 1 ) * searchRadius_h, gridYbegin, 0,
            gridSector_g2,
            "distributedGrid_g2.vtk" )


#====================================================================================================
