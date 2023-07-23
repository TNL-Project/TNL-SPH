#---------------------------------------------------------------------------#
#
# case: distributedDamBreak2D_WCSPH-DBC
# case destription: 2D dam break problem, using 1d partitioning and MPI to
#                   utilize multiple GPUs
#
#---------------------------------------------------------------------------#
### Parameters of the case necessary for case creation:

# Dimensions of the box (wall), and the block of fluid [m]:
boxL = 1.61
boxH = 0.8

fluidL = 0.6
fluidH = 0.3

# Initial particle distance (dp)[m]:
dp = 0.002

# Smoothing length coefitient:
# - smoothing length (h)[m] = smoothing length coef (Coef_h)[-] * initial particle distance (d_p)[m]
# ( h = Coef_h * dp )
smoothingLentghCoef = 2**0.5

# Referential density of the medium (rho0)[kg/m^3]:
rho0 = 1000.

# Number of boundary layers (n_layer)[-]:
numberOfBoundaryLayers = 3 #TODO: Compute this automatically

# Numerical speed of sound (c0)[m/s]:
speedOfSound = 34.3

# Initial time step (dtInit)[s].
# - in case, that initial time step is not defined, is computed automatically.
timeStep = 0.00002 #otherwise is obtained automatically

# CFL number (CFL)[-]:
CFLnumber = 0.2

# Setup number of subdomains (devices), to split the problem [-]:
numberOfSubdomains = 2

# Print information during case generation
printInfoString = False
#---------------------------------------------------------------------------#

import os
resultsPath = r'./results'
if not os.path.exists( resultsPath ):
    os.makedirs( resultsPath )

sourcesPath = r'./sources'
if not os.path.exists( sourcesPath ):
    os.makedirs( sourcesPath )

### Generate geometry
boxL_n = round( boxL / dp )
boxH_n = round( boxH / dp )

fluidL_n = round( fluidL / dp )
fluidH_n = round( fluidH / dp )

# Generate fluid particles
fluid_rx = []; fluid_ry = []; fluid_rz = []

for x in range( fluidL_n ):
    for z in range( fluidH_n ):
        fluid_rx.append( dp * ( x + 1 ) )
        fluid_ry.append( 0. ) #we use only 2D case
        fluid_rz.append( dp * ( z + 1 ) )

        #fluid_rx.append( boxL - fluidL + dp * x  )
        #fluid_ry.append( 0. ) #we use only 2D case
        #fluid_rz.append( dp * ( z + 1 ) )

# Generate boundary particles
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

x_last = box_rx[ -1 - ( numberOfBoundaryLayers - 1 ) ] #due to discretisation, we need to save last value of bottom wall

# right wall
for layer in range( numberOfBoundaryLayers ):
    for z in range( boxH_n - 1 ):
        box_rx.append( x_last + dp * layer )
        box_ry.append( 0. ) #we use only 2D case
        box_rz.append( ( z + 1 ) * dp )

### Write fluid particles
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
saveParticlesVTK.save_polydata( fluidToWrite, "sources/dambreak_fluid.vtk" )

r = np.array( ( box_rx, box_rz, box_ry ), dtype=float ).T #!!
v = np.zeros( ( len( box_rx ), 3 ) )
rho = rho0 * np.ones( len( box_rx ) )
p = np.zeros( len( box_rx ) )
ptype = np.ones( len( box_rx ) )

boxToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
saveParticlesVTK.save_polydata( boxToWrite, "sources/dambreak_boundary.vtk" )

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
gridXbegin = 1.005 * ( min( min( fluid_rx ), min( box_rx ) )  - searchRadius )
gridYbegin = 1.005 * ( min( min( fluid_rz ), min( box_rz ) )  - searchRadius )

gridXend = 1.005 * ( max( max( fluid_rx ), max( box_rx ) ) + searchRadius )
gridYend = 1.005 * ( max( max( fluid_rz ), max( box_rz ) ) + searchRadius )

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
with open( 'sources/SPHCaseConfig.h', 'w' ) as file:
  file.write( fileSPHConf )


# Read and write (with possible edit) measuretool config file.
with open( 'template/MeasuretoolConfig.h', 'r' ) as file :
  fileMeasuretoolConf = file.read()

with open( 'sources/MeasuretoolConfig.h', 'w' ) as file:
  file.write( fileMeasuretoolConf )

### Divide the domain into subdomains
numberOfPtcsTotal = len( fluid_rx )
searchRadius_h = round( smoothingLentgh * 2 , 7 )

numerOfFluidParticlesPerSubdomain = ( int )( numberOfPtcsTotal / numberOfSubdomains )
print( f'midleParticleNumber: {numerOfFluidParticlesPerSubdomain}' )

gridSplits = []
for subdomain in range( numberOfSubdomains - 1 ):
    gridSplits.append( math.ceil( fluid_rx[ numerOfFluidParticlesPerSubdomain * ( subdomain + 1 ) ] / searchRadius ) )
#np.sort( gridSplits )

print ( f'Grid splits: {gridSplits}' )

print( f'Grid global:' )
print( f'Grid global - size: [ {gridXsize}, {gridYsize} ]' )
print( f'Grid global - begin: [ {gridXbegin}, {gridYbegin} ]\n' )

gridSizes = []
gridOrigins = []; gridIndexOrigins = []
for subdomain in range( numberOfSubdomains ):
    if subdomain == 0:
        gridSizes.append( gridSplits[ subdomain ] - 0 )
        gridOrigins.append( gridXbegin )
        gridIndexOrigins.append( 0 )
    elif subdomain == numberOfSubdomains - 1:
        gridSizes.append( gridXsize - gridSplits[ subdomain - 1 ] )
        gridOrigins.append( gridXbegin + gridSplits[ subdomain - 1 ] * searchRadius )
        gridIndexOrigins.append( gridSplits[ subdomain - 1 ] )
    else:
        gridSizes.append( gridSplits[ subdomain - 1 ] - gridIndexOrigins[ subdomain -1 ] )
        gridOrigins.append( gridXbegin + gridSplits[ subdomain - 1 ] * searchRadius )
        gridIndexOrigins.append( gridSplits[ subdomain - 1 ] )

    print( f'Grid subdomain: {subdomain}' )
    print( f'Grid subdomain - size: [ {gridSizes[ subdomain ]}, {gridYsize} ]' )
    print( f'Grid subdomain - origin: [ {gridOrigins[ subdomain ]}, {gridYbegin} ]' )
    print( f'Grid subdomain - index origin: [ {gridIndexOrigins[ subdomain ]}, {gridYbegin} ]\n' )

print( f'Grid size: {gridXsize} gridSizeSplits: {np.sum( gridSizes )}\n' )

#=====================================================================================================
# Save subdomain grid function.

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

#=====================================================================================================

subdomainStringsArrays = []
subdomainStringTemplate = """
      //Subdomain #placeholderSubdomainNumber
      particlesParams[ #placeholderSubdomainNumber ].numberOfParticles = #placeholderFluidParticles;
      particlesParams[ #placeholderSubdomainNumber ].numberOfAllocatedParticles = #placeholderAllocatedFluidParticles;
      particlesParams[ #placeholderSubdomainNumber ].numberOfBoundaryParticles = #placeholderBoundaryParticles;
      particlesParams[ #placeholderSubdomainNumber ].numberOfAllocatedBoundaryParticles = #placeholderAllocatedBoundaryParticles;

      particlesParams[ #placeholderSubdomainNumber ].searchRadius = #placeholderSearchRadiusf * 1.001f;
      particlesParams[ #placeholderSubdomainNumber ].gridXsize = #placeholderGridXSize;
      particlesParams[ #placeholderSubdomainNumber ].gridYsize = #placeholderGridYSize;
      particlesParams[ #placeholderSubdomainNumber ].gridOrigin = { #placeholderGridXBeginf, #placeholderGridYBeginf };

      particlesParams[ #placeholderSubdomainNumber ].gridSize = { particlesParams[ #placeholderSubdomainNumber ].gridXsize, particlesParams[ #placeholderSubdomainNumber ].gridYsize };
      particlesParams[ #placeholderSubdomainNumber ].numberOfGridCells = particlesParams[ #placeholderSubdomainNumber ].gridXsize * particlesParams[ #placeholderSubdomainNumber ].gridYsize;

      //Subdomain #placeholderSubdomainNumber - Subdomain info
      subdomainParams[ #placeholderSubdomainNumber ].particleIdxStart = #placeholderParticleIdxStart;
      subdomainParams[ #placeholderSubdomainNumber ].particleIdxRealStart = #placeholderParticleIdxRealStart;
      subdomainParams[ #placeholderSubdomainNumber ].particleIdxEnd = #placeholderParticleIdxEnd;
      subdomainParams[ #placeholderSubdomainNumber ].particleIdxRealEnd = #placeholderParticleIdxRealEnd;

      subdomainParams[ #placeholderSubdomainNumber ].gridIdxOverlapStar = #placeholderGridIdxOverlapStart;
      subdomainParams[ #placeholderSubdomainNumber ].gridIdxStart = #placeholderGridIdxStart;
      subdomainParams[ #placeholderSubdomainNumber ].gridIdxOverlapEnd = #placeholderGridIdxOverlapEnd;
      subdomainParams[ #placeholderSubdomainNumber ].gridIdxEnd = #placeholderGridIdxEnd;
"""

def generateSubdomain( subdomain ):
    #Fields and variables to set
    subdomain_fluid_rx = []
    subdomain_fluid_ry = []
    subdomain_fluid_rz = []
    subdomain_fluid_ptype = []

    subdomain_box_rx = []
    subdomain_box_ry = []
    subdomain_box_rz = []
    subdomain_box_ptype = []

    subdomain_counter_nptcs = 0
    subdomain_counter_nptcsReal = 0

    subdomain_gridCoordinates_x = []
    subdomain_gridCoordinates_y = []
    subdomain_gridCoordinates_z = []
    subdomain_gridSector = []

    #Load the limits of current subdomain
    if subdomain == 0:
        lowerPositionLimit = gridOrigins[ subdomain ]
        upperPositionLimit = gridOrigins[ subdomain + 1 ]
    elif subdomain == numberOfSubdomains - 1:
        #lowerPositionLimit = gridOrigins[ subdomain ]
        #upperPositionLimit = gridXbegin + ( gridXsize + 1 ) * searchRadius_h
        lowerPositionLimit = gridOrigins[ subdomain ] - searchRadius_h #TODO
        upperPositionLimit = gridXbegin + ( gridXsize + 1 ) * searchRadius_h
    else:
        lowerPositionLimit = gridOrigins[ subdomain ]
        upperPositionLimit = gridOrigins[ subdomain + 1 ]

    #Prepare particle fields for given subdomain - fluid
    for i in range ( len( fluid_rx ) ):
        if( fluid_rx[ i ] > lowerPositionLimit ) and ( fluid_rx[ i ] <= upperPositionLimit ):
            subdomain_fluid_rx.append( fluid_rx[ i ] )
            subdomain_fluid_ry.append( fluid_ry[ i ] )
            subdomain_fluid_rz.append( fluid_rz[ i ] )
            subdomain_fluid_ptype.append( 0 )
            subdomain_counter_nptcs += 1
            subdomain_counter_nptcsReal += 1

        if( fluid_rx[ i ] > upperPositionLimit and fluid_rx[ i ] <= upperPositionLimit + searchRadius_h ):
            subdomain_fluid_rx.append( fluid_rx[ i ] )
            subdomain_fluid_ry.append( fluid_ry[ i ] )
            subdomain_fluid_rz.append( fluid_rz[ i ] )
            subdomain_fluid_ptype.append( 1 )
            subdomain_counter_nptcsReal += 1

    #Prepare particle fields for given subdomain - boundary
    for i in range ( len( box_rx ) ):
        if( box_rx[ i ] > lowerPositionLimit ) and ( box_rx[ i ] <= upperPositionLimit ):
            subdomain_box_rx.append( box_rx[ i ] )
            subdomain_box_ry.append( box_ry[ i ] )
            subdomain_box_rz.append( box_rz[ i ] )
            subdomain_box_ptype.append( 0 )

        if( box_rx[ i ] > upperPositionLimit and box_rx[ i ] <= upperPositionLimit + searchRadius_h ):
            subdomain_box_rx.append( box_rx[ i ] )
            subdomain_box_ry.append( box_ry[ i ] )
            subdomain_box_rz.append( box_rz[ i ] )
            subdomain_box_ptype.append( 1 )

    #Write particle for given subdomain - fluid
    r = np.array( ( subdomain_fluid_rx, subdomain_fluid_rz, subdomain_fluid_ry ), dtype=float ).T #!!
    v = np.zeros( ( len( subdomain_fluid_rx ), 3 ) )
    rho = rho0 * np.ones( len( subdomain_fluid_rx ) )
    p = np.zeros( len( subdomain_fluid_rx ) )
    ptype = np.array( ( subdomain_fluid_ptype ), dtype=float ).T

    fluidToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
    subdomain_fluid_outputname = "sources/dambreak_fluid_subdomain" + str( subdomain ) + '.vtk'
    print( f'Subdomain fluid ouputfilename: {subdomain_fluid_outputname}' )
    saveParticlesVTK.save_polydata( fluidToWrite, subdomain_fluid_outputname )

    #Write particle for given subdomain - boundary
    r = np.array( ( subdomain_box_rx, subdomain_box_rz, subdomain_box_ry ), dtype=float ).T #!!
    v = np.zeros( ( len( subdomain_box_rx ), 3 ) )
    rho = rho0 * np.ones( len( subdomain_box_rx ) )
    p = np.zeros( len( subdomain_box_rx ) )
    ptype = np.array( ( subdomain_box_ptype ), dtype=float ).T

    boxToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
    subdomain_boundary_outputname = "sources/dambreak_boundary_subdomain" + str( subdomain ) + '.vtk'
    print( f'Subdomain boundary ouputfilename: {subdomain_boundary_outputname}' )
    saveParticlesVTK.save_polydata( boxToWrite, subdomain_boundary_outputname )

    #Prepare informations about iven subdomain
    from copy import copy
    infoString = copy( subdomainStringTemplate )

    infoString = infoString.replace( '#placeholderSubdomainNumber', str( subdomain ) )

    infoString = infoString.replace( '#placeholderSearchRadius', str( searchRadius ) )
    infoString = infoString.replace( '#placeholderGridXSize', str( gridXsize ) )
    infoString = infoString.replace( '#placeholderGridYSize', str( gridYsize ) )
    infoString = infoString.replace( '#placeholderGridXBegin', str( round( gridXbegin, 9  ) ) )
    infoString = infoString.replace( '#placeholderGridYBegin', str( round( gridYbegin, 9  ) ) )

    infoString = infoString.replace( '#placeholderFluidParticles', str( len( subdomain_fluid_rx ) ) )
    infoString = infoString.replace( '#placeholderAllocatedFluidParticles', str( len( fluid_rx ) ) )
    infoString = infoString.replace( '#placeholderBoundaryParticles', str( len( subdomain_box_rx ) ) )
    infoString = infoString.replace( '#placeholderAllocatedBoundaryParticles', str( len( box_rx ) ) )

    infoString = infoString.replace( '#placeholderParticleIdxStart', str( 0 ) )
    infoString = infoString.replace( '#placeholderParticleIdxRealStart', str( 0 ) )
    infoString = infoString.replace( '#placeholderParticleIdxEnd', str( subdomain_counter_nptcs - 1  ) )
    infoString = infoString.replace( '#placeholderParticleIdxRealEnd', str( subdomain_counter_nptcsReal - 1) )

    if subdomain == 0:
        gridIdxOverlapStart = 0
        gridIdxStart = 0
        gridIdxOverlapEnd = gridIndexOrigins[ subdomain + 1 ]
        gridIdxEnd = gridIndexOrigins[ subdomain + 1 ] - 1
    elif subdomain == numberOfSubdomains - 1:
        gridIdxOverlapStart = gridIndexOrigins[ subdomain ] - 1
        gridIdxStart = gridIndexOrigins[ subdomain ]
        gridIdxOverlapEnd = gridXsize
        gridIdxEnd = gridXsize
    else:
        gridIdxOverlapStart = gridIndexOrigins[ subdomain ] - 1
        gridIdxStart = gridIndexOrigins[ subdomain ]
        gridIdxOverlapEnd = gridIndexOrigins[ subdomain + 1 ]
        gridIdxEnd = gridIndexOrigins[ subdomain + 1 ] - 1

    infoString = infoString.replace( '#placeholderGridIdxOverlapStart', str( gridIdxOverlapStart ) )
    infoString = infoString.replace( '#placeholderGridIdxStart', str( gridIdxStart ) )
    infoString = infoString.replace( '#placeholderGridIdxOverlapEnd', str( gridIdxOverlapEnd ) )
    infoString = infoString.replace( '#placeholderGridIdxEnd', str( gridIdxEnd ) )

    if( printInfoString ): print( f'Subdomain: {subdomain} info string:\n{subdomainsString}' )
    subdomainStringsArrays.append( infoString )

    #Generate grid for given Subdomain
    #TODO: Make this nicer
    for y in range ( gridYsize ):
        for x in range ( gridSizes[ subdomain ] + 1 ):

            #First subdomain:
            if subdomain == 0:
                if x < gridSizes[ 0 ]:
                    subdomain_gridSector.append( 0 )
                else:
                    subdomain_gridSector.append( 1 )

            #Last subdomain:
            elif subdomain == numberOfSubdomains - 1:
                if x == 0 :
                    subdomain_gridSector.append( subdomain - 1 )
                else:
                    subdomain_gridSector.append( subdomain )

            #Subdomain in the middle:
            else:
                if x == 0 :
                    subdomain_gridSector.append( subdomain - 1 )
                elif x == gridSizes[ subdomain ]:
                    subdomain_gridSector.append( subdomain )
                else:
                    subdomain_gridSector.append( subdomain + 1 )

    #TODO: Add local verlap to grid begin:
    gridXOriginWithOverlap = gridOrigins[ subdomain ]
    if subdomain > 0: gridXOriginWithOverlap -= searchRadius_h

    # Write local grid G1
    subdomain_grid_outputname = "sources/dambreak_grid_subdomain" + str( subdomain ) + '.vtk'
    DomainGrid( gridSizes[ subdomain ] + 1, gridYsize, 1,       # grid size
                #gridOrigins[ subdomain ], gridYbegin, 0,       # coordinates of grid origin
                gridXOriginWithOverlap, gridYbegin, 0,          # coordinates of grid origin
                subdomain_gridSector,                           # array with index of grid sector
                subdomain_grid_outputname )                     # outputfile name

#====================================================================================================

gridCoordinates_x = []
gridCoordinates_y = []
gridCoordinates_z = []

gridSector = []

#---------------------------------------------------------------------------#

subdomainsString = ''
inputFluidFilesString = ''
inputBoundaryFilesString = ''

for subdomain in range( numberOfSubdomains ):
    generateSubdomain( subdomain )
    subdomainsString += subdomainStringsArrays[ subdomain ]

    inputFluidFilesString += "\"sources/dambreak_fluid_subdomain" + str( subdomain ) + '.vtk\"'
    inputBoundaryFilesString += "\"sources/dambreak_boundary_subdomain" + str( subdomain ) + '.vtk\"'

    if subdomain < numberOfSubdomains - 1:
        inputFluidFilesString += ', '
        inputBoundaryFilesString += ', '

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

if( printInfoString ): print( f'Subdomain string:\n{subdomainsString}' )

with open( 'template/ParticlesConfig_template.h', 'r' ) as file :
  fileParticleConf = file.read()
  fileParticleConf = fileParticleConf.replace( '#placeholderSubdomainInfo', subdomainsString )
  fileParticleConf = fileParticleConf.replace( '#placeholderNumberOfSubdomains', str( numberOfSubdomains ) )

# Write the file out again
with open( 'sources/ParticlesConfig.h', 'w' ) as file:
  file.write( fileParticleConf )

#---------------------------------------------------------------------------#
