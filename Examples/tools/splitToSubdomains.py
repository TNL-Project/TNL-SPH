import math
import numpy as np

subdomainStringTemplate = """
      //Subdomain #placeholderSubdomainNumber
      particlesParams[ #placeholderSubdomainNumber ].numberOfParticles = #placeholderFluidParticles;
      particlesParams[ #placeholderSubdomainNumber ].numberOfAllocatedParticles = #placeholderAllocatedFluidParticles;
      particlesParams[ #placeholderSubdomainNumber ].numberOfBoundaryParticles = #placeholderBoundaryParticles;
      particlesParams[ #placeholderSubdomainNumber ].numberOfAllocatedBoundaryParticles = #placeholderAllocatedBoundaryParticles;

      particlesParams[ #placeholderSubdomainNumber ].searchRadius = #placeholderSearchRadiusf * 1.001f;
      particlesParams[ #placeholderSubdomainNumber ].gridXsize = #placeholderGridXSize;
      particlesParams[ #placeholderSubdomainNumber ].gridYsize = #placeholderGridYSize;
      particlesParams[ #placeholderSubdomainNumber ].gridZsize = #placeholderGridZSize;
      particlesParams[ #placeholderSubdomainNumber ].gridOrigin = { #placeholderGridXBeginf, #placeholderGridYBeginf, #placeholderGridZBeginf };

      particlesParams[ #placeholderSubdomainNumber ].gridSize = { particlesParams[ #placeholderSubdomainNumber ].gridXsize, particlesParams[ #placeholderSubdomainNumber ].gridYsize, particlesParams[ #placeholderSubdomainNumber ].gridYsize };
      particlesParams[ #placeholderSubdomainNumber ].numberOfGridCells = particlesParams[ #placeholderSubdomainNumber ].gridXsize * particlesParams[ #placeholderSubdomainNumber ].gridYsize * particlesParams[ #placeholderSubdomainNumber ].gridZsize;

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

class Subdomains:

    def __init__( self, numberOfSubdomains ):
        self.numberOfSubdomains = numberOfSubdomains

        self.gridSplits = []
        self.gridSizes = []
        self.gridOrigins = []
        self.gridIndexOrigins = []

        self.subdomainConfigString = ''

    def divideIntoSubdomains( self, searchRadius, points_fluid, gridSize, gridBegin ):
        self.gridSize_x = gridSize[ 0 ]; self.gridSize_y = gridSize[ 1 ]; self.gridSize_z = gridSize[ 2 ]
        self.gridBegin_x = gridBegin[ 0 ]; self.gridBegin_y = gridBegin[ 1 ]; self.gridBegin_z = gridBegin[ 2 ]
        self.searchRadius = searchRadius
        numberOfPtcsTotal = len( points_fluid )

        self.numberOfFluidParticlesPerSubdomain = ( int )( numberOfPtcsTotal / self.numberOfSubdomains )
        for subdomain in range( self.numberOfSubdomains - 1 ):
            self.gridSplits.append( math.ceil(
                points_fluid[ self.numberOfFluidParticlesPerSubdomain * ( subdomain + 1 ) ][ 0 ] / searchRadius ) )

        for subdomain in range( self.numberOfSubdomains ):
            if subdomain == 0:
                self.gridSizes.append( self.gridSplits[ subdomain ] - 0 )
                self.gridOrigins.append( self.gridBegin_x )
                self.gridIndexOrigins.append( 0 )
            elif subdomain == self.numberOfSubdomains - 1:
                self.gridSizes.append( self.gridSize_x - self.gridSplits[ subdomain - 1 ] )
                self.gridOrigins.append( self.gridBegin_x + self.gridSplits[ subdomain - 1 ] * searchRadius )
                self.gridIndexOrigins.append( self.gridSplits[ subdomain - 1 ] )
            else:
                self.gridSizes.append( self.gridSplits[ subdomain - 1 ] - self.gridIndexOrigins[ subdomain -1 ] )
                self.gridOrigins.append( self.gridBegin_x + self.gridSplits[ subdomain - 1 ] * searchRadius )
                self.gridIndexOrigins.append( self.gridSplits[ subdomain - 1 ] )

        print( f'[Subdomians] Divide into subdomains - limits of subdomains:\n' \
               f' - number of subdomains: {self.numberOfSubdomains}\n' \
               f' - number of particles per subdoman: {self.numberOfFluidParticlesPerSubdomain}\n' \
               f' - grid split in columns: {self.gridSplits}\n' \
               f' - subdomains sizes: {self.gridSizes}\n' \
               f' - subdomains origin: {self.gridOrigins}\n' \
               f' - subdomains index origins: {self.gridIndexOrigins}' )

    def generateSubdomains( self, points_fluid, fields_fluid, points_box, fields_box ):

        inputFluidFilesString = ''
        inputBoundaryFilesString = ''

        for subdomain in range( self.numberOfSubdomains ):
            self.subdomainConfigString += self.generateSubdomain ( subdomain, points_fluid, [], points_box, [] )

            inputFluidFilesString += "\"sources/dambreak_fluid_subdomain" + str( subdomain ) + '.vtk\"'
            inputBoundaryFilesString += "\"sources/dambreak_boundary_subdomain" + str( subdomain ) + '.vtk\"'

            if subdomain < self.numberOfSubdomains - 1:
                inputFluidFilesString += ', '
                inputBoundaryFilesString += ', '

        return self.subdomainConfigString, inputFluidFilesString, inputBoundaryFilesString

    def generateSubdomain( self,subdomain, points_fluid, fields_fluid, points_box, fields_box ):
        #Fields and variables to set
        subdomain_fluid_r = []
        subdomain_fluid_ptype = []

        subdomain_box_r = []
        subdomain_box_ptype = []

        subdomain_counter_nptcs = 0
        subdomain_counter_nptcsReal = 0

        subdomain_gridSector = []

        #Load the limits of current subdomain
        if subdomain == 0:
            lowerPositionLimit = self.gridOrigins[ subdomain ]
            upperPositionLimit = self.gridOrigins[ subdomain + 1 ]

            gridIdxOverlapStart = 0
            gridIdxStart = 0
            gridIdxOverlapEnd = self.gridIndexOrigins[ subdomain + 1 ]
            gridIdxEnd = self.gridIndexOrigins[ subdomain + 1 ] - 1

        elif subdomain == self.numberOfSubdomains - 1:
            lowerPositionLimit = self.gridOrigins[ subdomain ] - self.searchRadius
            upperPositionLimit = self.gridBegin_x + ( self.gridSize_x + 1 ) * self.searchRadius

            gridIdxOverlapStart = self.gridIndexOrigins[ subdomain ] - 1
            gridIdxStart = self.gridIndexOrigins[ subdomain ]
            gridIdxOverlapEnd = self.gridSize_x
            gridIdxEnd = self.gridSize_x
        else:
            lowerPositionLimit = self.gridOrigins[ subdomain ]
            upperPositionLimit = self.gridOrigins[ subdomain + 1 ]

            gridIdxOverlapStart = self.gridIndexOrigins[ subdomain ] - 1
            gridIdxStart = self.gridIndexOrigins[ subdomain ]
            gridIdxOverlapEnd = self.gridIndexOrigins[ subdomain + 1 ]
            gridIdxEnd = self.gridIndexOrigins[ subdomain + 1 ] - 1

        #Prepare particle fields for given subdomain - fluid
        for i in range ( len( points_fluid ) ):
            if points_fluid[ i ][ 0 ] > lowerPositionLimit and points_fluid[ i ][ 0 ] <= upperPositionLimit:
                subdomain_fluid_r.append( points_fluid[ i ] )
                subdomain_fluid_ptype.append( 0 )
                subdomain_counter_nptcs += 1
                subdomain_counter_nptcsReal += 1

            if points_fluid[ i ][ 0 ] > upperPositionLimit and points_fluid[ i ][ 0 ] <= upperPositionLimit + self.searchRadius:
                subdomain_fluid_r.append( points_fluid[ i ] )
                subdomain_fluid_ptype.append( 1 )
                subdomain_counter_nptcsReal += 1

        #Prepare particle fields for given subdomain - boundary
        for i in range ( len( points_box ) ):
            if points_box[ i ][ 0 ] > lowerPositionLimit and points_box[ i ][ 0 ] <= upperPositionLimit:
                subdomain_box_r.append( points_box[ i ] )
                subdomain_box_ptype.append( 0 )

            if points_box[ i ][ 0 ] > upperPositionLimit and points_box[ i ][ 0 ] <= upperPositionLimit + self.searchRadius:
                subdomain_box_r.append( points_box[ i ] )
                subdomain_box_ptype.append( 1 )

        #TODO: Remove
        rho0 = 1000.

        import saveParticlesVTK
        #Write particle for given subdomain - fluid
        r = np.array( subdomain_fluid_r, dtype=float )# .T #!!
        v = np.zeros( ( len( subdomain_fluid_r ), 3 ) )
        rho = rho0 * np.ones( len( subdomain_fluid_r ) )
        p = np.zeros( len( subdomain_fluid_r ) )
        ptype = np.array( ( subdomain_fluid_ptype ), dtype=float ).T

        fluidToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
        subdomain_fluid_outputname = "sources/dambreak_fluid_subdomain" + str( subdomain ) + '.vtk'
        print( f'Subdomain fluid ouputfilename: {subdomain_fluid_outputname}' )
        saveParticlesVTK.save_polydata( fluidToWrite, subdomain_fluid_outputname )

        #Write particle for given subdomain - boundary
        r = np.array( subdomain_box_r, dtype=float ) # .T #!!
        v = np.zeros( ( len( subdomain_box_r ), 3 ) )
        rho = rho0 * np.ones( len( subdomain_box_r ) )
        p = np.zeros( len( subdomain_box_r ) )
        ptype = np.array( ( subdomain_box_ptype ), dtype=float ).T

        boxToWrite = saveParticlesVTK.create_pointcloud_polydata( r, v, rho, p, ptype )
        subdomain_boundary_outputname = "sources/dambreak_boundary_subdomain" + str( subdomain ) + '.vtk'
        print( f'Subdomain boundary ouputfilename: {subdomain_boundary_outputname}' )
        saveParticlesVTK.save_polydata( boxToWrite, subdomain_boundary_outputname )

        #Prepare informations about iven subdomain
        from copy import copy
        infoString = copy( subdomainStringTemplate )

        infoString = infoString.replace( '#placeholderSubdomainNumber', str( subdomain ) )

        infoString = infoString.replace( '#placeholderSearchRadius', str( self.searchRadius ) )
        infoString = infoString.replace( '#placeholderGridXSize', str( self.gridSize_x ) )
        infoString = infoString.replace( '#placeholderGridYSize', str( self.gridSize_y ) )
        infoString = infoString.replace( '#placeholderGridZSize', str( self.gridSize_z ) )
        infoString = infoString.replace( '#placeholderGridXBegin', str( round( self.gridBegin_x, 9  ) ) )
        infoString = infoString.replace( '#placeholderGridYBegin', str( round( self.gridBegin_y, 9  ) ) )
        infoString = infoString.replace( '#placeholderGridZBegin', str( round( self.gridBegin_z, 9  ) ) )

        infoString = infoString.replace( '#placeholderFluidParticles', str( len( subdomain_fluid_r ) ) )
        infoString = infoString.replace( '#placeholderAllocatedFluidParticles', str( len( subdomain_fluid_r ) * 2  ) )
        infoString = infoString.replace( '#placeholderBoundaryParticles', str( len( subdomain_box_r ) ) )
        infoString = infoString.replace( '#placeholderAllocatedBoundaryParticles', str( len( subdomain_box_r ) * 2 ) )

        infoString = infoString.replace( '#placeholderParticleIdxStart', str( 0 ) )
        infoString = infoString.replace( '#placeholderParticleIdxRealStart', str( 0 ) )
        infoString = infoString.replace( '#placeholderParticleIdxEnd', str( subdomain_counter_nptcs - 1  ) )
        infoString = infoString.replace( '#placeholderParticleIdxRealEnd', str( subdomain_counter_nptcsReal - 1) )

        infoString = infoString.replace( '#placeholderGridIdxOverlapStart', str( gridIdxOverlapStart ) )
        infoString = infoString.replace( '#placeholderGridIdxStart', str( gridIdxStart ) )
        infoString = infoString.replace( '#placeholderGridIdxOverlapEnd', str( gridIdxOverlapEnd ) )
        infoString = infoString.replace( '#placeholderGridIdxEnd', str( gridIdxEnd ) )

        #subdomainStringsArrays.append( infoString )

        #Generate grid for given Subdomain
        for y in range ( self.gridSize_y ):
            for x in range ( self.gridSizes[ subdomain ] + 1 ):

                #First subdomain:
                if subdomain == 0:
                    if x < self.gridSizes[ 0 ]:
                        subdomain_gridSector.append( 0 )
                    else:
                        subdomain_gridSector.append( 1 )

                #Last subdomain:
                elif subdomain == self.numberOfSubdomains - 1:
                    if x == 0 :
                        subdomain_gridSector.append( subdomain - 1 )
                    else:
                        subdomain_gridSector.append( subdomain )

                #Subdomain in the middle:
                else:
                    if x == 0 :
                        subdomain_gridSector.append( subdomain - 1 )
                    elif x == self.gridSizes[ subdomain ]:
                        subdomain_gridSector.append( subdomain )
                    else:
                        subdomain_gridSector.append( subdomain + 1 )

        #TODO: Add local verlap to grid begin:
        gridXOriginWithOverlap = self.gridOrigins[ subdomain ]
        if subdomain > 0: gridXOriginWithOverlap -= self.searchRadius

        # Write local grid
        subdomain_grid_outputname = "sources/dambreak_grid_subdomain" + str( subdomain ) + '.vtk'
        self.DomainGrid( self.gridSizes[ subdomain ] + 1, self.gridSize_y, self.gridSize_z,   # grid size
                         gridXOriginWithOverlap, self.gridBegin_y, self.gridBegin_z,          # coordinates of grid origin
                         subdomain_gridSector,                                                # array with index of grid sector
                         self.searchRadius,
                         subdomain_grid_outputname )                                          # outputfile name

        return infoString


    def DomainGrid( self, gridSize_x, gridSize_y, gridSize_z, gridBegin_x, gridBegin_y, gridBegin_z, gridSector, searchRadius, name ):
        from contextlib import redirect_stdout
        with open( name, 'w' ) as f:
            with redirect_stdout(f):
                print( "# vtk DataFile Version 3.0" )
                print( "vtk output" )
                print( "ASCII" )
                #print( "DATASET STRUCTURED_GRID" )
                print( "DATASET STRUCTURED_POINTS" )
                print( "DIMENSIONS ", gridSize_x + 1 , " ", gridSize_y + 1, " ", 1 )
                print( "ASPECT_RATIO ", searchRadius , " ", searchRadius , " ",  searchRadius )
                print( "ORIGIN ", gridBegin_x , " ", gridBegin_y , " ",  0  )
                print( "CELL_DATA ",  gridSize_x * gridSize_y * 1  )
                print( "SCALARS GridSector int 1 ")
                print( "LOOKUP_TABLE default" )
                for i in range( gridSize_x * gridSize_y * 1 ):
                    print( gridSector[ i ] )
