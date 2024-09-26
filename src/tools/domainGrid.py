from contextlib import redirect_stdout

def DomainGrid( gridSize_x, gridSize_y, gridSize_z, gridBegin_x, gridBegin_y, gridBegin_z, gridSector, searchRadius, name ):
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

def DomainGrid3D( gridSize_x, gridSize_y, gridSize_z, gridBegin_x, gridBegin_y, gridBegin_z, gridSector, searchRadius, name ):
    with open( name, 'w' ) as f:
        with redirect_stdout(f):
            print( "# vtk DataFile Version 3.0" )
            print( "vtk output" )
            print( "ASCII" )
            #print( "DATASET STRUCTURED_GRID" )
            print( "DATASET STRUCTURED_POINTS" )
            print( "DIMENSIONS ", gridSize_x + 1 , " ", gridSize_y + 1, " ", gridSize_z + 1 )
            print( "ASPECT_RATIO ", searchRadius , " ", searchRadius , " ",  searchRadius )
            print( "ORIGIN ", gridBegin_x , " ", gridBegin_y , " ",  gridBegin_z  )
            print( "CELL_DATA ",  gridSize_x * gridSize_y * gridSize_z  )
            print( "SCALARS GridSector int 1 ")
            print( "LOOKUP_TABLE default" )
            for i in range( gridSize_x * gridSize_y * gridSize_z ):
                print( gridSector[ i ] )
