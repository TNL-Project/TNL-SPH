import math

def getDomainLimits( points_fluid, points_box, searchRadius, modifySize_x, modifySize_y, modifySize_z ):

    gridBegin_x = 1.005 * ( min( min( points_fluid[ : , 0 ] ), min( points_box[ :, 0 ] ) ) - searchRadius )
    gridBegin_y = 1.005 * ( min( min( points_fluid[ : , 1 ] ), min( points_box[ :, 1 ] ) ) - searchRadius )
    gridBegin_z = 1.005 * ( min( min( points_fluid[ : , 2 ] ), min( points_box[ : ,2 ] ) ) - searchRadius )

    gridEnd_x = 1.005 * ( max( max( points_fluid[ :, 0 ] ), max( points_box[ :, 0 ] ) ) + searchRadius )
    gridEnd_y = 1.005 * ( max( max( points_fluid[ :, 1 ] ), max( points_box[ :, 1 ] ) ) + searchRadius )
    gridEnd_z = 1.005 * ( max( max( points_fluid[ :, 2 ] ), max( points_box[ :, 2 ] ) ) + searchRadius )

    gridSize_x = math.ceil( modifySize_x * ( gridEnd_x - gridBegin_x ) / searchRadius )
    gridSize_y = math.ceil( modifySize_y * ( gridEnd_y - gridBegin_y ) / searchRadius )
    gridSize_z = math.ceil( modifySize_z * ( gridEnd_z - gridBegin_z ) / searchRadius )

    gridBegin = [ gridBegin_x, gridBegin_y, gridBegin_z ]
    gridEnd = [ gridEnd_x, gridEnd_y, gridEnd_z ]
    gridSize = [ gridSize_x, gridSize_y, gridSize_z ]

    print( f'[getDomainLimits] Size and limits of domain:\n - gridBegin: {gridBegin}\n - gridEnd: {gridEnd}\n - gridSize: {gridSize}' )

    return gridBegin, gridEnd, gridSize
