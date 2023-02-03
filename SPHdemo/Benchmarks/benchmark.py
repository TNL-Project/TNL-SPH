
resolution = 0.001

with open( 'dualSPHysics_resources/damBreak2D_WCSPH-DBC_template.xml', 'r' ) as file :
  file_dualSPHysics_Conf = file.read()

file_dualSPHysics_Conf = file_dualSPHysics_Conf.replace( 'resolutionPlaceholder', str( resolution ) )

with open( 'dualSPHysics_resources/damBreak2D_WCSPH-DBC.xml', 'w' ) as file:
  file.write( file_dualSPHysics_Conf )
