#include <iostream>
#include <fstream>

namespace TNL {
namespace Writers {

template< typename IndexVectorType, typename PointType, typename RealType >
void
writeBackgroundGrid( const std::string& outputFileName,
                     const IndexVectorType& gridDimensions,
                     const PointType& gridOrigin,
                     const RealType& edgeSize )
{
   std::ofstream gridFile( outputFileName );

   gridFile << "# vtk DataFile Version 3.0" << std::endl;
   gridFile << "vtk output" << std::endl;
   gridFile << "ASCII" << std::endl;
   gridFile << "DATASET STRUCTURED_POINTS" << std::endl;
   gridFile << "DIMENSIONS " << gridDimensions[ 0 ] + 1 << " " << gridDimensions[ 1 ] + 1 << " " << 1 << std::endl;
   gridFile << "ASPECT_RATIO " << edgeSize << " " << edgeSize << " " <<  edgeSize << std::endl;
   gridFile << "ORIGIN " << gridOrigin[ 0 ] << " " << gridOrigin[ 1 ] << " "<<  0  << std::endl;
   gridFile << "CELL_DATA " << gridDimensions[ 0 ] * gridDimensions[ 1 ] * 1  << std::endl;

   gridFile.close();
}

}
}

