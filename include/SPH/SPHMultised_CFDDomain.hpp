#include "SPHMultised_CFDDomain.h"

namespace TNL {
namespace SPH {

template< typename Model >
void
SPHMultiset_CFDDomain< Model >::init( TNL::Config::ParameterContainer& parameters, TNL::Logger& logger )
{
   logger.writeHeader( "SPH domain initialization" );

   // compute domain properetis
   const VectorType domainOrigin = parameters.getXyz< VectorType >( "domainOrigin" );
   const VectorType domainSize = parameters.getXyz< VectorType >( "domainSize" );
   const RealType searchRadius = parameters.getParameter< RealType >( "searchRadius" );
   const IndexVectorType gridSize = TNL::ceil( ( domainSize - domainOrigin ) / searchRadius );

   // initialize main particle sets
   fluid->initialize( parameters.getParameter< int >( "numberOfParticles" ),
                      parameters.getParameter< int >( "numberOfAllocatedParticles" ),
                      searchRadius,
                      gridSize,
                      domainOrigin );

   boundary->initialize( parameters.getParameter< int >( "numberOfBoundaryParticles" ),
                         parameters.getParameter< int >( "numberOfAllocatedBoundaryParticles" ),
                         searchRadius,
                         gridSize,
                         domainOrigin );

   // initialize overlap sets
   // FIXME: Distributed cannot depend on compile parameter.
   if constexpr( Model::ModelConfigType::Distributed == true || Model::ModelConfigType::Subdomains == true )
   {
      int overlapCellsCount = 0;
      //TODO: Consider whether the overlap is in all dimensions
      if constexpr( Model::ModelConfigType::spaceDimension == 2 ) {
         overlapCellsCount = ( gridSize[ 0 ] + gridSize[ 1 ] + 4 ) * 2;
      }
      else if constexpr( Model::ModelConfigType::spaceDimension == 3 ) {
         const int xy = ( gridSize[ 0 ] + 2 ) * ( gridSize[ 1 ] + 2 );
         const int xz = ( gridSize[ 0 ] + 2 ) * ( gridSize[ 2 ] );
         const int yz = ( gridSize[ 1 ] ) * ( gridSize[ 2 ] );
         overlapCellsCount = 2 * xy + 2 * xz + 2 * yz;
      }
      int numberOfParticlesPerCell = parameters.getParameter< int >( "numberOfParticlesPerCell" );

      fluidOverlap->initialize( 0,
                                numberOfParticlesPerCell * overlapCellsCount,
                                searchRadius,
                                gridSize,
                                domainOrigin );

      boundaryOverlap->initialize( 0,
                                   numberOfParticlesPerCell * overlapCellsCount,
                                   searchRadius,
                                   gridSize,
                                   domainOrigin );
   }
}

}  //namespace SPH
}  //namespace TNL
