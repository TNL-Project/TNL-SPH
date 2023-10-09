#include "GhostZone.h"

namespace TNL {
namespace ParticleSystem {

template< typename ParticleConfig >
template< typename CellIndexer >
void
ParticleZone< ParticleConfig >::assignCells( IndexVectorType startingPoint, IndexVectorType size )
{
   auto cellsInZone_view = this->cellsInZone.getView();

   if constexpr( ParticleConfig::spaceDimension == 2 )
   {
      auto init = [=] __cuda_callable__ ( const IndexArrayType i ) mutable
      {
         cellsInZone_view = CellIndexer::EvaluateCellIndex( startingPoint[ 0 ] + i[ 0 ], startingPoint[ 1 ] + i[ 1 ] );
      };
      IndexVectorType begin{ 0, 0 };
      Algorithms::parallelFor< DeviceType >( begin, size, init );
   }

   if constexpr( ParticleConfig::spaceDimension == 3 )
   {
      auto init = [=] __cuda_callable__ ( const IndexArrayType i ) mutable
      {
         cellsInZone_view = CellIndexer::EvaluateCellIndex( startingPoint[ 0 ] + i[ 0 ], startingPoint[ 1 ] + i[ 1 ], startingPoint[ 2 ] + i[ 2 ] );
      };
      IndexVectorType begin{ 0, 0, 0 };
      Algorithms::parallelFor< DeviceType >( begin, size, init );
   }
}

template< typename ParticleConfig >
template< typename Array >
void
ParticleZone< ParticleConfig >::assignCells( Array& inputCells )
{

}

template< typename ParticleConfig >
template< typename ParticlesPointer >
void
ParticleZone< ParticleConfig >::collectNumbersOfParticlesInCells( ParticlesPointer& particles )
{
   const auto firstLastParticle_view = this->getFirstLastParticle().getConstView();
   const auto cellsInZone_view = this->cellsInZone.getConstView();
   auto numberOfParticlesInCell_view = this->numberOfParticlesInCell.getView();

   auto collectParticlesCounts = [=] __cuda_callable__ ( int i ) mutable
   {
      const GlobalIndexType cell = cellsInZone[ i ];
      const PairIndexType firstAndLastParticleInCell = firstLastParticle_view[ cell ];
      numberOfParticlesInCell_view[ i ] = firstAndLastParticleInCell[ 1 ] - firstAndLastParticleInCell[ 0 ];
   };
   Algorithms::parallelFor< DeviceType >( 0, this->numberOfCellsInZone, collectParticlesCounts );
}

template< typename ParticleConfig >
template< typename ParticlesPointer >
void
ParticleZone< ParticleConfig >::buildParticleList( ParticlesPointer& particles )
{
   Algorithms::inplaceInclusiveScan( this->numberOfParticlesInCell );
   this->numberOfCellsInZone = this->numberOfParticlesInCell.getElement( numberOfParticlesInZone - 1 );

   const auto firstLastParticle_view = this->getFirstLastParticle().getConstView();
   auto particlesInZone_view = this->particlesInZone.getView();

   auto collectParticles = [=] __cuda_callable__ ( int i ) mutable
   {
      const GlobalIndexType cell = cellsInZone[ i ];
      const PairIndexType firstLastParticle = firstLastCellParticle_view[ cell ];
      const GlobalIndexType numberOfPtcsInThisCell = firstLastParticle[ 1 ] - firstLastParticle[ 0 ];

      for( int i = 0; i < numberOfParticlesInCell; i++ )
      {
         particlesInZone_view[ cell + i ] = firstLastParticle[ 0 ] + i;
      }

   };
   Algorithms::parallelFor< DeviceType >( 0, this->numberOfCellsInZone, collectParticles );
}


template< typename ParticleConfig >
template< typename ParticlesPointer >
void
ParticleZone< ParticleConfig >::updateParticlesInZone( ParticlesPointer& particles )
{
   this->collectNumbersOfParticlesInCells( particles );
   this->buildParticleList( particles );
}

} // Particles
} // TNL
