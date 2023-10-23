#include "GhostZone.h"

namespace TNL {
namespace ParticleSystem {

template< typename ParticleConfig >
template< typename CellIndexer >
void
ParticleZone< ParticleConfig >::assignCells( IndexVectorType startingPoint, IndexVectorType size, IndexVectorType gridSize )
{
   auto cellsInZone_view = this->cellsInZone.getView();

   if constexpr( ParticleConfig::spaceDimension == 2 )
   {
      auto init = [=] __cuda_callable__ ( const GlobalIndexType i ) mutable
      {
         cellsInZone_view[ i ] = CellIndexer::EvaluateCellIndex( startingPoint[ 0 ] + size[ 0 ] * i, startingPoint[ 1 ] + size[ 1 ] * i, gridSize );
      };
      Algorithms::parallelFor< DeviceType >( 0, this->numberOfCellsInZone, init );
   }

   //if constexpr( ParticleConfig::spaceDimension == 3 )
   //{
   //   auto init = [=] __cuda_callable__ ( const IndexVectorType i ) mutable
   //   {
   //      cellsInZone_view = CellIndexer::EvaluateCellIndex( startingPoint[ 0 ] + i[ 0 ], startingPoint[ 1 ] + i[ 1 ], startingPoint[ 2 ] + i[ 2 ], gridSize );
   //   };
   //   IndexVectorType begin{ 0, 0, 0 };
   //   Algorithms::parallelFor< DeviceType >( begin, size, init );
   //}
}

template< typename ParticleConfig >
template< typename Array >
void
ParticleZone< ParticleConfig >::assignCells( Array& inputCells )
{

}

template< typename ParticleConfig >
const typename ParticleZone< ParticleConfig >::IndexArrayType&
ParticleZone< ParticleConfig >::getCellsInZone() const
{
   return cellsInZone;
}

template< typename ParticleConfig >
typename ParticleZone< ParticleConfig >::IndexArrayType&
ParticleZone< ParticleConfig >::getCellsInZone()
{
   return cellsInZone;
}

template< typename ParticleConfig >
const typename ParticleZone< ParticleConfig >::IndexArrayType&
ParticleZone< ParticleConfig >::getParticlesInZone() const
{
   return particlesInZone;
}

template< typename ParticleConfig >
typename ParticleZone< ParticleConfig >::IndexArrayType&
ParticleZone< ParticleConfig >::getParticlesInZone()
{
   return particlesInZone;
}

template< typename ParticleConfig >
template< typename ParticlesPointer >
void
ParticleZone< ParticleConfig >::collectNumbersOfParticlesInCells( ParticlesPointer& particles )
{
   const auto firstLastParticle_view = particles->getCellFirstLastParticleList().getConstView();
   const auto cellsInZone_view = this->cellsInZone.getConstView();
   auto numberOfParticlesInCell_view = this->numberOfParticlesInCell.getView();
   std::cout << "flp: "  << firstLastParticle_view << std::endl;

   auto collectParticlesCounts = [=] __cuda_callable__ ( int i ) mutable
   {
      const GlobalIndexType cell = cellsInZone_view[ i ];
      const PairIndexType firstAndLastParticleInCell = firstLastParticle_view[ cell ];

      if( firstAndLastParticleInCell[ 0 ] != FLT_MAX )
      {
         numberOfParticlesInCell_view[ i ] = firstAndLastParticleInCell[ 1 ] - firstAndLastParticleInCell[ 0 ] + 1;
      }
   };
   Algorithms::parallelFor< DeviceType >( 0, this->numberOfCellsInZone, collectParticlesCounts );
}

template< typename ParticleConfig >
template< typename ParticlesPointer >
void
ParticleZone< ParticleConfig >::buildParticleList( ParticlesPointer& particles )
{
   std::cout << "Prefix field before scan: " << numberOfParticlesInCell << std::endl;

   Algorithms::inplaceExclusiveScan( this->numberOfParticlesInCell );
   this->numberOfCellsInZone = this->numberOfParticlesInCell.getElement( numberOfCellsInZone - 1 );

   std::cout << "Prefix field: " << numberOfParticlesInCell << std::endl;

   const auto firstLastCellParticle_view = particles->getCellFirstLastParticleList().getConstView();
   const auto cellsInZone_view = this->cellsInZone.getConstView();
   const auto numberOfParticlesInCell_view = this->numberOfParticlesInCell.getConstView();
   auto particlesInZone_view = this->particlesInZone.getView();

   auto collectParticles = [=] __cuda_callable__ ( int i ) mutable //TODO: This i is cell index, rename it
   {
      const GlobalIndexType cell = cellsInZone_view[ i ];
      const PairIndexType firstLastParticle = firstLastCellParticle_view[ cell ];

      if( firstLastParticle[ 0 ] != FLT_MAX )
      {
         const GlobalIndexType numberOfPtcsInThisCell = firstLastParticle[ 1 ] - firstLastParticle[ 0 ] + 1;
         const GlobalIndexType particleListCellPrefix = numberOfParticlesInCell_view[ i ];

         for( int j = 0; j < numberOfPtcsInThisCell; j++ )
         {
            particlesInZone_view[ particleListCellPrefix + j ] = firstLastParticle[ 0 ] + j;
         }
      }

   };
   Algorithms::parallelFor< DeviceType >( 0, this->numberOfCellsInZone, collectParticles );
}


template< typename ParticleConfig >
template< typename ParticlesPointer >
void
ParticleZone< ParticleConfig >::updateParticlesInZone( ParticlesPointer& particles )
{
   std::cout << "collect number of particles" << std::endl;
   this->collectNumbersOfParticlesInCells( particles );
   std::cout << "buildParticleList" << std::endl;
   this->buildParticleList( particles );
}

} // Particles
} // TNL
