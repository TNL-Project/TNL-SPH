#include "ParticlesLinkedList.h"

namespace TNL {
namespace ParticleSystem {

template< typename ParticleConfig, typename Device >
const typename ParticlesLinkedList< ParticleConfig, Device >::PairIndexArrayType&
ParticlesLinkedList< ParticleConfig, Device >::getCellFirstLastParticleList() const
{
   return firstLastCellParticle;
}

template< typename ParticleConfig, typename Device >
typename ParticlesLinkedList< ParticleConfig, Device >::PairIndexArrayType&
ParticlesLinkedList< ParticleConfig, Device >::getCellFirstLastParticleList()
{
   return firstLastCellParticle;
}

template < typename ParticleConfig, typename Device >
__cuda_callable__
const typename ParticlesLinkedList< ParticleConfig, Device >::CellIndexType&
ParticlesLinkedList< ParticleConfig, Device >::getParticleCellIndex( GlobalIndexType particleIndex ) const
{
   TNL_ASSERT_GE( particleIndex, 0, "invalid particle index" );
   TNL_ASSERT_LT( particleIndex, numberOfParticles, "invalid particle index" );
   return this->particleCellInidices[ particleIndex ];
}

template < typename ParticleConfig, typename Device >
__cuda_callable__
typename ParticlesLinkedList< ParticleConfig, Device >::CellIndexType&
ParticlesLinkedList< ParticleConfig, Device >::getParticleCellIndex( GlobalIndexType particleIndex )
{
   TNL_ASSERT_GE( particleIndex, 0, "invalid particle index" );
   TNL_ASSERT_LT( particleIndex, numberOfParticles, "invalid particle index" );
   return this->particleCellInidices[ particleIndex ];
}

template < typename ParticleConfig, typename Device >
void
ParticlesLinkedList< ParticleConfig, Device >::computeParticleCellIndices()
{
   GlobalIndexType _numberOfParticles = this->numberOfParticles;

   auto view = this->particleCellInidices.getView();
   auto view_points = this->points.getView();

   CellIndexer::ComputeParticleCellIndex( view, view_points, _numberOfParticles, this->gridDimension, this->gridOrigin, this->radius );
}

template < typename ParticleConfig, typename Device >
void
ParticlesLinkedList< ParticleConfig, Device >::sortParticles()
{

   GlobalIndexType numberOfParticle = this->getNumberOfParticles();
   auto view_particleCellIndices = this->particleCellInidices.getView();
   auto view_map = this->sortPermutations->getView();

   this->sortPermutations->forAllElements( [] __cuda_callable__ ( int i, int& value ) { value = i; } );
   thrust::sort_by_key( thrust::device, view_particleCellIndices.getArrayData(),
         view_particleCellIndices.getArrayData() + numberOfParticle, view_map.getArrayData() ); //TODO: replace thrust::device

   auto view_points = this->getPoints().getView();
   auto view_points_swap = this->points_swap.getView();
   thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticle,
         view_points.getArrayData(), view_points_swap.getArrayData() );
   this->getPoints().swap( this->points_swap );
}

//neighborSearch utilities
template< typename ParticleConfig, typename Device >
void
ParticlesLinkedList< ParticleConfig, Device >::resetListWithIndices
()
{
   auto view_firstLastCellParticle = this->firstLastCellParticle.getView();
   auto init = [=] __cuda_callable__ ( int i ) mutable
   {
      view_firstLastCellParticle[ i ] = INT_MAX ;
   };
   Algorithms::parallelFor< DeviceType >( 0, this->firstLastCellParticle.getSize(), init );
}

template< typename ParticleConfig, typename Device >
void
ParticlesLinkedList< ParticleConfig, Device >::particlesToCells
()
{
   const GlobalIndexType numberOfParticles = this->getNumberOfParticles();

   if( numberOfParticles == 0 ) //temp
      return;

   auto view_firstLastCellParticle = this->firstLastCellParticle.getView();
   const auto view_particleCellIndex = this->particleCellInidices.getView();

   if( numberOfParticles == 1 ) //temp
   {
      view_firstLastCellParticle.setElement( view_particleCellIndex.getElement( 0 ), { 0, 0 } );
      return;
   }

   //resolve first particle
   view_firstLastCellParticle.setElement( view_particleCellIndex.getElement( 0 ),
         { 0, ( view_particleCellIndex.getElement( 0 ) != view_particleCellIndex.getElement( 0+1 ) ) ? 0 : INT_MAX } ) ;

   auto init = [=] __cuda_callable__ ( int i ) mutable
   {
      if( view_particleCellIndex[ i ] != view_particleCellIndex[ i-1 ] )
         view_firstLastCellParticle[  view_particleCellIndex[ i ] ][ 0 ] = i ;
      if( view_particleCellIndex[ i ] != view_particleCellIndex[ i+1 ] )
         view_firstLastCellParticle[  view_particleCellIndex[ i ] ][ 1 ] =  i ;
   };
   Algorithms::parallelFor< DeviceType >( 1, numberOfParticles - 1, init );

   //resolve last partile
   //I think there is bug in the initial version. In case there are two particles in the last cell, the first particle in last cell is overwritten.
   /*
   view_firstLastCellParticle.setElement( view_particleCellIndex.getElement( numberOfParticles - 1 ),
         { ( view_particleCellIndex.getElement( numberOfParticles -1 ) != view_particleCellIndex.getElement( numberOfParticles-2 ) ) ? numberOfParticles-1 : INT_MAX, numberOfParticles - 1 } );
   */
   //Workaround
   PairIndexType lastActiveCellContains = view_firstLastCellParticle.getElement( view_particleCellIndex.getElement( numberOfParticles - 1 ) );
   view_firstLastCellParticle.setElement( view_particleCellIndex.getElement( numberOfParticles - 1 ),
         { ( view_particleCellIndex.getElement( numberOfParticles -1 ) != view_particleCellIndex.getElement( numberOfParticles-2 ) ) ? numberOfParticles-1 : lastActiveCellContains[ 0 ], numberOfParticles - 1 } );


}



} //namespace TNL
} //namespace Particles
