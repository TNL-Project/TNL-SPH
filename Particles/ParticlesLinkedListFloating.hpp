#include "ParticlesLinkedListFloating.h"

namespace TNL {
namespace ParticleSystem {

template < typename ParticleConfig, typename DeviceType >
const typename ParticlesLinkedList< ParticleConfig, DeviceType >::GlobalIndexType
ParticlesLinkedList< ParticleConfig, DeviceType >::getFirstActiveParticle() const
{
   return this->firstActiveParticle;
}

template < typename ParticleConfig, typename DeviceType >
const typename ParticlesLinkedList< ParticleConfig, DeviceType >::GlobalIndexType
ParticlesLinkedList< ParticleConfig, DeviceType >::getLastActiveParticle() const
{
   return this->lastActiveParticle;
}

template < typename ParticleConfig, typename DeviceType >
void
ParticlesLinkedList< ParticleConfig, DeviceType >::setFirstActiveParticle( GlobalIndexType firstActiveParticle )
{
   this->firstActiveParticle = firstActiveParticle;
}

template < typename ParticleConfig, typename DeviceType >
void
ParticlesLinkedList< ParticleConfig, DeviceType >::setLastActiveParticle( GlobalIndexType lastActiveParticle )
{
   this->lastActiveParticle = lastActiveParticle;
}

template < typename ParticleConfig, typename DeviceType >
const typename ParticlesLinkedList< ParticleConfig, DeviceType >::IndexVectorType
ParticlesLinkedList< ParticleConfig, DeviceType >::getGridSize() const
{
   return gridDimension;
}

template < typename ParticleConfig, typename DeviceType >
void
ParticlesLinkedList< ParticleConfig, DeviceType >::setGridSize( IndexVectorType gridSize )
{
   gridDimension = gridSize;
}

template < typename ParticleConfig, typename DeviceType >
const typename ParticlesLinkedList< ParticleConfig, DeviceType >::PointType
ParticlesLinkedList< ParticleConfig, DeviceType >::getGridOrigin() const
{
   return gridOrigin;
}

template < typename ParticleConfig, typename DeviceType >
void
ParticlesLinkedList< ParticleConfig, DeviceType >::setGridOrigin( PointType gridBegin )
{
   gridOrigin = gridBegin;
}

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
   TNL_ASSERT_LT( particleIndex, this->numberOfParticles, "invalid particle index" );
   return this->particleCellInidices[ particleIndex ];
}

template < typename ParticleConfig, typename Device >
__cuda_callable__
typename ParticlesLinkedList< ParticleConfig, Device >::CellIndexType&
ParticlesLinkedList< ParticleConfig, Device >::getParticleCellIndex( GlobalIndexType particleIndex )
{
   TNL_ASSERT_GE( particleIndex, 0, "invalid particle index" );
   TNL_ASSERT_LT( particleIndex, this->numberOfParticles, "invalid particle index" );
   return this->particleCellInidices[ particleIndex ];
}

template < typename ParticleConfig, typename Device >
void
ParticlesLinkedList< ParticleConfig, Device >::computeParticleCellIndices()
{
   GlobalIndexType _numberOfParticles = this->numberOfParticles;

   auto view = this->particleCellInidices.getView();
   auto view_points = this->points.getView();

   CellIndexer::ComputeParticleCellIndex(
         view, view_points, firstActiveParticle, lastActiveParticle, gridDimension, gridOrigin, this->radius );
}

template < typename ParticleConfig, typename Device >
void
ParticlesLinkedList< ParticleConfig, Device >::sortParticles()
{

   GlobalIndexType numberOfParticle = this->getNumberOfParticles();
   auto view_particleCellIndices = this->particleCellInidices.getView();
   auto view_map = this->sortPermutations->getView();
   this->sortPermutations->forAllElements( [] __cuda_callable__ ( int i, int& value ) { value = i; } );
   thrust::sort_by_key( thrust::device,
                        view_particleCellIndices.getArrayData() + firstActiveParticle,
                        view_particleCellIndices.getArrayData() + lastActiveParticle + 1,
                        view_map.getArrayData() ); //TODO: replace thrust::device

   auto view_points = this->getPoints().getView();
   auto view_points_swap = this->points_swap.getView();
   thrust::gather( thrust::device,
                   view_map.getArrayData() + firstActiveParticle,
                   view_map.getArrayData() + lastActiveParticle + 1,
                   view_points.getArrayData(),
                   view_points_swap.getArrayData() ); //TODO: replace thrust::device
   this->getPoints().swap( this->points_swap );
}

template< typename ParticleConfig, typename Device >
void
ParticlesLinkedList< ParticleConfig, Device >::resetListWithIndices()
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
ParticlesLinkedList< ParticleConfig, Device >::particlesToCells()
{
   const GlobalIndexType numberOfParticles = this->getNumberOfParticles();

   if( numberOfParticles == 0 ) //temp
      return;

   auto view_firstLastCellParticle = this->firstLastCellParticle.getView();
   const auto view_particleCellIndex = this->particleCellInidices.getView();

   if( numberOfParticles == 1 ) //temp
   {
      view_firstLastCellParticle.setElement( view_particleCellIndex.getElement( firstActiveParticle ), { 0, 0 } );
      return;
   }

   //resolve first particle
   view_firstLastCellParticle.setElement( view_particleCellIndex.getElement( firstActiveParticle ),
         { 0, ( view_particleCellIndex.getElement( firstActiveParticle ) != view_particleCellIndex.getElement( firstActiveParticle + 1 ) ) ? 0 : INT_MAX } ) ;

   auto init = [=] __cuda_callable__ ( int i ) mutable
   {
      if( view_particleCellIndex[ i ] != view_particleCellIndex[ i-1 ] )
         view_firstLastCellParticle[  view_particleCellIndex[ i ] ][ 0 ] = i ;
      if( view_particleCellIndex[ i ] != view_particleCellIndex[ i+1 ] )
         view_firstLastCellParticle[  view_particleCellIndex[ i ] ][ 1 ] =  i ;
   };
   Algorithms::parallelFor< DeviceType >( firstActiveParticle, lastActiveParticle, init ); // [0, N-1)

   //resolve last partile
   //I think there is bug in the initial version. In case there are two particles in the last cell, the first particle in last cell is overwritten.
   /*
   view_firstLastCellParticle.setElement( view_particleCellIndex.getElement( numberOfParticles - 1 ),
         { ( view_particleCellIndex.getElement( numberOfParticles -1 ) != view_particleCellIndex.getElement( numberOfParticles-2 ) ) ? numberOfParticles-1 : INT_MAX, numberOfParticles - 1 } );
   */
   //Workaround
   PairIndexType lastActiveCellContains = view_firstLastCellParticle.getElement( view_particleCellIndex.getElement( lastActiveParticle ) ); // N - 1
   view_firstLastCellParticle.setElement( view_particleCellIndex.getElement( lastActiveParticle ),
         { ( view_particleCellIndex.getElement( lastActiveParticle ) != view_particleCellIndex.getElement( lastActiveParticle - 1 ) ) ? lastActiveParticle : lastActiveCellContains[ 0 ], lastActiveParticle } );

}

//move to detail
template< typename ParticleConfig, typename Device >
typename ParticlesLinkedList< ParticleConfig, Device >::PairIndexType
ParticlesLinkedList< ParticleConfig, Device >::getFirstLastParticleInColumnOfCells( const GlobalIndexType& gridColumn )
{
   //static_assert( std::is_same< CellIndexer::, DeviceType >::value, "mismatched DeviceType of the array" );

   const GlobalIndexType indexOfFirstColumnCell = CellIndexer::EvaluateCellIndex( gridColumn, 1, gridDimension );
   const GlobalIndexType indexOfLastColumnCell = CellIndexer::EvaluateCellIndex(
         gridColumn, gridDimension[ 1 ] - 1, gridDimension );
   const auto view_firstLastCellParticle = firstLastCellParticle.getConstView( indexOfFirstColumnCell, indexOfLastColumnCell );

   auto fetch_vect = [=] __cuda_callable__ ( int i ) -> PairIndexType  { return view_firstLastCellParticle[ i ]; };
   auto reduction_vect = [=] __cuda_callable__ ( const PairIndexType& a, const PairIndexType& b ) -> PairIndexType
   { return { min( a[ 0 ], b[ 0 ] ), max( a[ 1 ], ( b[ 1 ] < INT_MAX ) ? b[ 1 ] : -1 ) }; };

   PairIndexType identity = { INT_MAX , INT_MIN };
   PairIndexType firstLastParticle = Algorithms::reduce< Devices::Cuda >(
         0, view_firstLastCellParticle.getSize(), fetch_vect, reduction_vect, identity );

   return firstLastParticle;
}



} //namespace TNL
} //namespace Particles
