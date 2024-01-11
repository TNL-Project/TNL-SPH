#include "ParticlesLinkedList.h"

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
   if constexpr ( ParticleConfig::spaceDimension == 2 )
      firstLastCellParticle.setSize( gridSize[ 0 ] * gridSize[ 1 ] );
   if constexpr ( ParticleConfig::spaceDimension == 3 )
      firstLastCellParticle.setSize( gridSize[ 0 ] * gridSize[ 1 ] * gridSize[ 2 ] );
}

template < typename ParticleConfig, typename DeviceType >
void
ParticlesLinkedList< ParticleConfig, DeviceType >::setSize( const GlobalIndexType& size )
{
   BaseType::setSize( size );
   this->particleCellInidices.setSize( size );
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

template< typename ParticleConfig, typename Device >
const typename ParticlesLinkedList< ParticleConfig, Device >::CellIndexArrayType&
ParticlesLinkedList< ParticleConfig, Device >::getParticleCellIndices() const
{
   return particleCellInidices;
}

template< typename ParticleConfig, typename Device >
typename ParticlesLinkedList< ParticleConfig, Device >::CellIndexArrayType&
ParticlesLinkedList< ParticleConfig, Device >::getParticleCellIndices()
{
   return particleCellInidices;
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
   using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< Device >;
   ThrustDeviceType thrustDevice;
   thrust::sort_by_key( thrustDevice,
                        view_particleCellIndices.getArrayData() + firstActiveParticle,
                        view_particleCellIndices.getArrayData() + lastActiveParticle + 1,
                        view_map.getArrayData() ); //TODO: replace thrust::device

   auto view_points = this->getPoints().getView();
   auto view_points_swap = this->points_swap.getView();
   thrust::gather( thrustDevice,
                   view_map.getArrayData(),
                   view_map.getArrayData() + numberOfParticle,
                   view_points.getArrayData() + firstActiveParticle,
                   view_points_swap.getArrayData() + firstActiveParticle ); //TODO: replace thrust::device
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
         { firstActiveParticle, ( view_particleCellIndex.getElement( firstActiveParticle ) != view_particleCellIndex.getElement( firstActiveParticle + 1 ) ) ? firstActiveParticle : INT_MAX } ) ; //careful with the firstActiveParticle instead of 0

   auto init = [=] __cuda_callable__ ( int i ) mutable
   {
      if( view_particleCellIndex[ i ] != view_particleCellIndex[ i-1 ] )
         view_firstLastCellParticle[  view_particleCellIndex[ i ] ][ 0 ] = i ;
      if( view_particleCellIndex[ i ] != view_particleCellIndex[ i+1 ] )
         view_firstLastCellParticle[  view_particleCellIndex[ i ] ][ 1 ] =  i ;
   };
   Algorithms::parallelFor< DeviceType >( firstActiveParticle + 1, lastActiveParticle, init ); // [1, N-1)

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

//move to detail
template< typename ParticleConfig, typename Device >
typename ParticlesLinkedList< ParticleConfig, Device >::PairIndexType
ParticlesLinkedList< ParticleConfig, Device >::getFirstLastParticleInBlockOfCells( const GlobalIndexType& gridBlock )
{
   PairIndexType firstLastParticle;

   //for( int j = 1; j < gridDimension[ 1 ]; j++ )
   for( int j = gridDimension[ 1 ]; j > 0; j-- )
   {
      const GlobalIndexType indexOfFirstColumnCell = CellIndexer::EvaluateCellIndex( gridBlock, j, 1, gridDimension );
      const GlobalIndexType indexOfLastColumnCell = CellIndexer::EvaluateCellIndex(
            gridBlock, j, gridDimension[ 1 ] - 1, gridDimension );
      const auto view_firstLastCellParticle = firstLastCellParticle.getConstView( indexOfFirstColumnCell, indexOfLastColumnCell );

      auto fetch_vect = [=] __cuda_callable__ ( int i ) -> PairIndexType  { return view_firstLastCellParticle[ i ]; };
      auto reduction_vect = [=] __cuda_callable__ ( const PairIndexType& a, const PairIndexType& b ) -> PairIndexType
      { return { min( a[ 0 ], b[ 0 ] ), max( a[ 1 ], ( b[ 1 ] < INT_MAX ) ? b[ 1 ] : -1 ) }; };

      PairIndexType identity = { INT_MAX , INT_MIN };
      PairIndexType firstLastParticleLocal = Algorithms::reduce< Devices::Cuda >(
            0, view_firstLastCellParticle.getSize(), fetch_vect, reduction_vect, identity );

      //if( firstLastParticleLocal[ 0 ] < INT_MAX ){
      //   firstLastParticle[ 0 ] = firstLastParticleLocal[ 0 ];
      //   std::cout << "[ Particles::getFirstLastParticleInBlockOfCells ] [ Rank: " << TNL::MPI::GetRank() << " ] firstLastParticle (loop for first element): " << firstLastParticle << std::endl;
      //   break;
      //}
      if( firstLastParticleLocal[ 0 ] < INT_MAX ){
         firstLastParticle[ 1 ] = firstLastParticleLocal[ 1 ];
         //std::cout << "[ Particles::getFirstLastParticleInBlockOfCells ] [ Rank: " << TNL::MPI::GetRank() << " ] firstLastParticle (loop for first element): " << firstLastParticle << std::endl;
         break;
      }
   }

   //for( int j = gridDimension[ 1 ]; j > 1; j-- )
   for( int j = 1; j < gridDimension[ 1 ]; j++ )
   {
      const GlobalIndexType indexOfFirstColumnCell = CellIndexer::EvaluateCellIndex( gridBlock, j, 1, gridDimension );
      const GlobalIndexType indexOfLastColumnCell = CellIndexer::EvaluateCellIndex(
            gridBlock, j, gridDimension[ 1 ] - 1, gridDimension );
      const auto view_firstLastCellParticle = firstLastCellParticle.getConstView( indexOfFirstColumnCell, indexOfLastColumnCell );

      auto fetch_vect = [=] __cuda_callable__ ( int i ) -> PairIndexType  { return view_firstLastCellParticle[ i ]; };
      auto reduction_vect = [=] __cuda_callable__ ( const PairIndexType& a, const PairIndexType& b ) -> PairIndexType
      { return { min( a[ 0 ], b[ 0 ] ), max( a[ 1 ], ( b[ 1 ] < INT_MAX ) ? b[ 1 ] : -1 ) }; };

      PairIndexType identity = { INT_MAX , INT_MIN };
      PairIndexType firstLastParticleLocal = Algorithms::reduce< Devices::Cuda >(
            0, view_firstLastCellParticle.getSize(), fetch_vect, reduction_vect, identity );

      //if( firstLastParticleLocal[ 1 ] > -1 ){
      //   firstLastParticle[ 1 ] = firstLastParticleLocal[ 1 ];
      //   std::cout << "[ Particles::getFirstLastParticleInBlockOfCells ] [ Rank: " << TNL::MPI::GetRank() << " ] firstLastParticle (loop for second element): " << firstLastParticle << std::endl;
      //   break;
      //}
      if( firstLastParticleLocal[ 1 ] > -1 ){
         firstLastParticle[ 0 ] = firstLastParticleLocal[ 0 ];
         //std::cout << "[ Particles::getFirstLastParticleInBlockOfCells ] [ Rank: " << TNL::MPI::GetRank() << " ] firstLastParticle (loop for second element): " << firstLastParticle << std::endl;
         break;
      }
   }

   //std::cout << "[ Particles::getFirstLastParticleInBlockOfCells ] [ Rank: " << TNL::MPI::GetRank() << " ] firstLastParticle (to return): " << firstLastParticle << std::endl;
   return firstLastParticle;
}

template < typename ParticleConfig, typename DeviceType >
void
ParticlesLinkedList< ParticleConfig, DeviceType >::writeProlog( TNL::Logger& logger ) const noexcept
{
   BaseType::writeProlog( logger );
   logger.writeParameter( "Grid dimensions:", this->gridDimension );
   logger.writeParameter( "Number of grid cells:", this->firstLastCellParticle.getSize() );
   logger.writeParameter( "Grid origin:", this->gridOrigin );
   logger.writeParameter( "First active particle index:", this->firstActiveParticle );
   logger.writeParameter( "Last active particle index:", this->lastActiveParticle );
}

} //namespace TNL
} //namespace Particles
