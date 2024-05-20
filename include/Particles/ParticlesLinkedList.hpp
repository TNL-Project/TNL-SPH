#include "ParticlesLinkedList.h"

namespace TNL {
namespace ParticleSystem {

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
template< typename UseWithDomainDecomposition, std::enable_if_t< !UseWithDomainDecomposition::value, bool > Enabled >
void
ParticlesLinkedList< ParticleConfig, Device >::computeParticleCellIndices()
{
   const RealType searchRadius = this->radius;
   const PointType gridOrigin = this->gridOrigin;
   const IndexVectorType gridDimension = this->gridDimension;
   auto view_particeCellIndices = this->particleCellInidices.getView();
   const auto points_view = this->points.getConstView();

   auto indexParticles = [ = ] __cuda_callable__( GlobalIndexType i ) mutable
   {
      const IndexVectorType cellCoords = TNL::floor( ( points_view[ i ] - gridOrigin ) / searchRadius );
      view_particeCellIndices[ i ] = CellIndexer::EvaluateCellIndex( cellCoords, gridDimension );
   };
   Algorithms::parallelFor< DeviceType >( 0, this->numberOfParticles, indexParticles );
}

template< typename ParticleConfig, typename Device >
template< typename UseWithDomainDecomposition, std::enable_if_t< UseWithDomainDecomposition::value, bool > Enabled >
void
ParticlesLinkedList< ParticleConfig, Device >::computeParticleCellIndices()
{
   auto view_particeCellIndices = this->particleCellInidices.getView();
   const auto view_points = this->points.getConstView();
   const RealType searchRadius = this->radius;
   //TODO: Allow to capture protected members so we can use globalGridOrigin and globalGridDimension directly
   //FIXME: Global grid origin should be shifted aswell with search radius
   //const PointType shifOriginDueToOverlaps = this->radius;
   const PointType globalGridOrigin = this->globalGridOrigin;
   //const IndexVectorType globalGridDimension = this->globalGridDimension;
   //FIXME: Resolve the mess with gridOrigin/gridDimension and gridOriginWithOverlap/gridDimensionWithOverlap
   const PointType gridOriginWithOverlap_ = gridOrigin;
   const IndexVectorType gridDimensionWithOverlap_ = gridDimension;
   const IndexVectorType gridOriginGlobalCoords = TNL::floor( ( gridOriginWithOverlap_ - globalGridOrigin ) / searchRadius );

   auto indexParticles = [=] __cuda_callable__ ( GlobalIndexType i ) mutable
   {
      const PointType point = view_points[ i ];
      if( view_points[ i ][ 0 ] == FLT_MAX || view_points[ i ][ 1 ] == FLT_MAX ){
         view_particeCellIndices[ i ] = INT_MAX;
      }
      else{
         const IndexVectorType cellGlobalCoords = TNL::floor( ( point - globalGridOrigin ) / searchRadius );
         const IndexVectorType cellCoords = cellGlobalCoords - gridOriginGlobalCoords;
         view_particeCellIndices[ i ] = CellIndexer::EvaluateCellIndex( cellCoords, gridDimensionWithOverlap_ );
      }

   };
   Algorithms::parallelFor< DeviceType >( firstActiveParticle, lastActiveParticle + 1, indexParticles );
}

template < typename ParticleConfig, typename Device >
void
ParticlesLinkedList< ParticleConfig, Device >::removeParitclesOutOfDomain()
{
   const PointType domainOrigin = this->gridInteriorOrigin;
   //const PointType domainSize = this->gridInteriorDimension * this->radius;
   const PointType domainSize = this->interiorSize;
   auto view_points = this->points.getView();
   auto checkParticlePosition = [=] __cuda_callable__ ( int i ) mutable
   {
      if( this->isInsideDomain( view_points[ i ], domainOrigin, domainSize ) ){
         return 0;
      }
      else {
         view_points[ i ] = FLT_MAX;
         return 1;
      }
   };
   const GlobalIndexType numberOfParticlesToRemove = Algorithms::reduce< DeviceType >( this->firstActiveParticle,
                                                                                       this->lastActiveParticle + 1,
                                                                                       checkParticlePosition,
                                                                                       TNL::Plus() );
   this->setNumberOfParticlesToRemove( this->getNumberOfParticlesToRemove() + numberOfParticlesToRemove );
}

template < typename ParticleConfig, typename Device >
const typename ParticlesLinkedList< ParticleConfig, Device >::GlobalIndexType
ParticlesLinkedList< ParticleConfig, Device >::getNumberOfParticlesToRemove() const
{
   return numberOfParticlesToRemove;
}

template < typename ParticleConfig, typename Device >
void
ParticlesLinkedList< ParticleConfig, Device >::setNumberOfParticlesToRemove( GlobalIndexType removeCount )
{
   this->numberOfParticlesToRemove = removeCount;
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
