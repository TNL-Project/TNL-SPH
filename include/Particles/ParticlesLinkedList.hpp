#include "ParticlesLinkedList.h"
#include "details/details.h"

namespace TNL {
namespace ParticleSystem {

template < typename ParticleConfig, typename DeviceType >
void
ParticlesLinkedList< ParticleConfig, DeviceType >::setGridDimensions( const IndexVectorType& dimensions  )
{
   this->gridDimension = dimensions;
   if constexpr ( ParticleConfig::spaceDimension == 2 )
      firstLastCellParticle.setSize( dimensions[ 0 ] * dimensions[ 1 ] );
   if constexpr ( ParticleConfig::spaceDimension == 3 )
      firstLastCellParticle.setSize( dimensions[ 0 ] * dimensions[ 1 ] * dimensions[ 2 ] );
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
   const RealType searchRadius = this->radius;
   const PointType gridReferentialOrigin = this->gridReferentialOrigin;
   const PointType gridOriginWithOverlap = this->getGridOriginWithOverlap();
   const IndexVectorType gridDimensionWithOverlap = this->getGridDimensionsWithOverlap();
   const IndexVectorType gridOriginGlobalCoords = TNL::floor( ( gridOriginWithOverlap - gridReferentialOrigin ) / searchRadius );
   auto view_particeCellIndices = this->particleCellInidices.getView();
   const auto view_points = this->points.getConstView();

   auto indexParticles = [=] __cuda_callable__ ( GlobalIndexType i ) mutable
   {
      const PointType point = view_points[ i ];
      if( detail::isInvalid( view_points[ i ] )){
         view_particeCellIndices[ i ] = INT_MAX;
      }
      else{
         const IndexVectorType cellGlobalCoords = TNL::floor( ( point - gridReferentialOrigin ) / searchRadius );
         const IndexVectorType cellCoords = cellGlobalCoords - gridOriginGlobalCoords;
         view_particeCellIndices[ i ] = CellIndexer::EvaluateCellIndex( cellCoords, gridDimensionWithOverlap );
      }
   };
   Algorithms::parallelFor< DeviceType >( 0, this->numberOfParticles, indexParticles );
}

template < typename ParticleConfig, typename Device >
void
ParticlesLinkedList< ParticleConfig, Device >::removeParitclesOutOfDomain()
{
   //const PointType domainOrigin = this->gridInteriorOrigin;
   ////const PointType domainSize = this->gridInteriorDimension * this->radius;
   //const PointType domainSize = this->interiorSize;
   const PointType domainOrigin = this->gridOrigin;
   const PointType domainSize = this->gridDimension * this->radius;
   //const PointType domainSize = this->interiorSize;
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
   const GlobalIndexType numberOfParticlesToRemove = Algorithms::reduce< DeviceType >( 0,
                                                                                       this->numberOfParticles,
                                                                                       checkParticlePosition,
                                                                                       TNL::Plus() );
   this->setNumberOfParticlesToRemove( this->getNumberOfParticlesToRemove() + numberOfParticlesToRemove );
}

template < typename ParticleConfig, typename Device >
void
ParticlesLinkedList< ParticleConfig, Device >::sortParticles()
{
   const GlobalIndexType numberOfParticle = this->getNumberOfParticles();
   auto view_particleCellIndices = this->particleCellInidices.getView();
   auto view_map = this->sortPermutations->getView();
   this->sortPermutations->forAllElements( [] __cuda_callable__ ( int i, int& value ) { value = i; } );
   using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< Device >;
   ThrustDeviceType thrustDevice;
   thrust::sort_by_key( thrustDevice,
                        view_particleCellIndices.getArrayData(),
                        view_particleCellIndices.getArrayData() + numberOfParticle,
                        view_map.getArrayData() );
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
      view_firstLastCellParticle.setElement( view_particleCellIndex.getElement( 0 ), { 0, 0 } );
      return;
   }

   //resolve first particle
   const GlobalIndexType firstParticleIdx = 0;
   view_firstLastCellParticle.setElement( view_particleCellIndex.getElement( firstParticleIdx ),
         { firstParticleIdx, ( view_particleCellIndex.getElement( firstParticleIdx ) != view_particleCellIndex.getElement( firstParticleIdx + 1 ) ) ? firstParticleIdx : INT_MAX } ) ;

   // resolve particles [1, N-2)
   auto init = [=] __cuda_callable__ ( int i ) mutable
   {
      if( view_particleCellIndex[ i ] != view_particleCellIndex[ i-1 ] )
         view_firstLastCellParticle[  view_particleCellIndex[ i ] ][ 0 ] = i ;
      if( view_particleCellIndex[ i ] != view_particleCellIndex[ i+1 ] )
         view_firstLastCellParticle[  view_particleCellIndex[ i ] ][ 1 ] =  i ;
   };
   Algorithms::parallelFor< DeviceType >( 0, this->numberOfParticles - 1, init );

   //resolve last partile
   const GlobalIndexType lastParticleIdx = numberOfParticles - 1;
   const PairIndexType lastActiveCellContains = view_firstLastCellParticle.getElement(
         view_particleCellIndex.getElement( lastParticleIdx ) ); // N - 1
   view_firstLastCellParticle.setElement( view_particleCellIndex.getElement( lastParticleIdx ),
         { ( view_particleCellIndex.getElement( lastParticleIdx ) != view_particleCellIndex.getElement( lastParticleIdx - 1 ) ) ? lastParticleIdx : lastActiveCellContains[ 0 ], lastParticleIdx } );

}

template< typename ParticleConfig, typename Device >
void
ParticlesLinkedList< ParticleConfig, Device >::searchForNeighbors()
{
   if( this->getNumberOfParticles() == 0 )
      return;
   resetListWithIndices();
   computeParticleCellIndices();
   sortParticles();
   this->reorderParticles();
   //update number of particles - removed particles with invalid positions are shifted at the end of the array
   if( this->getNumberOfParticles() != 0 )
      this->setNumberOfParticles( this->getNumberOfParticles() - this->getNumberOfParticlesToRemove() );
   particlesToCells();
}

template < typename ParticleConfig, typename DeviceType >
void
ParticlesLinkedList< ParticleConfig, DeviceType >::writeProlog( TNL::Logger& logger ) const noexcept
{
   BaseType::writeProlog( logger );
}

} //namespace TNL
} //namespace Particles
