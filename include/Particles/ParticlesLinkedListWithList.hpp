#include "ParticlesLinkedListWithList.h"
#include "ParticlesTraits.h"

namespace TNL {
namespace ParticleSystem {

template< typename ParticleConfig, typename DeviceType >
void
ParticlesLinkedListWithList< ParticleConfig, DeviceType >::setSize( const GlobalIndexType& size )
{
   BaseType::setSize( size );
   this->neighborListStorage.setSize( size * this->neighborsCountLimit );
}

template< typename ParticleConfig, typename DeviceType >
void
ParticlesLinkedListWithList< ParticleConfig, DeviceType >::setNeighborsCountLimit( const GlobalIndexType& limit )
{
   this->neighborsCountLimit = limit;
   this->neighborListStorage.setSize( this->getNumberOfAllocatedParticles() * limit );
}

template< typename ParticleConfig, typename DeviceType >
const typename ParticlesLinkedListWithList< ParticleConfig, DeviceType >::GlobalIndexType
ParticlesLinkedListWithList< ParticleConfig, DeviceType >::getNeighborsCountLimit() const
{
   return neighborsCountLimit;
}

template< typename ParticleConfig, typename DeviceType >
void
ParticlesLinkedListWithList< ParticleConfig, DeviceType >::setParticleSetLabel( const int& label )
{
   this->particleSetLabel = label;
}

template< typename ParticleConfig, typename DeviceType >
const int
ParticlesLinkedListWithList< ParticleConfig, DeviceType >::getParticleSetLabel() const
{
   return particleSetLabel;
}

template< typename ParticleConfig, typename DeviceType >
const typename ParticlesLinkedListWithList< ParticleConfig, DeviceType >::NeighborListType&
ParticlesLinkedListWithList< ParticleConfig, DeviceType >::getNeighborList() const
{
   return neighborList;
}

template< typename ParticleConfig, typename DeviceType >
typename ParticlesLinkedListWithList< ParticleConfig, DeviceType >::NeighborListType&
ParticlesLinkedListWithList< ParticleConfig, DeviceType >::getNeighborList()
{
   return neighborList;
}

template< typename ParticleConfig, typename DeviceType >
const typename ParticlesLinkedListWithList< ParticleConfig, DeviceType >::IndexArrayType&
ParticlesLinkedListWithList< ParticleConfig, DeviceType >::getNeighborListStorage() const
{
   return neighborListStorage;
}

template < typename ParticleConfig, typename DeviceType >
typename ParticlesLinkedListWithList< ParticleConfig, DeviceType >::IndexArrayType&
ParticlesLinkedListWithList< ParticleConfig, DeviceType >::getNeighborListStorage()
{
   return neighborListStorage;
}

template< typename ParticleConfig, typename Device >
typename ParticlesLinkedListWithList< ParticleConfig, Device >::NeighborsLoopParams
ParticlesLinkedListWithList< ParticleConfig, Device >::getCLLwLSearchToken()
{
   NeighborsLoopParams searchToken;

   searchToken.particleSetLabel = this->getParticleSetLabel();
   searchToken.numberOfParticles = this->getNumberOfParticles();
   searchToken.neighborListStorageView.bind( this->getNeighborListStorage().getView() );
   searchToken.neighborsCountLimit = this->getNeighborsCountLimit();

   return searchToken;
}

template< typename ParticleConfig, typename Device >
template< typename ParticlesPointerType >
typename ParticlesLinkedListWithList< ParticleConfig, Device >::NeighborsLoopParams
ParticlesLinkedListWithList< ParticleConfig, Device >::getCLLwLSearchToken( ParticlesPointerType& particlesToSearch )
{
   NeighborsLoopParams searchToken;

   searchToken.particleSetLabel = particlesToSearch->getParticleSetLabel();
   searchToken.numberOfParticles = this->getNumberOfParticles();
   searchToken.neighborListStorageView.bind( this->getNeighborListStorage().getView() );
   searchToken.neighborsCountLimit = this->getNeighborsCountLimit();

   return searchToken;
}

template< typename ParticleConfig, typename Device >
typename ParticlesLinkedListWithList< ParticleConfig, Device >::NeighborsLoopParams
ParticlesLinkedListWithList< ParticleConfig, Device >::getSearchToken()
{
   return this->getCLLwLSearchToken();
}

template< typename ParticleConfig, typename Device >
template< typename ParticlesPointerType >
typename ParticlesLinkedListWithList< ParticleConfig, Device >::NeighborsLoopParams
ParticlesLinkedListWithList< ParticleConfig, Device >::getSearchToken( ParticlesPointerType& particlesToSearch )
{
   return this->getCLLwLSearchToken( particlesToSearch );
}

template < typename ParticleConfig, typename DeviceType >
void
ParticlesLinkedListWithList< ParticleConfig, DeviceType >::resetNeighborList()
{
   this->neighborListStorage = 0;
}

template < typename ParticleConfig, typename DeviceType >
void
ParticlesLinkedListWithList< ParticleConfig, DeviceType >::makeSetSearchable()
{
   BaseType::searchForNeighbors();
}

template< typename ParticleConfig, typename DeviceType >
template< typename Layout,
          std::enable_if_t< std::is_same_v< Layout, NeighborListLayouts::NeighborMajorLinear >, bool > Enabled >
void
ParticlesLinkedListWithList< ParticleConfig, DeviceType >::buildParticleList()
{
   auto searchInSet = this->getCLLSearchToken();

   const RealType searchRadius = this->getSearchRadius();
   const auto view_points = this->getPoints().getConstView();
   auto neighborListStorage_view = this->neighborListStorage.getView();
   const GlobalIndexType neighborsCountLimit = this->neighborsCountLimit;
   const GlobalIndexType numberOfParticles = this->getNumberOfParticles();

   auto compareDistance = [=] __cuda_callable__ ( LocalIndexType i,
                                                  LocalIndexType j,
                                                  PointType& r_i,
                                                  GlobalIndexType& globIdxStart,
                                                  GlobalIndexType* neihgborsCount ) mutable
   {
      const PointType r_j = view_points[ j ];
      const PointType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         neighborListStorage_view[ numberOfParticles * ( *neihgborsCount + 1 ) + i ] = j;
         ( *neihgborsCount )++;
      }
   };

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const PointType r_i = view_points[ i ];
      const GlobalIndexType globIdxStart = i * neighborsCountLimit;
      GlobalIndexType neihgborsCount = 0;

      BaseType::NeighborsLoop::exec( i, r_i, searchInSet, compareDistance, globIdxStart, &neihgborsCount );
      neighborListStorage_view[ i ] = neihgborsCount;

   };
   this->forAll( particleLoop );
}

template< typename ParticleConfig, typename DeviceType >
template< typename ParticlesPointerType,
          typename Layout,
          std::enable_if_t< std::is_same_v< Layout, NeighborListLayouts::NeighborMajorLinear >, bool > Enabled >
void
ParticlesLinkedListWithList< ParticleConfig, DeviceType >::addToParticleList( ParticlesPointerType& particlesToSearch )
{
   auto searchInSet = particlesToSearch->getCLLSearchToken();

   const RealType searchRadius = this->getSearchRadius();
   const auto view_points = this->getPoints().getConstView();
   const auto view_pointsToSearch = particlesToSearch->getPoints().getConstView();
   auto neighborListStorage_view = this->neighborListStorage.getView();
   const GlobalIndexType neighborsCountLimit = this->neighborsCountLimit;
   const int particlesToSearchLabel = particlesToSearch->getParticleSetLabel();
   const GlobalIndexType numberOfParticles = this->getNumberOfParticles();

   auto compareDistance = [=] __cuda_callable__ ( LocalIndexType i,
                                                  LocalIndexType j,
                                                  PointType& r_i,
                                                  int& anotherSetOffset,
                                                  GlobalIndexType& globIdxStart,
                                                  GlobalIndexType* neihgborsCount ) mutable
   {
      const PointType r_j = view_pointsToSearch[ j ];
      const PointType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );

      if( drs <= searchRadius )
      {
         neighborListStorage_view[ numberOfParticles * ( *neihgborsCount + anotherSetOffset + 1 ) + i ] = j;
         ( *neihgborsCount )++;
      }
   };

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const PointType r_i = view_points[ i ];
      const GlobalIndexType globIdxStart = i * neighborsCountLimit;
      GlobalIndexType neihgborsCount = 0;

      // in case we search in different set, offset to store neighbors
      int anotherSetOffset = 0;
      for( int k = 0; k < particlesToSearchLabel; k++ )
         anotherSetOffset += neighborListStorage_view[ anotherSetOffset * numberOfParticles + i ] + 1;

      BaseType::NeighborsLoopAnotherSet::exec( i, r_i, searchInSet, compareDistance, anotherSetOffset, globIdxStart, &neihgborsCount );
      neighborListStorage_view[ numberOfParticles *  anotherSetOffset  + i ] = neihgborsCount;
   };
   this->forAll( particleLoop );
}

template< typename ParticleConfig, typename DeviceType >
template< typename Layout,
          std::enable_if_t< std::is_same_v< Layout, NeighborListLayouts::ParticleMajorLinear >, bool > Enabled >
void
ParticlesLinkedListWithList< ParticleConfig, DeviceType >::buildParticleList()
{
   auto searchInSet = this->getCLLSearchToken();

   const RealType searchRadius = this->getSearchRadius();
   const auto view_points = this->getPoints().getConstView();
   auto neighborListStorage_view = this->neighborListStorage.getView();
   const GlobalIndexType neighborsCountLimit = this->neighborsCountLimit;
   const GlobalIndexType numberOfParticles = this->getNumberOfParticles();

   auto compareDistance = [=] __cuda_callable__ ( LocalIndexType i,
                                                  LocalIndexType j,
                                                  PointType& r_i,
                                                  GlobalIndexType& globIdxStart,
                                                  GlobalIndexType* neihgborsCount ) mutable
   {
      const PointType r_j = view_points[ j ];
      const PointType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius ){
         const GlobalIndexType globalIdx = globIdxStart + ( *neihgborsCount );
         neighborListStorage_view[ globalIdx + 1 ] = j;
         ( *neihgborsCount )++;
      }
   };

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const PointType r_i = view_points[ i ];
      const GlobalIndexType globIdxStart = i * neighborsCountLimit;
      GlobalIndexType neihgborsCount = 0;

      BaseType::NeighborsLoop::exec( i, r_i, searchInSet, compareDistance, globIdxStart, &neihgborsCount );
      neighborListStorage_view[ i ] = neihgborsCount;

   };
   this->forAll( particleLoop );
}

template< typename ParticleConfig, typename DeviceType >
template< typename ParticlesPointerType,
          typename Layout,
          std::enable_if_t< std::is_same_v< Layout, NeighborListLayouts::ParticleMajorLinear >, bool > Enabled >
void
ParticlesLinkedListWithList< ParticleConfig, DeviceType >::buildParticleList( ParticlesPointerType& particlesToSearch )
{
   auto searchInSet = particlesToSearch->getCLLSearchToken();

   const RealType searchRadius = this->getSearchRadius();
   const auto view_points = this->getPoints().getConstView();
   const auto view_pointsToSearch = particlesToSearch->getPoints().getConstView();
   auto neighborListStorage_view = this->neighborListStorage.getView();
   const GlobalIndexType neighborsCountLimit = this->neighborsCountLimit;
   const int particlesToSearchLabel = particlesToSearch->getParticleSetLabel();
   const GlobalIndexType numberOfParticles = this->getNumberOfParticles();

   auto compareDistance = [=] __cuda_callable__ ( LocalIndexType i,
                                                  LocalIndexType j,
                                                  PointType& r_i,
                                                  int& anotherSetOffset,
                                                  GlobalIndexType& globIdxStart,
                                                  GlobalIndexType* neihgborsCount ) mutable
   {
      const PointType r_j = view_pointsToSearch[ j ];
      const PointType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         neighborListStorage_view[ numberOfParticles * ( *neihgborsCount ) + i + 1 ] = j;
         ( *neihgborsCount )++;
      }
   };

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const PointType r_i = view_points[ i ];
      const GlobalIndexType globIdxStart = i * neighborsCountLimit;
      GlobalIndexType neihgborsCount = 0;

      // in case we search in different set, offset to store neighbors
      int anotherSetOffset = 0;
      for( int k = 0; k < particlesToSearchLabel; k++ )
         anotherSetOffset += neighborListStorage_view[ anotherSetOffset + i ] + 1;

      BaseType::NeighborsLoop::exec( i, r_i, searchInSet, compareDistance, anotherSetOffset, globIdxStart, &neihgborsCount );
      //neighborListStorage_view[ numberOfParticles * ( neihgborsCount ) + i + 1 ] = INT_MAX;
      neighborListStorage_view[ i ] = neihgborsCount;

      //neighborListStorage_view[ anotherSetOffset + globIdxStart ] = neihgborsCount;
      // setup default number of neighbors for the particles from second set
      //neighborListStorage_view[ anotherSetOffset + globIdxStart + neihgborsCount + 1 ] = 0;

   };
   //this->forAll( particleLoop );
   Algorithms::parallelFor< DeviceType >( 0, this->getNumberOfParticles(), particleLoop );
}

template< typename ParticleConfig, typename DeviceType >
void
ParticlesLinkedListWithList< ParticleConfig, DeviceType >::searchForNeighbors()
{
   BaseType::searchForNeighbors();
   this->buildParticleList();
}

template < typename ParticleConfig, typename DeviceType >
void
ParticlesLinkedListWithList< ParticleConfig, DeviceType >::writeProlog( TNL::Logger& logger ) const noexcept
{
   BaseType::writeProlog( logger );
   logger.writeParameter( "Particle set label:", this->particleSetLabel );
   logger.writeParameter( "Max. number of neighbors:", this->neighborsCountLimit );

   //Debug
   logger.writeParameter( "D: neighborListStorage.getSize(): ", this->neighborListStorage.getSize() );
}

} //namespace TNL
} //namespace Particles
