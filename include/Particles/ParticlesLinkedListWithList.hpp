#include "ParticlesLinkedListWithList.h"

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

/*
template< typename ParticleConfig, typename DeviceType >
void
ParticlesLinkedListWithList< ParticleConfig, DeviceType >::buildParticleList()
{
   auto searchInSet = this->getSearchToken();

   const RealType searchRadius = this->getSearchRadius();
   const auto view_points = this->getPoints().getConstView();
   auto neighborListStorage_view = this->neighborListStorage.getView();
   const GlobalIndexType neighborsCountLimit = this->neighborsCountLimit;

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
      neighborListStorage_view[ globIdxStart ] = neihgborsCount;
      // setup default number of neighbors for the particles from second set
      neighborListStorage_view[ globIdxStart + neihgborsCount + 1 ] = 0;
   };
   this->forAll( particleLoop );
}
*/

template< typename ParticleConfig, typename DeviceType >
template< typename ParticlesPointerType >
void
ParticlesLinkedListWithList< ParticleConfig, DeviceType >::buildParticleList( ParticlesPointerType& particlesToSearch )
{
   auto searchInSet = particlesToSearch->getSearchToken();

   const RealType searchRadius = this->getSearchRadius();
   const auto view_points = this->getPoints().getConstView();
   const auto view_pointsToSearch = particlesToSearch->getPoints().getConstView();
   auto neighborListStorage_view = this->neighborListStorage.getView();
   const GlobalIndexType neighborsCountLimit = this->neighborsCountLimit;

   auto compareDistance = [=] __cuda_callable__ ( LocalIndexType i,
                                                  LocalIndexType j,
                                                  PointType& r_i,
                                                  GlobalIndexType& globIdxStart,
                                                  GlobalIndexType* neihgborsCount ) mutable
   {
      const PointType r_j = view_pointsToSearch[ j ];
      const PointType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
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
      neighborListStorage_view[ globIdxStart ] = neihgborsCount;
      // setup default number of neighbors for the particles from second set
      neighborListStorage_view[ globIdxStart + neihgborsCount + 1 ] = 0;
   };
   this->forAll( particleLoop );
}

template< typename ParticleConfig, typename DeviceType >
template< typename ParticlesPointerType >
void
ParticlesLinkedListWithList< ParticleConfig, DeviceType >::searchForNeighbors( ParticlesPointerType& particlesToSearch )
{
   BaseType::searchForNeighbors();
   this->buildParticleList( particlesToSearch );
}

} //namespace TNL
} //namespace Particles
