#include "Particles.h"

namespace TNL {
namespace ParticleSystem {

template< typename ParticleConfig, typename DeviceType >
constexpr int
Particles< ParticleConfig, DeviceType>::getParticleDimension()
{
   return spaceDimension;
}

template< typename ParticleConfig, typename DeviceType >
__cuda_callable__
typename Particles< ParticleConfig, DeviceType >::GlobalIndexType
Particles< ParticleConfig, DeviceType>::getNumberOfParticles()
{
   return numberOfParticles;
}

template< typename ParticleConfig, typename DeviceType >
__cuda_callable__
const typename Particles< ParticleConfig, DeviceType >::GlobalIndexType
Particles< ParticleConfig, DeviceType>::getNumberOfParticles() const
{
   return numberOfParticles;
}

template< typename ParticleConfig, typename DeviceType >
__cuda_callable__
const typename Particles< ParticleConfig, DeviceType >::GlobalIndexType
Particles< ParticleConfig, DeviceType>::getNumberOfAllocatedParticles() const
{
   return numberOfAllocatedParticles;
}

template< typename ParticleConfig, typename DeviceType >
void
Particles< ParticleConfig, DeviceType >::setNumberOfParticles( const GlobalIndexType& newNumberOfParticles )
{
   this->numberOfParticles = newNumberOfParticles;
}

template< typename ParticleConfig, typename DeviceType >
void
Particles< ParticleConfig, DeviceType >::setSize( const GlobalIndexType& size )
{
   numberOfAllocatedParticles = size;
   this->points.setSize( size );
   this->points_swap.setSize( size );
   this->sortPermutations->setSize( size ); //TODO: sort permulation should not be pointer.
}

template< typename ParticleConfig, typename DeviceType >
__cuda_callable__
const typename Particles< ParticleConfig, DeviceType >::RealType
Particles< ParticleConfig, DeviceType>::getSearchRadius() const
{
   return radius;
}

template< typename ParticleConfig, typename DeviceType >
void
Particles< ParticleConfig, DeviceType>::setSearchRadius( const RealType& searchRadius )
{
   this->radius = searchRadius;
}

template < typename ParticleConfig, typename DeviceType >
const typename Particles< ParticleConfig, DeviceType >::ParticleTraitsType::PointArrayType&
Particles< ParticleConfig, DeviceType >::getPoints() const
{
   return points;
}

template < typename ParticleConfig, typename DeviceType >
typename Particles< ParticleConfig, DeviceType >::ParticleTraitsType::PointArrayType&
Particles< ParticleConfig, DeviceType >::getPoints()
{
   return points;
}

template < typename ParticleConfig, typename DeviceType >
const typename Particles< ParticleConfig, DeviceType >::ParticleTraitsType::PointArrayType&
Particles< ParticleConfig, DeviceType >::getPointsSwap() const
{
   return points_swap;
}

template < typename ParticleConfig, typename DeviceType >
typename Particles< ParticleConfig, DeviceType >::ParticleTraitsType::PointArrayType&
Particles< ParticleConfig, DeviceType >::getPointsSwap()
{
   return points_swap;
}

template < typename ParticleConfig, typename DeviceType >
__cuda_callable__
const typename Particles< ParticleConfig, DeviceType >::PointType&
Particles< ParticleConfig, DeviceType >::getPoint(GlobalIndexType particleIndex) const
{
   TNL_ASSERT_GE( particleIndex, 0, "invalid particle index" );
   TNL_ASSERT_LT( particleIndex, numberOfParticles, "invalid particle index" );
   return this->points[ particleIndex ];
}

template < typename ParticleConfig, typename DeviceType >
__cuda_callable__
typename Particles< ParticleConfig, DeviceType >::PointType&
Particles< ParticleConfig, DeviceType >::getPoint(GlobalIndexType particleIndex)
{
   TNL_ASSERT_GE( particleIndex, 0, "invalid particle index" );
   TNL_ASSERT_LT( particleIndex, numberOfParticles, "invalid particle index" );
   return this->points[ particleIndex ];
}

template < typename ParticleConfig, typename DeviceType >
__cuda_callable__
void
Particles<ParticleConfig, DeviceType>::setPoint(GlobalIndexType particleIndex, PointType point)
{
   this->points[ particleIndex ] = point;
}

template < typename ParticleConfig, typename DeviceType >
const typename Particles< ParticleConfig, DeviceType >::IndexArrayTypePointer&
Particles< ParticleConfig, DeviceType >::getSortPermutations() const
{
   return this->sortPermutations;
}

template < typename ParticleConfig, typename DeviceType >
typename Particles< ParticleConfig, DeviceType >::IndexArrayTypePointer&
Particles< ParticleConfig, DeviceType >::getSortPermutations()
{
   return this->sortPermutations;
}

template < typename ParticleConfig, typename DeviceType >
void
Particles< ParticleConfig, DeviceType >::writeProlog( TNL::Logger& logger ) const noexcept
{
   logger.writeParameter( "Number of particles:", this->numberOfParticles );
   logger.writeParameter( "Number of allocated particles:", this->numberOfAllocatedParticles );
   logger.writeParameter( "Search radius:", this->radius );
}

} //namespace TNL
} //namespace Particles
