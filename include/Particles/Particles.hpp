#include "Particles.h"

namespace TNL {
namespace ParticleSystem {

template< typename ParticleConfig, typename DeviceType >
constexpr int
Particles< ParticleConfig, DeviceType>::getParticlesDimension()
{
   return spaceDimension;
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
const typename Particles< ParticleConfig, DeviceType >::PointType
Particles< ParticleConfig, DeviceType >::getGridReferentialOrigin() const
{
   return gridReferentialOrigin;
}

template < typename ParticleConfig, typename DeviceType >
void
Particles< ParticleConfig, DeviceType >::setGridReferentialOrigin( const PointType& origin )
{
   gridReferentialOrigin = origin;
}

template < typename ParticleConfig, typename DeviceType >
const typename Particles< ParticleConfig, DeviceType >::PointType
Particles< ParticleConfig, DeviceType >::getGridOrigin() const
{
   return gridOrigin;
}

template < typename ParticleConfig, typename DeviceType >
void
Particles< ParticleConfig, DeviceType >::setGridOrigin( const PointType& origin )
{
   gridOrigin = origin;
}

template < typename ParticleConfig, typename DeviceType >
const typename Particles< ParticleConfig, DeviceType >::PointType
Particles< ParticleConfig, DeviceType >::getGridOriginWithOverlap() const
{
   const PointType originShift = this->overlapWidth * this->radius;
   const PointType shiftedGridOrigin = gridOrigin - originShift ;
   return shiftedGridOrigin;
}

template < typename ParticleConfig, typename DeviceType >
const typename Particles< ParticleConfig, DeviceType >::IndexVectorType
Particles< ParticleConfig, DeviceType >::getGridDimensions() const
{
   return gridDimension;
}


template < typename ParticleConfig, typename DeviceType >
void
Particles< ParticleConfig, DeviceType >::setGridDimensions( const IndexVectorType& dimensions )
{
   gridDimension = dimensions;
}

template < typename ParticleConfig, typename DeviceType >
const typename Particles< ParticleConfig, DeviceType >::IndexVectorType
Particles< ParticleConfig, DeviceType >::getGridDimensionsWithOverlap() const
{
   const IndexVectorType resizeGridDimensions = 2 * this->overlapWidth;
   const IndexVectorType resizedGridDimensions = gridDimension + resizeGridDimensions;
   return resizedGridDimensions;
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

//FIXME: I don't like this getPoint part, it is used only in VTK writer. Get this out.
//template < typename ParticleConfig, typename DeviceType >
//__cuda_callable__
//const typename Particles< ParticleConfig, DeviceType >::PointType&
//Particles< ParticleConfig, DeviceType >::getPoint(GlobalIndexType particleIndex) const
//{
//   TNL_ASSERT_GE( particleIndex, 0, "invalid particle index" );
//   TNL_ASSERT_LT( particleIndex, numberOfParticles, "invalid particle index" );
//   return this->points[ particleIndex ];
//}
//
//template < typename ParticleConfig, typename DeviceType >
//__cuda_callable__
//typename Particles< ParticleConfig, DeviceType >::PointType&
//Particles< ParticleConfig, DeviceType >::getPoint(GlobalIndexType particleIndex)
//{
//   TNL_ASSERT_GE( particleIndex, 0, "invalid particle index" );
//   TNL_ASSERT_LT( particleIndex, numberOfParticles, "invalid particle index" );
//   return this->points[ particleIndex ];
//}
//
//template < typename ParticleConfig, typename DeviceType >
//__cuda_callable__
//void
//Particles<ParticleConfig, DeviceType>::setPoint(GlobalIndexType particleIndex, PointType point)
//{
//   this->points[ particleIndex ] = point;
//}

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

template < typename ParticleConfig, typename Device >
const typename Particles< ParticleConfig, Device >::GlobalIndexType
Particles< ParticleConfig, Device >::getNumberOfParticlesToRemove() const
{
   return numberOfParticlesToRemove;
}

template < typename ParticleConfig, typename Device >
void
Particles< ParticleConfig, Device >::setNumberOfParticlesToRemove( const GlobalIndexType& removeCount )
{
   this->numberOfParticlesToRemove = removeCount;
}

template< typename ParticleConfig, typename Device >
template< typename Device2, typename Func >
void
Particles< ParticleConfig, Device >::forAll( Func f ) const
{
   const PointType domainOrigin = this->gridOrigin;
   const PointType domainSize = this->gridDimension * this->radius;
   const auto view_points = this->points.getConstView();
   auto wrapper = [=] __cuda_callable__( GlobalIndexType i ) mutable
   {
      if( this->isInsideDomain( view_points[ i ], domainOrigin, domainSize ) )
         f( i );
   };
   Algorithms::parallelFor< Device2 >( 0, numberOfParticles, wrapper );
}

template < typename ParticleConfig, typename Device >
__cuda_callable__
bool
Particles< ParticleConfig, Device >::isInsideDomain( const PointType& point,
                                                               const PointType& domainOrigin,
                                                               const PointType& domainSize ) const
{
   //FIXME: These two lines produces different results
   //if( ( point > domainOrigin ) && ( point < ( domainOrigin + domainSize ) ) )
   if( ( point[ 0 ] >= domainOrigin[ 0 ] ) && ( point[ 0 ] < ( domainOrigin[ 0 ] + domainSize[ 0 ] ) ) ) // >=, <= vs >, <
      return true;
   return false;
}

template < typename ParticleConfig, typename Device >
void
Particles< ParticleConfig, Device >::reorderParticles()
{
   const GlobalIndexType numberOfParticle = this->getNumberOfParticles();
   using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< Device >;
   ThrustDeviceType thrustDevice;
   auto view_map = this->sortPermutations->getView();
   auto view_points = this->getPoints().getView();
   auto view_points_swap = this->points_swap.getView();
   thrust::gather( thrustDevice,
                   view_map.getArrayData(),
                   view_map.getArrayData() + numberOfParticle,
                   view_points.getArrayData(),
                   view_points_swap.getArrayData() );
   this->getPoints().swap( this->points_swap );
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
