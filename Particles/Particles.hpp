#include "Particles.h"

namespace TNL {
namespace ParticleSystem {

/* PARTICLE RELATED TOOLS */

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
Particles< ParticleConfig, DeviceType >::setNumberOfParticles( GlobalIndexType newNumberOfParticles )
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
__cuda_callable__
const typename Particles< ParticleConfig, DeviceType >::IndexArrayTypePointer&
Particles< ParticleConfig, DeviceType >::getSortPermutations() const
{
   return this->sortPermutations;
}

template < typename ParticleConfig, typename DeviceType >
__cuda_callable__
typename Particles< ParticleConfig, DeviceType >::IndexArrayTypePointer&
Particles< ParticleConfig, DeviceType >::getSortPermutations()
{
   return this->sortPermutations;
}

//template < typename ParticleConfig, typename DeviceType >
//void
//Particles< ParticleConfig, DeviceType >::computeParticleCellIndices()
//{
//   GlobalIndexType _numberOfParticles = this->numberOfParticles;
//
//   auto view = this->particleCellInidices.getView();
//   auto view_points = this->points.getView();
//
//   CellIndexer::ComputeParticleCellIndex( view, view_points, _numberOfParticles, gridDimension, gridOrigin, radius );
//}

//template < typename ParticleConfig, typename DeviceType >
//void
//Particles< ParticleConfig, DeviceType >::sortParticles()
//{
//   auto view_particleCellIndices = this->particleCellInidices.getView();
//   auto view_points = this->points.getView();
//   Algorithms::sort< DeviceType, GlobalIndexType >(
//       0, this->numberOfParticles,
//       [=] __cuda_callable__ ( int i, int j ) -> bool {
//         return view_particleCellIndices[ i ] < view_particleCellIndices[ j ]; },
//       [=] __cuda_callable__ ( int i, int j ) mutable {
//         swap( view_particleCellIndices[ i ], view_particleCellIndices[ j ] );
//         swap( view_points[ i ], view_points[ j ] ); } );
//}

//template < typename ParticleConfig, typename DeviceType >
//void
//Particles< ParticleConfig, DeviceType >::sortParticles()
//{
//
//   GlobalIndexType numberOfParticle = this->getNumberOfParticles();
//   auto view_particleCellIndices = this->getParticleCellIndices().getView();
//   auto view_map = this->sortPermutations->getView();
//
//   sortPermutations->forAllElements( [] __cuda_callable__ ( int i, int& value ) { value = i; } );
//   thrust::sort_by_key( thrust::device, view_particleCellIndices.getArrayData(),
//         view_particleCellIndices.getArrayData() + numberOfParticle, view_map.getArrayData() ); //TODO: replace thrust::device
//
//   auto view_points = this->getPoints().getView();
//   auto view_points_swap = this->points_swap.getView();
//   thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticle,
//         view_points.getArrayData(), view_points_swap.getArrayData() );
//   this->getPoints().swap( this->points_swap );
//}

/* PARTICLE RELATED TEMP TOOLS */

template < typename ParticleConfig, typename DeviceType >
void
Particles< ParticleConfig, DeviceType >::generateRandomParticles()
{
   std::cout << " Doesn't work at this moment. " << std::endl;
}

/* GRID RELATED TOOLS */
template < typename ParticleConfig, typename DeviceType >
__cuda_callable__ //TODO: Comment.
const typename Particles< ParticleConfig, DeviceType >::IndexVectorType
Particles< ParticleConfig, DeviceType >::getGridSize() const
{
   return gridDimension;
}

template < typename ParticleConfig, typename DeviceType >
void
Particles< ParticleConfig, DeviceType >::setGridSize( IndexVectorType gridSize )
{
   gridDimension = gridSize;
}

template < typename ParticleConfig, typename DeviceType >
const typename Particles< ParticleConfig, DeviceType >::PointType
Particles< ParticleConfig, DeviceType >::getGridOrigin() const
{
   return gridOrigin;
}

template < typename ParticleConfig, typename DeviceType >
void
Particles< ParticleConfig, DeviceType >::setGridOrigin( PointType gridBegin )
{
   gridOrigin = gridBegin; //FIXME: Names.
}

} //namespace TNL
} //namespace Particles
