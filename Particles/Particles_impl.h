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
Particles< ParticleConfig, DeviceType>::getNumberOfParticles() const
{
  return numberOfParticles;
}

template< typename ParticleConfig, typename DeviceType >
__cuda_callable__
typename Particles< ParticleConfig, DeviceType >::RealType
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
const typename Particles< ParticleConfig, DeviceType >::ParticleTraitsType::CellIndexArrayType&
Particles< ParticleConfig, DeviceType >::getParticleCellIndices() const
{
  return particleCellInidices;
}

template < typename ParticleConfig, typename DeviceType >
typename Particles< ParticleConfig, DeviceType >::ParticleTraitsType::CellIndexArrayType&
Particles< ParticleConfig, DeviceType >::getParticleCellIndices()
{
  return particleCellInidices;
}

template < typename ParticleConfig, typename DeviceType >
__cuda_callable__
const typename Particles< ParticleConfig, DeviceType >::CellIndexType&
Particles< ParticleConfig, DeviceType >::getParticleCellIndex(GlobalIndexType particleIndex) const
{
  TNL_ASSERT_GE( particleIndex, 0, "invalid particle index" );
  TNL_ASSERT_LT( particleIndex, numberOfParticles, "invalid particle index" );
  return this->particleCellInidices[ particleIndex ];
}

template < typename ParticleConfig, typename DeviceType >
__cuda_callable__
typename Particles< ParticleConfig, DeviceType >::CellIndexType&
Particles< ParticleConfig, DeviceType >::getParticleCellIndex(GlobalIndexType particleIndex)
{
  TNL_ASSERT_GE( particleIndex, 0, "invalid particle index" );
  TNL_ASSERT_LT( particleIndex, numberOfParticles, "invalid particle index" );
  return this->particleCellInidices[ particleIndex ];
}

template < typename ParticleConfig, typename DeviceType >
void
Particles< ParticleConfig, DeviceType >::computeParticleCellIndices()
{
  GlobalIndexType _numberOfParticles = this->numberOfParticles;
  GlobalIndexType _numberOfCells = ParticleConfig::gridYsize;

  auto view = this->particleCellInidices.getView();
  auto view_points = this->points.getView();

  CellIndexer::ComputeParticleCellIndex( view, view_points, _numberOfParticles );

}

template < typename ParticleConfig, typename DeviceType >
void
Particles< ParticleConfig, DeviceType >::sortParticles()
{
  auto view_particleCellIndices = this->particleCellInidices.getView();
  auto view_points = this->points.getView();
  Algorithms::sort< DeviceType, GlobalIndexType >(
      0, this->numberOfParticles,
      [=] __cuda_callable__ ( int i, int j ) -> bool {
        return view_particleCellIndices[ i ] < view_particleCellIndices[ j ]; },
      [=] __cuda_callable__ ( int i, int j ) mutable {
        swap( view_particleCellIndices[ i ], view_particleCellIndices[ j ] );
        swap( view_points[ i ], view_points[ j ] ); } );
}

/* PARTICLE RELATED TEMP TOOLS */

template < typename ParticleConfig, typename DeviceType >
void
Particles< ParticleConfig, DeviceType >::generateRandomParticles()
{
  typename Particles< ParticleConfig, DeviceType>::PointArrayType aux_points(this->numberOfParticles);
  aux_points.forAllElements( [=] __cuda_callable__ ( LocalIndexType i, PointType& value )
      { value[0] = static_cast< RealType >( std::rand() / static_cast< RealType >( RAND_MAX / ( 8 - 0 ) ));
        value[1] = static_cast< RealType >( std::rand() / static_cast< RealType >( RAND_MAX / ( 8 - 0 ) )); });
  this->points = aux_points;

}

/* GRID RELATED TOOLS */


template < typename ParticleConfig, typename DeviceType >
void
Particles< ParticleConfig, DeviceType >::SetupGrid()
{

}

template < typename ParticleConfig, typename DeviceType >
const typename Particles< ParticleConfig, DeviceType >::ParticleTraitsType::CellIndexArrayType&
Particles< ParticleConfig, DeviceType >::getGridCellIndices() const
{
  return gridCellIndices;
}

template < typename ParticleConfig, typename DeviceType >
typename Particles< ParticleConfig, DeviceType >::ParticleTraitsType::CellIndexArrayType&
Particles< ParticleConfig, DeviceType >::getGridCellIndices()
{
  return gridCellIndices;
}

template < typename ParticleConfig, typename DeviceType >
__cuda_callable__
const typename Particles< ParticleConfig, DeviceType >::CellIndexType&
Particles< ParticleConfig, DeviceType >::getGridCellIndex(GlobalIndexType cellIndex) const
{
  TNL_ASSERT_GE( cellIndex, 0, "invalid cell index" );
  TNL_ASSERT_LT( cellIndex, numberOfParticles, "invalid cell index" ); //CIDX
  return this->gridCellIndices[ cellIndex ];
}

template < typename ParticleConfig, typename DeviceType >
__cuda_callable__
typename Particles< ParticleConfig, DeviceType >::CellIndexType&
Particles< ParticleConfig, DeviceType >::getGridCellIndex(GlobalIndexType cellIndex)
{
  TNL_ASSERT_GE( cellIndex, 0, "invalid cell index" );
  TNL_ASSERT_LT( cellIndex, numberOfParticles, "invalid cell index" ); //CIDX
  return this->gridCellIndices[ cellIndex ];
}

template < typename ParticleConfig, typename DeviceType >
void
Particles< ParticleConfig, DeviceType >::computeGridCellIndices()
{

  auto view = this->gridCellIndices.getView();
  auto view_points = this->points.getView();

  //CellIndexView viewVIEW = this->gridCellIndices.getView(); //or this way?
  //PointsView view_pointsVIEW = this->points.getView(); //or this way?

  CellIndexer::ComputeCellIndex(view, view_points);
}

/* general */
template < typename ParticleConfig, typename DeviceType >
void
Particles< ParticleConfig, DeviceType>::GetParticlesInformations()
{
  std::cout << "Number of particles: " << numberOfParticles << std::endl;
  std::cout << "Search radius: " << radius << std::endl;
  std::cout << "Grid details:\n" << *grid << std::endl;

  std::cout << "Neighbor list: " << neighborsList << std::endl;
}

/* NEIGHBOR LIST RELATED TOOLS */

template < typename ParticleConfig, typename DeviceType >
__cuda_callable__
const typename Particles< ParticleConfig, DeviceType >::GlobalIndexType&
Particles< ParticleConfig, DeviceType >::getNeighbor(GlobalIndexType i, GlobalIndexType j) const
{
  return neighbors[(ParticleConfig::maxOfNeigborsPerParticle + 1)*i + j];
}

template < typename ParticleConfig, typename DeviceType >
__cuda_callable__
typename Particles< ParticleConfig, DeviceType >::GlobalIndexType&
Particles< ParticleConfig, DeviceType >::getNeighbor(GlobalIndexType i, GlobalIndexType j)
{
  return neighbors[(ParticleConfig::maxOfNeigborsPerParticle + 1)*i + j];
}

template < typename ParticleConfig, typename DeviceType >
__cuda_callable__
const typename Particles< ParticleConfig, DeviceType >::LocalIndexType&
Particles< ParticleConfig, DeviceType >::getNeighborsCount(GlobalIndexType i) const
{
  return neighbors[(ParticleConfig::maxOfNeigborsPerParticle)*i];
}


template < typename ParticleConfig, typename DeviceType >
__cuda_callable__
typename Particles< ParticleConfig, DeviceType >::LocalIndexType&
Particles< ParticleConfig, DeviceType >::getNeighborsCount(GlobalIndexType i)
{
  return neighbors[(ParticleConfig::maxOfNeigborsPerParticle)*i];
}

template < typename ParticleConfig, typename DeviceType >
__cuda_callable__
void
Particles< ParticleConfig, DeviceType >::setNeighbor(GlobalIndexType i, GlobalIndexType j)
{
  TNL_ASSERT_LT( neighbors[(ParticleConfig::maxOfNeigborsPerParticle)*i]++, ParticleConfig::maxOfNeigborsPerParticle, "number of neighbors is reached" );

  neighbors[(ParticleConfig::maxOfNeigborsPerParticle)*i]++;
  neighbors[(ParticleConfig::maxOfNeigborsPerParticle)*i + neighbors[(ParticleConfig::maxOfNeigborsPerParticle)*i]] = j;
}


} //namespace TNL
} //namespace Particles
