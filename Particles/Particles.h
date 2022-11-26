#pragma once

#include <TNL/Algorithms/ParallelFor.h>
#include <TNL/Pointers/SharedPointer.h>
#include <TNL/Algorithms/sort.h>

#include "ParticlesTraits.h"
#include "GenerateCellIndex.h"


namespace TNL {
namespace ParticleSystem {

template < typename ParticleConfig, typename DeviceType >
class Particles{
public:

  using Device = DeviceType;
  using Config = ParticleConfig;
  using ParticleTraitsType = ParticlesTraits< Config, DeviceType >;

  static constexpr int spaceDimension = 2;

  /* common */
  using GlobalIndexType = typename ParticleTraitsType::GlobalIndexType;
  using LocalIndexType = typename ParticleTraitsType::LocalIndexType;
  using RealType = typename ParticleTraitsType::RealType;
  using CellIndexType = typename ParticleTraitsType::CellIndexType;
  using CellIndexArrayType = typename ParticleTraitsType::CellIndexArrayType; //nn
  using NeighborsCountArrayType = typename ParticleTraitsType::NeighborsCountArrayType;
  using NeighborsArrayType = typename ParticleTraitsType::NeighborsArrayType;
  using NeighborListType = typename ParticleTraitsType::NeighborListType;

  using CellIndexer = typename Config::CellIndexerType;

  /* particle related */
  using PointType = typename ParticleTraitsType::PointType;
  using PointArrayType = typename ParticleTraitsType::PointArrayType; //nn

  /* grid related */
  using GridType = typename ParticleTraitsType::GridType;
  using GridPointer = typename ParticleTraitsType::GridPointer;

  //temp
  //old using myCell = Meshes::GridEntity< GridType, 2, Meshes::GridEntityCrossStencilStorage<  > >;
  using myCell = Meshes::GridEntity< GridType, 2 >;

  /**
   * Constructors.
   */
  Particles() = default;

  Particles(GlobalIndexType size)
  : numberOfParticles(size), points(size)
  {

  }

  Particles(GlobalIndexType size, RealType radius)
  : numberOfParticles(size), points(size), radius(radius), particleCellInidices(size), gridCellIndices(Config::gridXsize*Config::gridYsize), neighborsCount(size, 0), neighbors(size*Config::maxOfNeigborsPerParticle, 0)
  {
    //grid->setSpaceSteps( { Config::searchRadius, Config::searchRadius } ); //removed
    grid->setDimensions( Config::gridXsize, Config::gridYsize );
    //grid->setDimensions( Config::gridXsize, Config::gridYsize ); //ignore boundary cells

    //gridCellIndices(grid->template getEntitiesCount<2>());
    //particleCellInidices = 0; //reset

    //grid->setOrigin( { Config::gridXbegin, Config::gridYbegin } ); //removed
    //grid->setOrigin( { Config::gridXbegin , Config::gridYbegin } ); //ignore boundary cells

    neighborsList.setSegmentsSizes(size, Config::maxOfNeigborsPerParticle);

  }

  /* PARTICLE RELATED TOOLS */

  /**
   * Get dimension of particle system.
   */
  static constexpr int
  getParticleDimension();

  /**
   * Get search radius.
   */
  __cuda_callable__
  RealType
  getSearchRadius() const;

  /**
   * Get number of particles in particle system.
   */
  __cuda_callable__
  GlobalIndexType
  getNumberOfParticles() const;

  /**
   * Get particle (i.e. point) positions.
   */
  const typename ParticleTraitsType::PointArrayType& // -> using..
  getPoints() const;

  typename ParticleTraitsType::PointArrayType& // -> using..
  getPoints();

  /**
   * Get position of given particle.
   */
  __cuda_callable__
  const PointType&
  getPoint( GlobalIndexType particleIndex ) const;

  __cuda_callable__
  PointType&
  getPoint( GlobalIndexType particleIndex );

  /**
   * Set position of given particle.
   */
  __cuda_callable__
  void
  setPoint( GlobalIndexType particleIndex, PointType point);

  /**
   * Get particle cell indices.
   */
  const typename ParticleTraitsType::CellIndexArrayType& // -> using..
  getParticleCellIndices() const;

  typename ParticleTraitsType::CellIndexArrayType& // -> using..
  getParticleCellIndices();

  /**
   * Get cell index of given partile.
   */
  __cuda_callable__
  const CellIndexType&
  getParticleCellIndex( GlobalIndexType particleIndex ) const;

  __cuda_callable__
  CellIndexType&
  getParticleCellIndex( GlobalIndexType particleIndex );

  /**
   * Get cell index of given partile.
   */
  void
  computeParticleCellIndices();

  /**
   * Sort particles by its cell index.
   */
  void sortParticles();

  /* PARTICLE RELATED TEMP TOOLS */

  void generateRandomParticles(); //used only for tests

  /* GRID RELATED TOOLS */

  void // -> move this into constructor..
  SetupGrid();

  const typename ParticleTraitsType::CellIndexArrayType& // -> using..
  getGridCellIndices() const;

  typename ParticleTraitsType::CellIndexArrayType& // -> using..
  getGridCellIndices();

  __cuda_callable__
  const CellIndexType&
  getGridCellIndex( GlobalIndexType cellIndex ) const;

  __cuda_callable__
  CellIndexType&
  getGridCellIndex( GlobalIndexType cellIndex );

  void computeGridCellIndices();

  /* general */
  void GetParticlesInformations();

  /* NEIGHBOR LIST RELATED TOOL */

  /**
   * Return jth neighbor of particle i.
   */
  __cuda_callable__
  const GlobalIndexType&
  getNeighbor( GlobalIndexType i, GlobalIndexType j ) const;

  __cuda_callable__
  GlobalIndexType&
  getNeighbor( GlobalIndexType i, GlobalIndexType j );

  /**
   * Return list with numbers of particles.
   */
  const NeighborsCountArrayType& // -> using..
  getNeighborsCountList() const;

  NeighborsCountArrayType& // -> using..
  getNeighborsCountList();

  /**
   * Return number of neighbors of particle i.
   */
  __cuda_callable__
  const LocalIndexType&
  getNeighborsCount( GlobalIndexType particleIndex ) const;

  __cuda_callable__
  LocalIndexType&
  getNeighborsCount( GlobalIndexType particleIndex );

  /**
   * Set j as neighbor for particle i.
   */
  __cuda_callable__
  void
  setNeighbor( GlobalIndexType i, GlobalIndexType j );

  /**
   * Remove all neighbors and clear the neighbor list.
   */
  void
  resetNeighborList();

  /* gird related */
  GridPointer grid; //temporarily moved outside protected

  /* NEIGHBOR LIST RELATED TEMP TOOL */
  //temp - off
  __cuda_callable__
  myCell&
  getCell( GlobalIndexType i );

	//temp
  __cuda_callable__
	GlobalIndexType
	getCountOfGridCells();

  /**
   * Print/save neighbor whole neighbor list.
   */
  void
  saveNeighborList(std::string neigborListFile);


protected:

  /* particle related*/
  GlobalIndexType numberOfParticles;
  GlobalIndexType gridSize;

  RealType radius;

  NeighborsCountArrayType neighborsCount;
  NeighborsArrayType neighbors;

  NeighborListType neighborsList; //not used

  //PointArrayType* points = nullptr;
  PointArrayType points;
  CellIndexArrayType particleCellInidices;


  CellIndexArrayType gridCellIndices;

};



} //namespace Particles

} //namespace TNL

#include "Particles_impl.h"

