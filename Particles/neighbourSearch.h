#pragma once

#include <limits> //UINT_MAX
#include <utility> //std::forward

//#define UINT_MAX = std::numeric_limits<unsigned int>::max()

#include "Particles.h"

namespace TNL {
namespace ParticleSystem {

template< typename ParticleConfig, typename ParticleSystem >
class NeighborSearch
{
public:

  /* common */
  using DeviceType = typename ParticleSystem::Device; //mh
	using ParticlePointer = typename Pointers::SharedPointer< ParticleSystem, DeviceType >;

  using LocalIndexType = typename ParticleSystem::LocalIndexType;
  using GlobalIndexType = typename ParticleSystem::GlobalIndexType;
  using RealType = typename ParticleSystem::RealType;

	using CellIndexType = typename ParticleSystem::CellIndexType;
  using CellIndexArrayType = typename ParticleSystem::CellIndexArrayType;
	using CellIndexArrayView = typename Containers::ArrayView< typename ParticleSystem::CellIndexType, DeviceType >;

	/* Is this good idea? */
	using PointTypeArrayView = typename Containers::ArrayView< typename ParticleSystem::PointType, DeviceType >;
	/* Another temp stuff, just make it works, refucktor it away after. */
	using NeighborsCountView = typename Containers::ArrayView< typename ParticleSystem::LocalIndexType, DeviceType >;
	using NeighborsView = typename Containers::ArrayView< typename ParticleSystem::GlobalIndexType, DeviceType >;

  /* grid related */
  using GridType = typename ParticleSystem::GridType;
  using GridPointer = typename ParticleSystem::GridPointer;

  using myCell = typename ParticleSystem::myCell;
  //temp
  using GridCell = typename ParticleSystem::GridType::Cell;

  /**
   * Constructors.
   */
  NeighborSearch(ParticlePointer& particles, GlobalIndexType cellCount)
  : particles(particles), firstCellParticle(cellCount), lastCellParticle(cellCount)
  {
    //firstCellParticle = -1;
    //lastCellParticle = -1;

    firstCellParticle = INT_MAX;
    lastCellParticle = INT_MAX;

  }

  /**
   * Get list of first particle in cells.
   */
  const typename ParticleSystem::CellIndexArrayType& // -> using..
  getCellFirstParticleList() const;

  typename ParticleSystem::CellIndexArrayType& // -> using..
  getCellFirstParticleList();

  /**
   * Get list of last particle in cells.
   */
  const typename ParticleSystem::CellIndexArrayType& // -> using..
  getCellLastParticleList() const;

  typename ParticleSystem::CellIndexArrayType& // -> using..
  getCellLastParticleList();

  /**
   * Assign to each cell index of first contained particle.
   */
  __cuda_callable__
  void
  particlesToCells();

  /**
   * Run cycle over cells to search for neighbors.
   */
  __cuda_callable__
  void
  //getNeighborsFromTwoCells(Cell centralCell, Cell neighborCell);
  //getNeighborsFromTwoCells(LocalIndexType centralCell, LocalIndexType neighborCell, CellIndexArrayView view_firstCellParticle , CellIndexArrayView view_lastCellParticle, PointTypeArrayView view_particles );
  getNeighborsFromTwoCells(LocalIndexType centralCell, LocalIndexType neighborCell, CellIndexArrayView view_firstCellParticle , CellIndexArrayView view_lastCellParticle, PointTypeArrayView view_particles, NeighborsCountView view_neighborsCount, NeighborsView view_neighbors, RealType searchRadius);

  /**
   * Test particles in two neighbor cells.
   */
  __cuda_callable__
  void
  getAllNeighborCells();

	/**
	 * Particle negihbor loop.
	 */
	template< typename Function, typename... FunctionArgs >
  __cuda_callable__
  void
	neighborParticleLoop( Function f, FunctionArgs... args  ); //rename this

	template< typename Function, typename... FunctionArgs >
  __cuda_callable__
  void
	//onlyNeighborParticleLoop( GlobalIndexType i, Function f, FunctionArgs... args  ); //rename this
	onlyNeighborParticleLoop( const GlobalIndexType i, const GlobalIndexType numberOfParticles, const CellIndexArrayView view_firstCellParticle, const CellIndexArrayView view_particleCellIndex, Function f, FunctionArgs... args );

  /**
   * Test particles in two neighbor cells.
   */
  //__cuda_callable__
  void
  runCycleOverGrid();

  /**
   * Run the cycle to search for neighbors.
   */
  //__cuda_callable__
  void
  searchForNeighbors();

  /**
   * New version of neighbor search.
   */
  void
  searchForNeighborsFull(); //rename this

  void
  searchForNeighborsOnly(); //rename this

  /**
   * Run the cycle to search for neighbors.
   */
	void
  resetListWithIndices();


protected:

  //const ParticleSystem& particles;
  ParticlePointer particles;

  CellIndexArrayType firstCellParticle;
  CellIndexArrayType lastCellParticle;

};


} // ParticleSystem
} // TNL

#include "neighbourSearch_impl.h"

