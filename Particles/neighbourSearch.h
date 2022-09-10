#pragma once

#include "Particles.h"

namespace TNL {
namespace ParticleSystem {

template< typename ParticleConfig, typename ParticleSystem >
class NeighborSearch
{
public:

  /* common */
  using DeviceType = typename ParticleSystem::Device; //mh

  using LocalIndexType = typename ParticleSystem::LocalIndexType;
  using GlobalIndexType = typename ParticleSystem::GlobalIndexType;

  using CellIndexArrayType = typename ParticleSystem::CellIndexArrayType;

  /* grid related */
  using GridType = typename ParticleSystem::GridType;
  using GridPointer = typename ParticleSystem::GridPointer;

  using myCell = typename ParticleSystem::myCell;
  //temp
  using GridCell = typename ParticleSystem::GridType::Cell;

  /**
   * Constructors.
   */
  NeighborSearch(ParticleSystem& particles, GlobalIndexType cellCount)
  : particles(particles), firstCellParticle(cellCount), lastCellParticle(cellCount)
  {
    firstCellParticle = -1;
    lastCellParticle = -1;
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
  getNeighborsFromTwoCells(LocalIndexType centralCell, LocalIndexType neighborCell);

  /**
   * Test particles in two neighbor cells.
   */
  __cuda_callable__
  void
  getAllNeighborCells();

  /**
   * Test particles in two neighbor cells.
   */
  __cuda_callable__
  void
  runCycleOverGrid();

  /**
   * Run the cycle to search for neighbors.
   */
  __cuda_callable__
  void
  searchForNeighbors();


protected:

  //const ParticleSystem& particles;
  ParticleSystem& particles;

  CellIndexArrayType firstCellParticle;
  CellIndexArrayType lastCellParticle;


};


} // ParticleSystem
} // TNL

#include "neighbourSearch_impl.h"

