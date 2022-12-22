#pragma once

#include <limits> //UINT_MAX
#include <utility> //std::forward

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

   /* bucketing */
   using PairIndexType = Containers::StaticVector< 2, GlobalIndexType >;
   using PairIndexArrayType = Containers::Array< PairIndexType, DeviceType, GlobalIndexType >;
   using PairIndexArrayView = typename Containers::ArrayView< PairIndexType, DeviceType >;

   /**
    * Constructors.
    */
   NeighborSearch(ParticlePointer& particles, GlobalIndexType cellCount)
   : particles(particles), firstLastCellParticle(cellCount)
   {
      firstLastCellParticle = INT_MAX;
   }

   /**
    * Get list of first and last particle in cells.
    */
   const PairIndexArrayType& // -> using..
   getCellFirstLastParticleList() const;

   PairIndexArrayType& // -> using..
   getCellFirstLastParticleList();

   /**
    * Assign to each cell index of first contained particle.
    */
   __cuda_callable__
   void
   particlesToCells();

   /**
    * Runs cycle over all particles, for each particles run cycle over ints negihbors
    * and exec given function f for each pair.
    */
   template< typename Function, typename... FunctionArgs >
   __cuda_callable__
   void
   loopOverParticlesAndNeighbors( Function f, FunctionArgs... args  ); //rename this

   /**
    * Loop that for given particle i runs cycle over all its neighbors
    * and exec given function f for each pair.
    * TODO: Improve the interface and argument handling.
    */
   template< typename Function, typename... FunctionArgs >
   __cuda_callable__
   void
   loopOverNeighbors( const GlobalIndexType i, const GlobalIndexType& numberOfParticles, const PairIndexArrayView& view_firstLastCellParticle, const CellIndexArrayView& view_particleCellIndex, Function f, FunctionArgs... args );

   /**
    * For all particles run loop over neighbors and assemble neighbor list.
     * TEMP: This function is actually not used, but show the structure of negihbor loop.
    */
   void
   searchForNeighborsExampleLoop();

   /**
    * Runs all necessary functions and fills neighbor list.
    */
   void
   searchForNeighbors();

   /**
    * Search for neighbors and assemble neighbor list.
     * TEMP: This function has some limitations, but show how the neighbor loop should work.
    */
   void
   searchForNeighborsWithForAll();

   /**
    * Search for neighbors and assemble neighbor list.
     */
   void
   searchForNeighborsWithForEach();

   /**
    * Reset the list with first and last particle in cell.
    */
   void
   resetListWithIndices(); //protected?

protected:

   //int numberOfCells
   ParticlePointer particles;

   PairIndexArrayType firstLastCellParticle;

};

} // ParticleSystem
} // TNL

#include "neighbourSearch_impl.h"

