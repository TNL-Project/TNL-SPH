#pragma once

#include <limits> //UINT_MAX
#include <utility> //std::forward

#include "Particles.h"

namespace TNL {
namespace ParticleSystem {

template< typename NeighborSearch >
class NeighborsLoopParams
{
public:
   using DeviceType = typename NeighborSearch::DeviceType;
   using GlobalIndexType = typename NeighborSearch::GlobalIndexType;
   using PairIndexType = Containers::StaticVector< 2, GlobalIndexType >;
   using CellIndexArrayView = typename Containers::ArrayView< typename NeighborSearch::CellIndexType, DeviceType >;
   using PairIndexArrayView = typename Containers::ArrayView< PairIndexType, DeviceType >;
   using NeighborSearchPointerType = typename Pointers::SharedPointer< NeighborSearch, DeviceType >;

   NeighborsLoopParams( NeighborSearchPointerType& neighborSearch )
   : numberOfParticles( neighborSearch->getParticles()->getNumberOfParticles() ),
     gridSize( neighborSearch->getParticles()->getGridSize() ),
     view_firstLastCellParticle( neighborSearch->getCellFirstLastParticleList().getView() ) {}

   GlobalIndexType i;
   Containers::StaticVector< 2, GlobalIndexType > gridIndex;

   const GlobalIndexType numberOfParticles;
   const Containers::StaticVector< 2, GlobalIndexType > gridSize;
   const PairIndexArrayView view_firstLastCellParticle;
};

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
   using CellIndexArrayView = typename Containers::ArrayView< typename ParticleSystem::CellIndexType, DeviceType >;
   /* bucketing */
   using PairIndexType = Containers::StaticVector< 2, GlobalIndexType >;
   using PairIndexArrayType = Containers::Array< PairIndexType, DeviceType, GlobalIndexType >;
   using PairIndexArrayView = typename Containers::ArrayView< PairIndexType, DeviceType >;

   /* args */
   using NeighborsLoopParams = NeighborsLoopParams< NeighborSearch< ParticleConfig, ParticleSystem > >;

   /**
    * Constructors.
    */
   NeighborSearch( ParticlePointer& particles, GlobalIndexType cellCount )
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
   void
   particlesToCells();


   /**
    * Get particle system pointer.
    * TODO: Move.
    */
   ParticlePointer&
   getParticles()
   {
      return this->particles;
   }

   const ParticlePointer&
   getParticles() const
   {
      return this->particles;
   }

   /**
    * Loop over neighbor 2D.
    */
   template< typename Function, typename... FunctionArgs >
   __cuda_callable__
   void
   loopOverNeighbors( const GlobalIndexType i,
                      const GlobalIndexType& numberOfParticles,
                      const Containers::StaticVector< 2, GlobalIndexType >& gridIndex,
                      const Containers::StaticVector< 2, GlobalIndexType >& gridSize,
                      const PairIndexArrayView& view_firstLastCellParticle,
                      Function f, FunctionArgs... args );

   template< typename Function, typename... FunctionArgs >
   __cuda_callable__
   void
   loopOverNeighbors( const NeighborsLoopParams& params,
                      Function f, FunctionArgs... args );

   /**
    * Loop over neighbor 2D - loop over another set.
    *
    * This is basically the same as in the previous case, but this function
    * is able to perform interaction even for i == j, which is neccessary if
    * the iteration goes over another set.
    */
   template< typename Function, typename... FunctionArgs >
   __cuda_callable__
   void
   loopOverNeighborsAnotherSet( const GlobalIndexType i,
                                const GlobalIndexType& numberOfParticles,
                                const Containers::StaticVector< 2, GlobalIndexType >& gridIndex,
                                const Containers::StaticVector< 2, GlobalIndexType >& gridSize,
                                const PairIndexArrayView& view_firstLastCellParticle,
                                Function f, FunctionArgs... args );

   /**
    * Loop over neighbor 3D.
    */
   template< typename Function, typename... FunctionArgs >
   __cuda_callable__
   void
   loopOverNeighbors( const GlobalIndexType i,
                      const GlobalIndexType& numberOfParticles,
                      const Containers::StaticVector< 3, GlobalIndexType >& gridIndex,
                      const Containers::StaticVector< 3, GlobalIndexType >& gridSize,
                      const PairIndexArrayView& view_firstLastCellParticle,
                      Function f, FunctionArgs... args );

   /**
    * Loop over neighbor 3D - loop over another set.
    *
    * This is basically the same as in the previous case, but this function
    * is able to perform interaction even for i == j, which is neccessary if
    * the iteration goes over another set.
    */
   template< typename Function, typename... FunctionArgs >
   __cuda_callable__
   void
   loopOverNeighborsAnotherSet( const GlobalIndexType i,
                                const GlobalIndexType& numberOfParticles,
                                const Containers::StaticVector< 3, GlobalIndexType >& gridIndex,
                                const Containers::StaticVector< 3, GlobalIndexType >& gridSize,
                                const PairIndexArrayView& view_firstLastCellParticle,
                                Function f, FunctionArgs... args );

   //with vector
   template< typename Function, typename... FunctionArgs >
   __cuda_callable__
   void
   loopOverNeighborsBlocks( const GlobalIndexType i,
                            const GlobalIndexType& numberOfParticles,
                            const Containers::StaticVector< 3, GlobalIndexType >& gridIndex,
                            const Containers::StaticVector< 3, GlobalIndexType >& gridSize,
                            const PairIndexArrayView& view_firstLastCellParticle,
                            const CellIndexArrayView& view_particleCellIndex, //TODO: Remove
                            Function f, FunctionArgs... args );

   /**
    * Runs all necessary functions and fills neighbor list.
    */
   void
   searchForNeighbors();

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

