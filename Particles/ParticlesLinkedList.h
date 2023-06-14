#pragma once

#include <TNL/Algorithms/parallelFor.h>
#include <TNL/Pointers/SharedPointer.h>
#include <TNL/Algorithms/sort.h>

#include <limits> //UINT_MAX

#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#include <thrust/gather.h>

#include "ParticlesTraits.h"
#include "GenerateCellIndex.h"
#include "Particles.h"

namespace TNL {
namespace ParticleSystem {

template< typename ParticleSystem >
class NeighborsLoopParams
{
public:
   using DeviceType = typename ParticleSystem::DeviceType;
   using GlobalIndexType = typename ParticleSystem::GlobalIndexType;
   using PairIndexType = Containers::StaticVector< 2, GlobalIndexType >;
   using CellIndexArrayView = typename Containers::ArrayView< typename ParticleSystem::CellIndexType, DeviceType >;
   using PairIndexArrayView = typename Containers::ArrayView< PairIndexType, DeviceType >;
   using ParticlesPointerType = typename Pointers::SharedPointer< ParticleSystem, DeviceType >;
   using PointType = typename ParticleSystem::PointType;
   using CellIndexer = typename ParticleSystem::CellIndexer;
   using IndexVectorType = typename ParticleSystem::IndexVectorType;
   using RealType = typename ParticleSystem::RealType;

   NeighborsLoopParams( ParticlesPointerType& particles )
   : numberOfParticles( particles->getNumberOfParticles() ),
     gridSize( particles->getGridSize() ),
     gridOrigin( particles->getGridOrigin() ),
     searchRadius( particles->getSearchRadius() ),
     view_firstLastCellParticle( particles->getCellFirstLastParticleList().getView() ) {}

   GlobalIndexType i;
   Containers::StaticVector< 2, GlobalIndexType > gridIndex;

   const GlobalIndexType numberOfParticles;
   const Containers::StaticVector< 2, GlobalIndexType > gridSize;
   //const typename ParticleSystem::IndexVectorType gridSize;
   const PairIndexArrayView view_firstLastCellParticle;
   const PointType gridOrigin;
   const RealType searchRadius;
};

template < typename ParticleConfig, typename Device >
class ParticlesLinkedList : public Particles< ParticleConfig, Device > {
public:

   using DeviceType = Device;
   using Config = ParticleConfig;
   using ParticleTraitsType = ParticlesTraits< Config, DeviceType >;

   /* common */
   using GlobalIndexType = typename ParticleTraitsType::GlobalIndexType;
   using LocalIndexType = typename ParticleTraitsType::LocalIndexType;
   using RealType = typename ParticleTraitsType::RealType;
   using CellIndexType = typename ParticleTraitsType::CellIndexType;
   using CellIndexArrayType = typename ParticleTraitsType::CellIndexArrayType; //nn
   using NeighborsCountArrayType = typename ParticleTraitsType::NeighborsCountArrayType;
   using NeighborsArrayType = typename ParticleTraitsType::NeighborsArrayType;
   using NeighborListType = typename ParticleTraitsType::NeighborListType;

   using IndexArrayType = typename ParticleTraitsType::CellIndexArrayType; //TODO: Clean up the types.
   using IndexArrayTypePointer = typename Pointers::SharedPointer< IndexArrayType, Device >;

   using CellIndexer = typename Config::CellIndexerType;

   /* particle related */
   using PointType = typename ParticleTraitsType::PointType;
   using PointArrayType = typename ParticleTraitsType::PointArrayType; //nn

   /* grid related */
   using IndexVectorType = typename ParticleTraitsType::IndexVectorType;
   using GridType = typename ParticleTraitsType::GridType;
   using GridPointer = typename ParticleTraitsType::GridPointer;

   /* neighbor list related */
   using PairIndexType = Containers::StaticVector< 2, GlobalIndexType >;
   using PairIndexArrayType = Containers::Array< PairIndexType, DeviceType, GlobalIndexType >;
   using PairIndexArrayView = typename Containers::ArrayView< PairIndexType, DeviceType >;

   /* args */
   using NeighborsLoopParams = NeighborsLoopParams< ParticlesLinkedList< ParticleConfig, DeviceType > >;

   /**
    * Constructors.
    */
   ParticlesLinkedList() = default;

   ParticlesLinkedList( GlobalIndexType size, GlobalIndexType sizeAllocated, RealType radius, GlobalIndexType cellCount )
   : Particles< ParticleConfig, DeviceType >( size, sizeAllocated, radius ),
     particleCellInidices( sizeAllocated ),
     firstLastCellParticle( cellCount )
   {
      firstLastCellParticle = INT_MAX;
   }

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
    * Get list of first and last particle in cells.
    */
   const PairIndexArrayType&
   getCellFirstLastParticleList() const;

   PairIndexArrayType&
   getCellFirstLastParticleList();

   /**
    * Get cell index of given partile.
    */
   void
   computeParticleCellIndices();

   /**
    * Sort particles by its cell index.
    */
   void sortParticles();

   //NEIGHBOR SEARCH UTILITIES
   /**
    * Reset the list with first and last particle in cell.
    */
   void
   resetListWithIndices(); //protected?

   /**
    * Assign to each cell index of first contained particle.
    */
   void
   particlesToCells();


protected:

   //neighborsearch related;
   CellIndexArrayType particleCellInidices;
   PairIndexArrayType firstLastCellParticle;

};

} //namespace Particles
} //namespace TNL

#include "ParticlesLinkedList.hpp"

