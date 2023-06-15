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
   using BaseType = Particles< ParticleConfig, Device >;

   /* common */
   using GlobalIndexType = typename BaseType::GlobalIndexType;
   using LocalIndexType = typename BaseType::LocalIndexType;
   using RealType = typename BaseType::RealType;
   using PointType = typename BaseType::PointType;
   using PointArrayType = typename BaseType::PointArrayType;

   /* neighbor search related */
   using CellIndexer = typename Config::CellIndexerType;

   using IndexVectorType = typename ParticleTraitsType::IndexVectorType;
   using CellIndexType = typename ParticleTraitsType::CellIndexType;
   using CellIndexArrayType = typename ParticleTraitsType::CellIndexArrayType;
   using PairIndexType = typename ParticleTraitsType::PairIndexType;
   using PairIndexArrayType = typename ParticleTraitsType::PairIndexArrayType;

   /* args for neighbors loop */
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

   __cuda_callable__ //TODO: Comment.
   const IndexVectorType
   getGridSize() const;

   void
   setGridSize( IndexVectorType gridSize );

   /**
    * Get and set grid origin.
    * The grid here is just implicit.
    */
   const PointType
   getGridOrigin() const;

   void
   setGridOrigin( PointType gridOrigin );

   /**
    * Get particle cell indices.
    */
   const CellIndexArrayType& // -> using..
   getParticleCellIndices() const;

   CellIndexArrayType& // -> using..
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

   //related to implicit grid
   PointType gridOrigin;
   IndexVectorType gridDimension;

   //related to search for neighbors
   CellIndexArrayType particleCellInidices;
   PairIndexArrayType firstLastCellParticle;

};

} //namespace Particles
} //namespace TNL

#include "ParticlesLinkedList.hpp"

