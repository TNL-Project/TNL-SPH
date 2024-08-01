#pragma once

#include <TNL/Algorithms/parallelFor.h>
#include <TNL/Pointers/SharedPointer.h>
#include <TNL/Algorithms/sort.h>

#include <limits> //UINT_MAX

#include "../SPH/shared/thrustExecPolicySelector.h"
#include <string_view>
#include <thrust/sort.h>
#include <thrust/gather.h>
#include <type_traits>

#include "ParticlesTraits.h"
#include "CellIndexer.h"
#include "ParticlesLinkedList.h"
#include "neighborSearchLoopCLLWithList.h"

namespace TNL {
namespace ParticleSystem {

template< typename ParticleSystem >
class NeighborsLoopParamsCLLWithList
{
public:

   using DeviceType = typename ParticleSystem::DeviceType;
   using GlobalIndexType = typename ParticleSystem::GlobalIndexType;
   using PointType = typename ParticleSystem::PointType;
   using NeighborListView = typename ParticleSystem::NeighborListType::ViewType;
   using IndexArrayView = typename Containers::ArrayView< GlobalIndexType, DeviceType >;
   using ParticlesPointerType = typename Pointers::SharedPointer< ParticleSystem, DeviceType >;

   NeighborsLoopParamsCLLWithList( ParticlesPointerType& particles )
   : neighborsCountLimit( particles->getNeighborsCountLimit() ),
     neighborListStorageView( particles->getNeighborListStorage().getView() ) {}

   GlobalIndexType neighborsCountLimit;
   IndexArrayView neighborListStorageView;
};

template < typename ParticleConfig, typename Device >
class ParticlesLinkedListWithList : public ParticlesLinkedList< ParticleConfig, Device > {
public:

   using DeviceType = Device;
   using Config = ParticleConfig;
   using ParticleTraitsType = ParticlesTraits< Config, DeviceType >;
   using BaseType = ParticlesLinkedList< ParticleConfig, Device >;

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

   using IndexArrayType = typename ParticleTraitsType::IndexArrayType;
   using NeighborListType = typename ParticleTraitsType::NeighborListType;

   /* args for neighbors loop */
   using NeighborsLoopParams = NeighborsLoopParamsCLLWithList< ParticlesLinkedListWithList< ParticleConfig, DeviceType > >;
   using NeighborsLoop = NeighborsLoopCellLinkedListWithList;
   using NeighborsLoopAnotherSet = NeighborsLoopCellLinkedListWithListAnitherSet;

   /**
    * Constructors.
    */
   //ParticlesLinkedList() : particleCellInidices( 0 ), firstLastCellParticle( 0 ) {}

   static std::string
   writeModelType()
   {
      return "TNL::ParticleSystem::ParticlesCellLinkedListWithList";
   }

   void
   setSize( const GlobalIndexType& size );

   void
   setNeighborsCountLimit( const GlobalIndexType& limit );

   const GlobalIndexType
   getNeighborsCountLimit() const;

   /**
    * Get neighbor list.
    */
   const NeighborListType&
   getNeighborList() const;

   NeighborListType&
   getNeighborList();

   /**
    * Get neighbor list.
    */
   const IndexArrayType&
   getNeighborListStorage() const;

   IndexArrayType&
   getNeighborListStorage();

   /**
    * Run all procedures required to perform neighbor search.
    */
   void
   buildParticleList();

   /**
    * Run all procedures required to perform neighbor search.
    */
   void
   searchForNeighbors();

   //void
   //writeProlog( TNL::Logger& logge ) const noexcept;

protected:

   /**
    * Defines maximum number of possible neighbors.
    */
   GlobalIndexType neighborsCountLimit = 100;

   /**
    * Sparse format explicitly storing particle neighbors. Particle pairs are stored
    * in form of sparse matrices.
    */
   NeighborListType neighborList;

   /**
    * Storage for the actual neighbor list.
    */
   IndexArrayType neighborListStorage;

};

} //namespace Particles
} //namespace TNL

#include "ParticlesLinkedListWithList.hpp"

