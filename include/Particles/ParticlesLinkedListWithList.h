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

   NeighborsLoopParamsCLLWithList() = default;

   NeighborsLoopParamsCLLWithList( ParticlesPointerType& particles )
   : particleSetLabel( particles->getParticleSetLabel() ),
     numberOfParticles( particles->getNumberOfParticles() ),
     neighborsCountLimit( particles->getNeighborsCountLimit() ),
     neighborListStorageView( particles->getNeighborListStorage().getView() ) {}

   int particleSetLabel;
   GlobalIndexType neighborsCountLimit;
   IndexArrayView neighborListStorageView;

   GlobalIndexType numberOfParticles;
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
   using NeighborListLayout = NeighborListLayouts::NeighborMajorLinear;

   /**
    * Constructors.
    */
   //ParticlesLinkedList() : particleCellInidices( 0 ), firstLastCellParticle( 0 ) {}

   static std::string
   writeModelType()
   {
      return "TNL::ParticleSystem::ParticlesCellLinkedListWithList";
   }

   static constexpr bool
   specifySearchedSetExplicitly()
   {
      return true;
   }

   void
   setSize( const GlobalIndexType& size );

   void
   setNeighborsCountLimit( const GlobalIndexType& limit );

   const GlobalIndexType
   getNeighborsCountLimit() const;

   void
   setParticleSetLabel( const int& label );

   const int
   getParticleSetLabel() const;

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

   NeighborsLoopParams
   getCLLwLSearchToken();

   template< typename ParticlesPointerType >
   NeighborsLoopParams
   getCLLwLSearchToken( ParticlesPointerType& particlesToSearch );

   NeighborsLoopParams
   getSearchToken();

   template< typename ParticlesPointerType >
   NeighborsLoopParams
   getSearchToken( ParticlesPointerType& particlesToSearch );

   /**
    * \brief Reset list with neighbor.
    */
   void
   resetNeighborList();

   /**
    * \bief Run all the procedures required before collecting the neighbors
    */
   void
   makeSetSearchable();

   /**
    * Run all procedures required to perform neighbor search.
    *
    * Using Neihgbor Major Linear Layout
    */
   template< typename Layout = NeighborListLayout,
             std::enable_if_t< std::is_same_v< Layout, NeighborListLayouts::NeighborMajorLinear >, bool > Enabled = true >
   void
   buildParticleList();

   template< typename ParticlesPointerType,
             typename Layout = NeighborListLayout,
             std::enable_if_t< std::is_same_v< Layout, NeighborListLayouts::NeighborMajorLinear >, bool > Enabled = true >
   void
   buildParticleList( ParticlesPointerType& particlesToSearch );

   template< typename ParticlesPointerType,
             typename Layout = NeighborListLayout,
             std::enable_if_t< std::is_same_v< Layout, NeighborListLayouts::NeighborMajorLinear >, bool > Enabled = true >
   void
   addToParticleList( ParticlesPointerType& particlesToSearch );

   /**
    * Run all procedures required to perform neighbor search.
    *
    * Using Particles Major Linear Layout
    */
   template< typename Layout = NeighborListLayout,
             std::enable_if_t< std::is_same_v< Layout, NeighborListLayouts::ParticleMajorLinear >, bool > Enabled = true >
   void
   buildParticleList();

   template< typename ParticlesPointerType,
             typename Layout = NeighborListLayout,
             std::enable_if_t< std::is_same_v< Layout, NeighborListLayouts::ParticleMajorLinear >, bool > Enabled = true >
   void
   buildParticleList( ParticlesPointerType& particlesToSearch );

   template< typename ParticlesPointerType,
             typename Layout = NeighborListLayout,
             std::enable_if_t< std::is_same_v< Layout, NeighborListLayouts::ParticleMajorLinear >, bool > Enabled = true >
   void
   addToParticleList( ParticlesPointerType& particlesToSearch );

   /**
    * Run all procedures required to perform neighbor search.
    */
   void
   searchForNeighbors();

   void
   writeProlog( TNL::Logger& logger ) const noexcept;

protected:

   /**
    * Defines maximum number of possible neighbors.
    */
   GlobalIndexType neighborsCountLimit = 256;

   /**
    * Sparse format explicitly storing particle neighbors. Particle pairs are stored
    * in form of sparse matrices.
    */
   NeighborListType neighborList;

   /**
    * Storage for the actual neighbor list.
    */
   IndexArrayType neighborListStorage;

   /**
    * Particle set number.
    *
    * This variable is used to add neighbors from different particle sets.
    * All particle sets have to be labeled in ascending order integers starting
    * from zero. This allows add and iterate over neighbors from another set.
    */
   int particleSetLabel = 0;

};

} //namespace Particles
} //namespace TNL

#include "ParticlesLinkedListWithList.hpp"

