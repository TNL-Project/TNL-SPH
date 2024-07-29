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
#include "Particles.h"

namespace TNL {
namespace ParticleSystem {

template< typename ParticleSystem >
class NeighborsLoopParams
{
public:

   using DeviceType = typename ParticleSystem::DeviceType;
   using GlobalIndexType = typename ParticleSystem::GlobalIndexType;
   using RealType = typename ParticleSystem::RealType;
   using PointType = typename ParticleSystem::PointType;
   using IndexVectorType = typename ParticleSystem::IndexVectorType;
   using PairIndexType = typename ParticleSystem::PairIndexType;
   using PairIndexArrayView = typename Containers::ArrayView< PairIndexType, DeviceType >;
   using ParticlesPointerType = typename Pointers::SharedPointer< ParticleSystem, DeviceType >;

   using CellIndexer = typename ParticleSystem::CellIndexer;

   NeighborsLoopParams( ParticlesPointerType& particles )
   : numberOfParticles( particles->getNumberOfParticles() ),
     gridSize( particles->getGridDimensions() ),
     gridOrigin( particles->getGridOrigin() ),
     searchRadius( particles->getSearchRadius() ),
     view_firstLastCellParticle( particles->getCellFirstLastParticleList().getView() ) {}

   const GlobalIndexType numberOfParticles;
   const IndexVectorType gridSize;
   const PointType gridOrigin;
   const RealType searchRadius;

   const PairIndexArrayView view_firstLastCellParticle;
};

template < typename ParticleConfig, typename Device >
class ParticlesLinkedListWithList : public Particles< ParticleConfig, Device > {
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

   using NeighborListType = typename ParticleTraitsType::NeighborListType;

   /* args for neighbors loop */
   using NeighborsLoopParams = NeighborsLoopParams< ParticlesLinkedListWithList< ParticleConfig, DeviceType > >;

   /**
    * Constructors.
    */
   ParticlesLinkedList() : particleCellInidices( 0 ), firstLastCellParticle( 0 ) {}

   static std::string
   writeModelType()
   {
      return "TNL::ParticleSystem::ParticlesLinkedList";
   }

   /**
    * Run all procedures required to perform neighbor search.
    */
   void
   buildParticleList();

   void
   writeProlog( TNL::Logger& logge ) const noexcept;

protected:

   /**
    * Defines maximum number of possible neighbors.
    */
   GlobalIndexType neighborsCountLimit;

   /**
    * Sparse format explicitly storing particle neighbors. Particle pairs are stored
    * in form of sparse matrices.
    */
   NeighborListType neighborList;

};

} //namespace Particles
} //namespace TNL

#include "ParticlesLinkedList.hpp"

