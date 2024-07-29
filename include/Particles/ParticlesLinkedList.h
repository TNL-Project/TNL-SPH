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
#include "neighborSearchLoop.h"

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
     gridSize( particles->getGridDimensionsWithOverlap() ),
     gridOrigin( particles->getGridOriginWithOverlap() ),
     searchRadius( particles->getSearchRadius() ),
     view_firstLastCellParticle( particles->getCellFirstLastParticleList().getView() ) {}

   const GlobalIndexType numberOfParticles;
   const IndexVectorType gridSize;
   const PointType gridOrigin;
   const RealType searchRadius;

   const PairIndexArrayView view_firstLastCellParticle;
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
   using NeighborsLoop = NeighborsLoopCellLinkedList;
   using NeighborsLoopAnotherSet = NeighborsLoopCellLinkedListAnotherSet;
   using NeighborsLoopParams = NeighborsLoopParams< ParticlesLinkedList< ParticleConfig, DeviceType > >;

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
    * \brief Set size of the particle system, i. e. number of
    * allocated particles.
    *
    * \param size Number of allocated particles.
    */
   void
   setSize( const GlobalIndexType& size );

   void
   setGridDimensions( const IndexVectorType& dimensions );

   /**
    * \brief Set width of overlap expressed in number of cells.
    *
    * \param Integer expressing the overlap width in number of cells
    */
   void
   setOverlapWidth( const GlobalIndexType width );

   /**
    * Get particle cell indices.
    */
   const CellIndexArrayType&
   getParticleCellIndices() const;

   CellIndexArrayType&
   getParticleCellIndices();

   /**
    * Get list of first and last particle in cells.
    */
   const PairIndexArrayType&
   getCellFirstLastParticleList() const;

   PairIndexArrayType&
   getCellFirstLastParticleList();

   /**
    * Reset the list with first and last particle in cell.
    */
   void
   resetListWithIndices();

   /**
    * Get cell index of given partile - with enabled domain decomposition.
    */
   template< typename UseWithDomainDecomposition = typename Config::UseWithDomainDecomposition,
             std::enable_if_t< UseWithDomainDecomposition::value, bool > Enabled = true >
   void
   computeParticleCellIndices();

   /**
    * Get cell index of given partile - without enabled domain decomposition.
    */
   template< typename UseWithDomainDecomposition = typename Config::UseWithDomainDecomposition,
             std::enable_if_t< !UseWithDomainDecomposition::value, bool > Enabled = true >
   void
   computeParticleCellIndices();

   /**
    * Sort particles by its cell index.
    */
   void
   sortParticles();

   /**
    * Assign to each cell index of first contained particle.
    */
   void
   particlesToCells();

   /**
    * FIXME: Here or in base?
    * Start remove procedure for all particles out of interior region.
    */
   void
   removeParitclesOutOfDomain();

   /**
    * Run all procedures required to perform neighbor search.
    */
   void
   searchForNeighbors();

   void
   writeProlog( TNL::Logger& logge ) const noexcept;

protected:

   /**
    * Array with sice of allocated particles indices of corresponding cells
    * computed based on particle position.
    */
   CellIndexArrayType particleCellInidices;

   /**
    * Array with size of number of cells pairs storing index of first and last
    * contained in each cell. We assume that particles are continuously sorted.
    */
   PairIndexArrayType firstLastCellParticle;

};

} //namespace Particles
} //namespace TNL

#include "ParticlesLinkedList.hpp"

