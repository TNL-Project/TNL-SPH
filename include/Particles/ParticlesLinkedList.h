#pragma once

#include <TNL/Algorithms/parallelFor.h>
#include <TNL/Pointers/SharedPointer.h>
#include <TNL/Algorithms/sort.h>

#include <limits> //UINT_MAX

#include "../SPH/shared/thrustExecPolicySelector.h"
#include <string_view>
#include <thrust/sort.h>
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
   using RealType = typename ParticleSystem::RealType;
   using PointType = typename ParticleSystem::PointType;
   using IndexVectorType = typename ParticleSystem::IndexVectorType;
   using PairIndexType = typename ParticleSystem::PairIndexType;
   using PairIndexArrayView = typename Containers::ArrayView< PairIndexType, DeviceType >;
   using ParticlesPointerType = typename Pointers::SharedPointer< ParticleSystem, DeviceType >;

   using CellIndexer = typename ParticleSystem::CellIndexer;

   NeighborsLoopParams( ParticlesPointerType& particles )
   : numberOfParticles( particles->getLastActiveParticle() + 1 ),
     gridSize( particles->getGridSize() ),
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
   ParticlesLinkedList() : particleCellInidices( 0 ), firstActiveParticle( 0 ) {}

   ParticlesLinkedList( GlobalIndexType size, GlobalIndexType sizeAllocated, RealType radius, GlobalIndexType cellCount )
   : Particles< ParticleConfig, DeviceType >( size, sizeAllocated, radius ),
     firstActiveParticle( 0 ),
     lastActiveParticle( size - 1 ),
     particleCellInidices( sizeAllocated ),
     firstLastCellParticle( cellCount )
   {
      firstLastCellParticle = INT_MAX;
   }

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

   /**
    * \brief Get index of first active particle in particle set.
    * For simulation with varying number of particles, particle array is
    * allocated with bigger size then is the actual number of particles of
    * active particles. Within this approach, first particle doesn't have to be
    * necessary the particle at position zero.
    */
   const GlobalIndexType
   getFirstActiveParticle() const;

   /**
    * \brief Get and set index of last active particle in particle set.
    * For simulation with varying number of particles, particle array is
    * allocated with bigger size then is the actual number of particles of
    * active particles. Within this approach, first particle doesn't have to be
    * necessary the particle at position zero.
    */
   const GlobalIndexType
   getLastActiveParticle() const;

   /**
    * \brief Set and set index of first active particle in particle set.
    * For simulation with varying number of particles, particle array is
    * allocated with bigger size then is the actual number of particles of
    * active particles. Within this approach, first particle doesn't have to be
    * necessary the particle at position zero.
    *
    * \param firstActiveParticle Position of first active particle.
    */
   void
   setFirstActiveParticle( GlobalIndexType firstActiveParticle );

   /**
    * \brief Set and set index of last active particle in particle set.
    * For simulation with varying number of particles, particle array is
    * allocated with bigger size then is the actual number of particles of
    * active particles. Within this approach, first particle doesn't have to be
    * necessary the particle at position zero.
    *
    * \param lastActiveParticle Position of last active particle.
    */
   void
   setLastActiveParticle( GlobalIndexType lastActiveParticle );

   /**
    * \brief Returns dimensions of the implicit linked list grid.
    *
    * \return Coordinate vector with number of edges along each axis.
    */
   const IndexVectorType
   getGridSize() const;

   /**
    * \brief Set dimensions of the implicit linked list grid.
    *
    * \param gridSize grid dimensions given in a form of coordinate vector.
    */
   void
   setGridSize( IndexVectorType gridSize );

   /**
    * \brief Returns origin of the implicit linked list grid.
    *
    * \return the origin of the grid.
    */
   const PointType
   getGridOrigin() const;

   /**
    * \brief Set origin of the implicit linked list grid.
    *
    * \param the origin of the grid.
    */
   void
   setGridOrigin( PointType gridOrigin );

   //TODO: Following lines need to be think through, added due to overlaps
   const PointType
   getGridInteriorOrigin() const;

   const IndexVectorType
   getGridInteriorDimension() const;

   void
   setGridInteriorOrigin( PointType gridInteriorOrigin );

   void
   setGridInteriorDimension( IndexVectorType gridInteriorDimension );

   /**
    * Get particle cell indices.
    */
   const CellIndexArrayType&
   getParticleCellIndices() const;

   CellIndexArrayType&
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
    * \brief Test if is particle with given index located inside domain
    *
    * \param index of particle to test
    *
    * \return True or False based on the presence of particles in given domain.
    */
   __cuda_callable__
   bool
   isInsideDomain( const GlobalIndexType& index );

   /**
    * \brief Test if is particle with given index located inside domain
    *
    * \param position of particle to test
    *
    * \return True or False based on the presence of particles in given domain.
    */
   __cuda_callable__
   bool
   isInsideDomain( const PointType& point );

   /**
    */
   const GlobalIndexType
   getNumberOfParticlesToRemove() const;

   /**
    */
   void
   setNumberOfParticlesToRemove( GlobalIndexType removeCount );

   /**
    * \brief Execute function \e f in parallel for all particles
    *
    * The function \e f is executed as `f(i)`, where `GlobalIndexType i` is the global index of the
    * particle to be processed. The particle set itself is not passed to the function `f`, it is the user's
    * responsibility to ensure proper access to the mesh if needed, e.g. by the means of lambda capture
    * and/or using a \ref TNL::Pointers::SharedPointer "SharedPointer".
    */
   template< typename Device2 = DeviceType, typename Func >
   void
   forAll( Func f ) const;


   /**
    * Sort particles by its cell index.
    */
   void
   sortParticles();

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

   /**
    * Find first and last particle in grid column.
    */
   PairIndexType
   getFirstLastParticleInColumnOfCells( const GlobalIndexType& gridColumn );

   /**
    * Find first and last particle in grid block of column.
    */
   PairIndexType
   getFirstLastParticleInBlockOfCells( const GlobalIndexType& gridBlock );

   void
   writeProlog( TNL::Logger& logge ) const noexcept;

protected:

   //related to implicit grid
   PointType gridOrigin;
   IndexVectorType gridDimension;

   GlobalIndexType firstActiveParticle;
   GlobalIndexType lastActiveParticle;

   //related to search for neighbors
   CellIndexArrayType particleCellInidices;
   PairIndexArrayType firstLastCellParticle;

   //number of particles to remove after sort
   GlobalIndexType numberOfParticlesToRemove;

   //number of additional cell layers to from domain overlap
   GlobalIndexType overlapWidthInCells = 1;

   //NOTE: This is probably not necersary and should be removedd
   PointType gridInteriorOrigin;
   IndexVectorType gridInteriorDimension;


};

} //namespace Particles
} //namespace TNL

#include "ParticlesLinkedList.hpp"

