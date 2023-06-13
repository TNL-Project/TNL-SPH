#pragma once

#include <TNL/Algorithms/parallelFor.h>
#include <TNL/Pointers/SharedPointer.h>
#include <TNL/Algorithms/sort.h>

#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#include <thrust/gather.h>

#include "ParticlesTraits.h"
#include "GenerateCellIndex.h"

namespace TNL {
namespace ParticleSystem {

template < typename ParticleConfig, typename DeviceType >
class Particles{
public:

   using Device = DeviceType;
   using Config = ParticleConfig;
   using ParticleTraitsType = ParticlesTraits< Config, DeviceType >;

   static constexpr int spaceDimension = 2;

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

   /**
    * Constructors.
    */
   Particles() = default;

   Particles( GlobalIndexType size, GlobalIndexType sizeAllocated )
   : numberOfParticles( size ), numberOfAllocatedParticles( sizeAllocated ), points( sizeAllocated ) { }

   Particles( GlobalIndexType size, GlobalIndexType sizeAllocated, RealType radius )
   : numberOfParticles( size ),
     numberOfAllocatedParticles( sizeAllocated ),
     points( sizeAllocated ),
     points_swap( sizeAllocated ),
     sortPermutations( sizeAllocated ),
     radius( radius ),
     particleCellInidices( sizeAllocated )
   {
      //grid->setSpaceSteps( { Config::searchRadius, Config::searchRadius } ); //removed
      //3dto grid->setDimensions( Config::gridXsize, Config::gridYsize );
      //grid->setOrigin( { Config::gridXbegin, Config::gridYbegin } ); //removed
      //neighborsList.setSegmentsSizes( sizeAllocated, Config::maxOfNeigborsPerParticle ); DeactivatedAtm
   }

   /* PARTICLE RELATED TOOLS */

   /**
    * Get dimension of particle system.
    */
   static constexpr int
   getParticleDimension();

   /**
    * Get search radius.
    */
   __cuda_callable__
   const RealType
   getSearchRadius() const;

   /**
    * Get number of particles in particle system.
    */
   __cuda_callable__
   GlobalIndexType
   getNumberOfParticles();

   __cuda_callable__
   const GlobalIndexType
   getNumberOfParticles() const;

   __cuda_callable__
   const GlobalIndexType
   getNumberOfAllocatedParticles() const;

   void
   setNumberOfParticles( GlobalIndexType newNumberOfParticles );

   /**
    * Get particle (i.e. point) positions.
    */
   const typename ParticleTraitsType::PointArrayType& // -> using..
   getPoints() const;

   typename ParticleTraitsType::PointArrayType& // -> using..
   getPoints();

   /**
    * Get position of given particle.
    */
   __cuda_callable__
   const PointType&
   getPoint( GlobalIndexType particleIndex ) const;

   __cuda_callable__
   PointType&
   getPoint( GlobalIndexType particleIndex );

   /**
    * Set position of given particle.
    */
   __cuda_callable__
   void
   setPoint( GlobalIndexType particleIndex, PointType point);

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

   const IndexArrayTypePointer&
   getSortPermutations() const;

   IndexArrayTypePointer&
   getSortPermutations();

   /**
    * Get cell index of given partile.
    */
   void
   computeParticleCellIndices();

   /**
    * Sort particles by its cell index.
    */
   void sortParticles();

   /* PARTICLE RELATED TEMP TOOLS */

   void generateRandomParticles(); //used only for tests

   /* GRID RELATED TOOLS */

   /**
    * Get and set grid dimension.
    * The grid here is just implicit.
    */
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

   const typename ParticleTraitsType::CellIndexArrayType& // -> using..
   getGridCellIndices() const;

   typename ParticleTraitsType::CellIndexArrayType& // -> using..
   getGridCellIndices();

   __cuda_callable__
   const CellIndexType&
   getGridCellIndex( GlobalIndexType cellIndex ) const;

   __cuda_callable__
   CellIndexType&
   getGridCellIndex( GlobalIndexType cellIndex );

   void computeGridCellIndices();

   /* general */
   void GetParticlesInformations();

   /* NEIGHBOR LIST RELATED TOOL */

   /**
    * Return list with neighbor particles.
    */
   const NeighborsArrayType& // -> using..
   getNeighborsList() const;

   NeighborsArrayType& // -> using..
   getNeighborsList();

   /**
    * Return jth neighbor of particle i.
    */
   __cuda_callable__
   const GlobalIndexType&
   getNeighbor( GlobalIndexType i, GlobalIndexType j ) const;

   __cuda_callable__
   GlobalIndexType&
   getNeighbor( GlobalIndexType i, GlobalIndexType j );

   /**
    * Return list with numbers of particles.
    */
   const NeighborsCountArrayType& // -> using..
   getNeighborsCountList() const;

   NeighborsCountArrayType& // -> using..
   getNeighborsCountList();

   /**
    * Return number of neighbors for particle i.
    */
   __cuda_callable__
   const LocalIndexType&
   getNeighborsCount( GlobalIndexType particleIndex ) const;

   __cuda_callable__
   LocalIndexType&
   getNeighborsCount( GlobalIndexType particleIndex );

   /**
    * Set j as neighbor for particle i.
    */
   __cuda_callable__
   void
   setNeighbor( GlobalIndexType i, GlobalIndexType j );

   /**
    * Remove all neighbors and clear the neighbor list.
    */
   void
   resetNeighborList();

   /**
    * Print/save neighbor whole neighbor list.
    */
   void
   saveNeighborList(std::string neigborListFile);

protected:

   /* particle related*/
   GlobalIndexType numberOfAllocatedParticles;
   GlobalIndexType numberOfParticles;
   GlobalIndexType gridSize;

   /* grid related */
   PointType gridOrigin; //atm we use only implicit grid.
   IndexVectorType gridDimension; //atm we use only implicit grid.

   GridPointer grid; // not used now

   RealType radius;

   NeighborsCountArrayType neighborsCount;
   NeighborsArrayType neighbors;

   NeighborListType neighborsList; //not used now

   PointArrayType points;
   PointArrayType points_swap; //avoid a inplace sort

   IndexArrayTypePointer sortPermutations;

   CellIndexArrayType particleCellInidices;

   CellIndexArrayType gridCellIndices;

};

} //namespace Particles
} //namespace TNL

#include "Particles_impl.h"

