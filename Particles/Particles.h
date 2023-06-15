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

   using GlobalIndexType = typename ParticleTraitsType::GlobalIndexType;
   using LocalIndexType = typename ParticleTraitsType::LocalIndexType;
   using RealType = typename ParticleTraitsType::RealType;

   using IndexArrayType = typename ParticleTraitsType::CellIndexArrayType;
   using IndexArrayTypePointer = typename Pointers::SharedPointer< IndexArrayType, Device >;
   using PointType = typename ParticleTraitsType::PointType;
   using PointArrayType = typename ParticleTraitsType::PointArrayType;

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
     radius( radius ) { }

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

   /**
    * Set number of particles in particle system.
    */
   void
   setNumberOfParticles( GlobalIndexType newNumberOfParticles );

   /**
    * Get particle (i.e. point) positions.
    */
   const PointArrayType&
   getPoints() const;

   PointArrayType&
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
    * Get list of prermutation for particles reordering.
    */
   const IndexArrayTypePointer&
   getSortPermutations() const;

   IndexArrayTypePointer&
   getSortPermutations();

protected:

   //information about particle system
   GlobalIndexType numberOfAllocatedParticles;
   GlobalIndexType numberOfParticles;
   RealType radius;

   //actual points
   PointArrayType points;
   PointArrayType points_swap; //avoid a inplace sort
   IndexArrayTypePointer sortPermutations;

};

} //namespace Particles
} //namespace TNL

#include "Particles.hpp"

