#pragma once

#include <TNL/Containers/Vector.h>
#include <TNL/Algorithms/reduce.h>

#include "../Particles/Particles.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Model, typename ParticleSystem, typename NeighborSearch >
class SPHSimulation
{
public:

   using DeviceType = typename ParticleSystem::Device;

   using LocalIndexType = typename ParticleSystem::LocalIndexType;
   using GlobalIndexType = typename ParticleSystem::GlobalIndexType;
   using PointType = typename ParticleSystem::PointType; //remove
   using RealType = typename ParticleSystem::RealType;

   using ParticlePointer = typename Pointers::SharedPointer< ParticleSystem, DeviceType >;
   using NeighborSearchPointer = typename Pointers::SharedPointer< NeighborSearch, DeviceType >;
   using ModelPointer = typename Pointers::SharedPointer< Model, DeviceType >;

   SPHSimulation() = default;

   SPHSimulation( GlobalIndexType size, RealType h, GlobalIndexType numberOfCells )
   : particles( size, h ), neighborSearch( particles, numberOfCells ), model( size, particles ) {};

   /**
    * Perform neighbors search and fill neighborsList in Particle system variable.
    */
   void PerformNeighborSearch( GlobalIndexType step );

   /**
    * Perform interaction for all particles, i.e. for all types.
    */
   template< typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm >
   void Interact();

   //protected: (or private?)

   ParticlePointer particles;
   NeighborSearchPointer neighborSearch;
   ModelPointer model;

};

} // SPH
} // ParticleSystem
} // TNL

#include "SPH_impl.h"

