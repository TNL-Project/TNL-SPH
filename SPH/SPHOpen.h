#pragma once

#include <TNL/Containers/Vector.h>
#include <TNL/Algorithms/reduce.h>

#include "../Particles/Particles.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Model, typename ParticleSystem, typename NeighborSearch >
class SPHOpenSystem
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

   using IntegratorPointer = typename Pointers::SharedPointer< typename Model::Integrator, DeviceType >; // draft

   SPHOpenSystem() = default;

   SPHOpenSystem( GlobalIndexType size, GlobalIndexType sizeAllocated,
         GlobalIndexType size_bound, GlobalIndexType sizeAllocated_bound,
         GlobalIndexType size_buffer, GlobalIndexType sizeAllocated_buffer,
         RealType h, GlobalIndexType numberOfCells )
   : particles( size, sizeAllocated, h ), neighborSearch( particles, numberOfCells ),
     particles_bound( size_bound, sizeAllocated_bound, h ), neighborSearch_bound( particles_bound, numberOfCells ),
     particles_buffer( size_buffer, sizeAllocated_buffer, h ), neighborSearch_buffer( particles_bound, numberOfCells ),
     model( sizeAllocated, particles, sizeAllocated_bound, particles_bound, sizeAllocated_buffer, particles_buffer ), integrator( model, sizeAllocated, sizeAllocated_bound, sizeAllocated_buffer ) {};


   /**
    * Perform neighbors search and fill neighborsList in Particle system variable.
    */
   void PerformNeighborSearch( GlobalIndexType step, TNL::Timer& timer_reset, TNL::Timer& timer_cellIndices, TNL::Timer& timer_sort, TNL::Timer& timer_toCells );

   /**
    * Perform interaction for all particles, i.e. for all types.
    */
   template< typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm >
   void Interact();

   template< typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm, typename EOS >
   void InteractModel();

   template< typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm, typename EOS, typename RiemannSolver >
   void InteractModel();

   //protected: (or private?)

   ParticlePointer particles;
   ParticlePointer particles_bound;
   ParticlePointer particles_buffer;

   NeighborSearchPointer neighborSearch;
   NeighborSearchPointer neighborSearch_bound;
   NeighborSearchPointer neighborSearch_buffer;

   ModelPointer model;

   IntegratorPointer integrator;

};

} // SPH
} // ParticleSystem
} // TNL

#include "SPHOpen_impl.h"

