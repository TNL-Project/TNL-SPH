#pragma once

#include <TNL/Containers/Vector.h>
#include <TNL/Algorithms/reduce.h>

#include "../Particles/Particles.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Model, typename ParticleSystem, typename NeighborSearch >
class SPHSimpleFluid
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

   SPHSimpleFluid() = default;

   SPHSimpleFluid( GlobalIndexType size, GlobalIndexType size_bound, RealType h, GlobalIndexType numberOfCells )
   : particles( size, h ), neighborSearch( particles, numberOfCells ), model( size, particles, size_bound, particles_bound ),
     particles_bound( size_bound, h ), neighborSearch_bound( particles_bound, numberOfCells ),
     integrator( model, size, size_bound ) {};


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
   NeighborSearchPointer neighborSearch;
   NeighborSearchPointer neighborSearch_bound;

   ModelPointer model;

   IntegratorPointer integrator;

};

} // SPH
} // ParticleSystem
} // TNL

#include "SPH_impl.h"

