#pragma once

#include <TNL/Containers/Vector.h>
#include <TNL/Algorithms/reduce.h>
#include <memory> //shared_ptr

#include "../Particles/Particles.h"
#include "OpenBoundaryBuffers.h"

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

   using SPHPointer = typename Pointers::SharedPointer< SPHOpenSystem, DeviceType >; // draft

   using OpenBoundaryPatch = OpenBoundaryBuffer_orthogonal< ParticleSystem, NeighborSearch, typename Model::SPHConfig >;
   using OpenBoudaryPatchPointer = Pointers::SharedPointer< OpenBoundaryPatch, DeviceType >;
   using OpenBoudaryPatchesPointerArray = Containers::Array< OpenBoudaryPatchPointer, DeviceType >;

   SPHOpenSystem() = default;

   SPHOpenSystem( GlobalIndexType size, GlobalIndexType sizeAllocated,
         GlobalIndexType size_bound, GlobalIndexType sizeAllocated_bound,
         GlobalIndexType size_buffer, GlobalIndexType sizeAllocated_buffer,
         RealType h, GlobalIndexType numberOfCells, GlobalIndexType numberOfInlets )
   : particles( size, sizeAllocated, h ), neighborSearch( particles, numberOfCells ),
     particles_bound( size_bound, sizeAllocated_bound, h ), neighborSearch_bound( particles_bound, numberOfCells ),
     model( sizeAllocated,
            particles,
            sizeAllocated_bound,
            particles_bound,
            sizeAllocated_buffer,
            openBoundaryPatch),
     integrator( model,
                 sizeAllocated,
                 sizeAllocated_bound,
                 sizeAllocated_buffer ),
     openBoundaryPatch( size_buffer,
                        sizeAllocated_buffer,
                        h,
                        numberOfCells ){};

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

   //OpenBoudaryPatchesPointer openBoundaryPatches;
   OpenBoudaryPatchPointer openBoundaryPatch;

   OpenBoudaryPatchesPointerArray openBoundaryPatches;

   ModelPointer model;

   IntegratorPointer integrator;


   //using SPHPointer = typename Pointers::SharedPointer< SPHOpenSystem, DeviceType >; // draft
   //using SPHPointer = typename Pointers::DevicePointer< SPHOpenSystem, DeviceType >; // draft
   //using SPHPointer = std::shared_ptr< std::add_const_t< SPHOpenSystem > >;

   //SPHOpenSystem* self()
   //{ return this; }

   //const SPHOpenSystem* self() const
   //{ return this; }


   //SPHOpenSystem& self()
   //{ return *this; }

   //const SPHOpenSystem& self() const
   //{ return *this; }

};

} // SPH
} // ParticleSystem
} // TNL

#include "SPHOpen_impl.h"

