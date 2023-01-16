#pragma once

#include <TNL/Containers/Vector.h>
#include <TNL/Algorithms/reduce.h>
#include <memory> //shared_ptr

#include "../Particles/Particles.h"


#include "Fluid.h"
#include "Boundary.h"
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
   using SPHConfig = typename Model::SPHConfig;
   using Variables = typename Model::ModelVariables;
   using IntegratorVariables = typename Model::IntegratorVariables;

   using Fluid = Fluid< ParticleSystem, NeighborSearch, SPHConfig, Variables, IntegratorVariables >;
   using FluidPointer = Pointers::SharedPointer< Fluid, DeviceType >;

   using Boundary = Boundary< ParticleSystem, NeighborSearch, SPHConfig, Variables, IntegratorVariables >;
   using BoundaryPointer = Pointers::SharedPointer< Boundary, DeviceType >;

   using OpenBoundaryPatch = OpenBoundaryBuffer_orthogonal< ParticleSystem, NeighborSearch, SPHConfig, Variables >;
   using OpenBoudaryPatchPointer = Pointers::SharedPointer< OpenBoundaryPatch, DeviceType >;


   //using OpenBoudaryPatchesPointerArray = Containers::Array< OpenBoudaryPatchPointer, DeviceType >;

   SPHOpenSystem() = default;

   SPHOpenSystem( GlobalIndexType size, GlobalIndexType sizeAllocated,
         GlobalIndexType size_bound, GlobalIndexType sizeAllocated_bound,
         GlobalIndexType size_buffer, GlobalIndexType sizeAllocated_buffer,
         RealType h, GlobalIndexType numberOfCells, GlobalIndexType numberOfInlets )
   : model(
         sizeAllocated,
         sizeAllocated_bound,
         sizeAllocated_buffer ),
     integrator(
         model,
         sizeAllocated,
         sizeAllocated_bound,
         sizeAllocated_buffer ),
     openBoundaryPatch(
         size_buffer,
         sizeAllocated_buffer,
         h,
         numberOfCells ),
      fluid(
         size,
         sizeAllocated,
         h,
         numberOfCells ),
      boundary(
         size_bound,
         sizeAllocated_bound,
         h,
         numberOfCells ){};

   /**
    * Perform neighbors search and fill neighborsList in Particle system variable.
    */
   void PerformNeighborSearch( GlobalIndexType step, TNL::Timer& timer_reset, TNL::Timer& timer_cellIndices, TNL::Timer& timer_sort, TNL::Timer& timer_toCells );

   /**
    * Perform interaction for all particles, i.e. for all types.
    */
   template< typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm, typename EOS >
   void Interact();

//protected:

   FluidPointer fluid;
   BoundaryPointer boundary;
   OpenBoudaryPatchPointer openBoundaryPatch;
   //OpenBoudaryPatchesPointerArray openBoundaryPatches;

   ModelPointer model;

   IntegratorPointer integrator;

};

} // SPH
} // ParticleSystem
} // TNL

#include "SPHOpen_impl.h"

