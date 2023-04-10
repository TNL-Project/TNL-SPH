#pragma once

#include <TNL/Containers/Vector.h>
#include <TNL/Algorithms/reduce.h>
#include <memory> //shared_ptr

#include "../Particles/Particles.h"
#include "SPH.h"

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
   using OpenBoudaryPatchesPointerArray = Containers::Array< OpenBoudaryPatchPointer, DeviceType >;

   SPHOpenSystem() = default;

   //SPHOpenSystem( GlobalIndexType size, GlobalIndexType sizeAllocated,
   //      GlobalIndexType size_bound, GlobalIndexType sizeAllocated_bound,
   //      RealType h, GlobalIndexType numberOfCells, GlobalIndexType numberOfInlets )
   //: model(), integrator(), fluid( size, sizeAllocated, h, numberOfCells ),
   //  boundary( size_bound, sizeAllocated_bound, h, numberOfCells ){};

   SPHOpenSystem( SPHSimpleFluidConfig< typename ParticleSystem::Config > sphConfig )
   :  model(),
      integrator(),
      fluid( sphConfig.sizeFluid, sphConfig.sizeAllocatedFluid, sphConfig.searchRadius, sphConfig.gridNumberOfCells ),
      boundary( sphConfig.sizeBoundary, sphConfig.sizeAllocatedBoundary, sphConfig.searchRadius, sphConfig.gridNumberOfCells )
   {
      fluid->particles->setGridSize( sphConfig.gridSize );
      fluid->particles->setGridOrigin( sphConfig.gridOrigin );
      boundary->particles->setGridSize( sphConfig.gridSize );
      boundary->particles->setGridOrigin( sphConfig.gridOrigin );
   };

   /**
    * Perform neighbors search and fill neighborsList in Particle system variable.
    */
   void PerformNeighborSearch( GlobalIndexType step, TNL::Timer& timer_reset, TNL::Timer& timer_cellIndices, TNL::Timer& timer_sort, TNL::Timer& timer_toCells );

   /**
    * Perform interaction for all particles, i.e. for all types.
    */
   template< typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm, typename EOS >
   void Interact();

   /**
    * TODO: I don't like this.
    */
   void addOpenBoundaryPatch( GlobalIndexType size_buffer, GlobalIndexType sizeAllocated_buffer, RealType h, GlobalIndexType numberOfCells );

   void
   writeProlog( TNL::Logger& logger ) const noexcept;

//protected:

   FluidPointer fluid;
   BoundaryPointer boundary;

   //OpenBoudaryPatchesPointerArray openBoundaryPatches;
   std::vector< OpenBoudaryPatchPointer > openBoundaryPatches;

   ModelPointer model;

   IntegratorPointer integrator;

};

} // SPH
} // ParticleSystem
} // TNL

template< typename Model, typename ParticleSystem, typename NeighborSearch >
std::ostream&
operator<<( std::ostream& str, const SPH::SPHOpenSystem< Model, ParticleSystem, NeighborSearch >& sphSimulation )
{
   TNL::Logger logger( 100, str );

   sphSimulation.writeProlog( logger );

   return str;
}

#include "SPHOpen_impl.h"

