#pragma once

#include <TNL/Logger.h>
#include <TNL/Containers/Vector.h>
#include <TNL/Algorithms/reduce.h>

#include "../Particles/ParticlesTraits.h"

#include "Fluid.h"
#include "Boundary.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Model >
class SPHSimpleFluid
{
public:

   using DeviceType = typename Model::DeviceType;
   using ParticleSystemType = typename Model::ParticlesType;; //Added due to measure tool.
   using ParticleSystem = typename Model::ParticlesType;; //Added due to measure tool.
   using ModelType = Model; //Added due to distributed simulation.

   using GlobalIndexType = typename ParticleSystem::GlobalIndexType;

   using ParticlePointer = typename Pointers::SharedPointer< ParticleSystem, DeviceType >;
   using ModelPointer = typename Pointers::SharedPointer< Model, DeviceType >;
   using IntegratorPointer = typename Pointers::SharedPointer< typename Model::Integrator, DeviceType >; // draft
   using SPHConfig = typename Model::SPHConfig;
   using IntegratorVariables = typename Model::IntegratorVariables;

   using FluidVariables = typename Model::FluidVariables;
   using Fluid = Fluid< ParticleSystem, SPHConfig, FluidVariables, IntegratorVariables >;
   using FluidPointer = Pointers::SharedPointer< Fluid, DeviceType >;

   using BoundaryVariables = typename Model::BoundaryVariables;
   using Boundary = Boundary< ParticleSystem, SPHConfig, BoundaryVariables, IntegratorVariables >;
   using BoundaryPointer = Pointers::SharedPointer< Boundary, DeviceType >;

   SPHSimpleFluid() = default;

   template< typename SPHSimpleFluidInit >
   SPHSimpleFluid( SPHSimpleFluidInit sphConfig )
   :  model(),
      integrator(),
      fluid( sphConfig.numberOfParticles, sphConfig.numberOfAllocatedParticles, sphConfig.searchRadius, sphConfig.numberOfGridCells ),
      boundary( sphConfig.numberOfBoundaryParticles, sphConfig.numberOfAllocatedBoundaryParticles, sphConfig.searchRadius, sphConfig.numberOfGridCells )
   {
      fluid->particles->setGridSize( sphConfig.gridSize );
      fluid->particles->setGridOrigin( sphConfig.gridOrigin );
      boundary->particles->setGridSize( sphConfig.gridSize );
      boundary->particles->setGridOrigin( sphConfig.gridOrigin );
   };

   /**
    * Perform neighbors search and fill neighborsList in Particle system variable.
    */
   void
   PerformNeighborSearch( GlobalIndexType step, TNL::Timer& timer_reset, TNL::Timer& timer_cellIndices, TNL::Timer& timer_sort, TNL::Timer& timer_toCells );

   /**
    * Perform interaction for all particles, i.e. for all types.
    */
   //delta-WCSPH models
   template< typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm, typename EOS, typename SPHState >
   void
   interact( SPHState& sphState );

   template< typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm, typename EOS, typename SPHState, typename TimeStepping >
   void
   interact( SPHState& sphState, TimeStepping& timeStepping );

   //RSPH models
   template< typename SPHKernelFunction, typename RiemannSolver, typename EOS, typename SPHState >
   void
   interact( SPHState& sphState );

   template< typename Writer >
   void
   save( const std::string& outputFilename, const int step, bool writeParticleCellIndex = false );

   void
   writeProlog( TNL::Logger& logger ) const noexcept;

//protected:

   FluidPointer fluid;
   BoundaryPointer boundary;

   ModelPointer model;

   IntegratorPointer integrator;

};

} // SPH
} // ParticleSystem
} // TNL

template< typename Model >
std::ostream&
operator<<( std::ostream& str, const SPH::SPHSimpleFluid< Model >& sphSimulation )
{
   TNL::Logger logger( 100, str );

   sphSimulation.writeProlog( logger );

   return str;
}


#include "SPH_impl.h"

