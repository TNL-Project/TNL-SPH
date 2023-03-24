#pragma once

#include <TNL/Containers/Vector.h>
#include <TNL/Algorithms/reduce.h>

#include "../Particles/Particles.h"
#include "../Particles/ParticlesTraits.h"

#include "Fluid.h"
#include "Boundary.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename ParticleConfig >
struct SPHSimpleFluidConfig
{
   using ParticleTraitsType = ParticlesTraits< ParticleConfig >;
   using GlobalIndexType = typename ParticleTraitsType::GlobalIndexType;
   using IndexVectorType = typename ParticleTraitsType::IndexVectorType;
   using RealType = typename ParticleTraitsType::RealType;
   using PointType = typename ParticleTraitsType::PointType;

   GlobalIndexType sizeFluid;
   GlobalIndexType sizeAllocatedFluid;

   GlobalIndexType sizeBoundary;
   GlobalIndexType sizeAllocatedBoundary;

   IndexVectorType gridSize;
   PointType gridOrigin;

   RealType searchRadius;

   template< typename Config >
   void loadParameters()
   {
      this->sizeFluid = Config::numberOfParticles;
      this->sizeAllocatedFluid = Config::numberOfAllocatedParticles;

      this->sizeBoundary = Config::numberOfBoundaryParticles;
      this->sizeAllocatedBoundary = Config::numberOfAllocatedBoundaryParticles;


      if constexpr( ParticleConfig::spaceDimension == 2 )
      {
         this->gridSize = { Config::gridXsize, Config::gridYsize };
         this->gridOrigin = { Config::gridXorigin, Config::gridYorigin };
      }
      else if constexpr( ParticleConfig::spaceDimension == 3 )
      {
         this->gridSize = { Config::gridXsize, Config::gridYsize, Config::gridZsize };
         this->gridOrigin = { Config::gridXorigin, Config::gridYorign, Config::gridZsize };
      }
      //else
      //   std::runtime_error( "SPHSimpleFluidConfig error: invalid number of dimension." );

      this->searchRadius = Config::searchRadius;
   }
};

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
   using SPHConfig = typename Model::SPHConfig;
   using Variables = typename Model::ModelVariables;
   using IntegratorVariables = typename Model::IntegratorVariables;

   using Fluid = Fluid< ParticleSystem, NeighborSearch, SPHConfig, Variables, IntegratorVariables >;
   using FluidPointer = Pointers::SharedPointer< Fluid, DeviceType >;

   using Boundary = Boundary< ParticleSystem, NeighborSearch, SPHConfig, Variables, IntegratorVariables >;
   using BoundaryPointer = Pointers::SharedPointer< Boundary, DeviceType >;

   SPHSimpleFluid() = default;

   SPHSimpleFluid( GlobalIndexType size, GlobalIndexType sizeAllocated,
         GlobalIndexType size_bound, GlobalIndexType sizeAllocated_bound,
         RealType h, GlobalIndexType numberOfCells, GlobalIndexType numberOfInlets )
   :  model(), integrator(), fluid( size, sizeAllocated, h, numberOfCells ),
     boundary( size_bound, sizeAllocated_bound, h, numberOfCells ){};

   /**
    * Perform neighbors search and fill neighborsList in Particle system variable.
    */
   void PerformNeighborSearch( GlobalIndexType step, TNL::Timer& timer_reset, TNL::Timer& timer_cellIndices, TNL::Timer& timer_sort, TNL::Timer& timer_toCells );

   /**
    * Perform interaction for all particles, i.e. for all types.
    */
   template< typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm, typename EOS >
   void Interact();

   template< typename SPHKernelFunction, typename RiemannSolver, typename EOS >
   void Interact();

//protected:

   FluidPointer fluid;
   BoundaryPointer boundary;

   ModelPointer model;

   IntegratorPointer integrator;

};

} // SPH
} // ParticleSystem
} // TNL

#include "SPH_impl.h"

