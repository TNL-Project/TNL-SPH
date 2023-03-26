#pragma once

#include <TNL/Logger.h>
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
   GlobalIndexType gridNumberOfCells;
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
         this->gridNumberOfCells = Config::gridXsize * Config::gridYsize;
      }
      else if constexpr( ParticleConfig::spaceDimension == 3 )
      {
         this->gridSize = { Config::gridXsize, Config::gridYsize, Config::gridZsize };
         this->gridOrigin = { Config::gridXorigin, Config::gridYorign, Config::gridZsize };
         this->gridNumberOfCells = Config::gridXsize * Config::gridYsize * Config::gridZsize;
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

   SPHSimpleFluid( SPHSimpleFluidConfig< typename ParticleSystem::Config > sphConfig )
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
   void
   PerformNeighborSearch( GlobalIndexType step, TNL::Timer& timer_reset, TNL::Timer& timer_cellIndices, TNL::Timer& timer_sort, TNL::Timer& timer_toCells );

   /**
    * Perform interaction for all particles, i.e. for all types.
    */
   template< typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm, typename EOS >
   void
   Interact();

   template< typename SPHKernelFunction, typename RiemannSolver, typename EOS >
   void
   Interact();

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

template< typename Model, typename ParticleSystem, typename NeighborSearch >
std::ostream&
operator<<( std::ostream& str, const SPH::SPHSimpleFluid< Model, ParticleSystem, NeighborSearch >& sphSimulation )
{
   TNL::Logger logger( 100, str );

   sphSimulation.writeProlog( logger );

   return str;
}


#include "SPH_impl.h"

