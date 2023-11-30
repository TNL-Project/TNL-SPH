#pragma once

#include <TNL/Containers/Vector.h>
#include <TNL/Algorithms/reduce.h>
#include <memory> //shared_ptr

#include "../Particles/ParticlesTraits.h"
#include "SPH.h"

#include "Fluid.h"
#include "Boundary.h"
#include "OpenBoundaryBuffers.h"
#include "OpenBoundaryConfig.h"
#include "TNL/Functional.h"
#include "TNL/Logger.h"
#include "TimeMeasurement.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Model >
class SPHOpenSystem
{
public:

   using DeviceType = typename Model::DeviceType;
   using ParticleSystemType = typename Model::ParticlesType;; //Added due to measure tool
   using ParticleSystem = typename Model::ParticlesType;; //Added due to measure tool

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

   using OpenBoundaryVariables = typename Model::OpenBoundaryVariables;
   using OpenBoundary = OpenBoundary< ParticleSystem, SPHConfig, OpenBoundaryVariables, IntegratorVariables >;
   using OpenBoundaryPointer = Pointers::SharedPointer< OpenBoundary, DeviceType >;

   //using OpenBoundaryConfigType = OpenBoundaryConfig< TNL::ParticleSystem::SPH::None, SPHConfig >;
   using OpenBoundaryConfigType = typename Model::OpenBoundaryConfig;

   SPHOpenSystem() = default;

   template< typename SPHOpenSystemInit >
   SPHOpenSystem( SPHOpenSystemInit sphConfig )
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
    * TODO: I don't like this.
    */
   template< typename SPHOpenSystemInit >
   void
   addOpenBoundaryPatch( SPHOpenSystemInit sphConfig );

   template< typename SPHOpenSystemInit >
   void
   addOpenBoundaryPatch( SPHOpenSystemInit sphConfig, std::vector< OpenBoundaryConfigType >& bufferParams );

   /**
    * Perform neighbors search and fill neighborsList in Particle system variable.
    */
   void
   performNeighborSearch( GlobalIndexType step, TimerMeasurement& timeMeasurement, TNL::Logger& log);

   template< typename PhysicalObjectPointer >
   void
   performNeighborSearchForObject( const GlobalIndexType& step,
                                   PhysicalObjectPointer& objectPointer,
                                   TimerMeasurement& timeMeasurement,
                                   TNL::Logger& log);

   template< typename PhysicalObjectPointer >
   void
   performNeighborSearchForOpenBoundaryPatches( const GlobalIndexType& step,
                                                TimerMeasurement& timeMeasurement,
                                                TNL::Logger& log );
   /**
    * \brief Perform interaction between all particles and all particle objects
    * in the simulation.
    */
   template< typename SPHKernelFunction, typename EOS, typename SPHState >
   void
   extrapolateOpenBC( SPHState& sphState, std::vector< OpenBoundaryConfigType >& bufferParams );

   /**
    * \brief Perform interaction between all particles and all particle objects
    * in the simulation.
    */
   template< typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm, typename EOS, typename SPHState >
   void
   interact( SPHState& sphState );

   /**
    * \brief Save all particle object to vtk files. Automatically saves all
    * available fileds.
    */
   template< typename Writer >
   void
   save( const std::string& outputFilename, const int step, bool writeParticleCellIndex = false  );

   void
   writeProlog( TNL::Logger& logger ) const noexcept;

//protected:

   FluidPointer fluid;
   BoundaryPointer boundary;
   std::vector< OpenBoundaryPointer > openBoundaryPatches;
   std::vector< OpenBoundaryConfigType > openBoundaryPatchesConfigs;

   //std::map< std::string, int > openBoundaryPatchesMap;

   ModelPointer model;

   IntegratorPointer integrator;

};

} // SPH
} // ParticleSystem
} // TNL

template< typename Model >
std::ostream&
operator<<( std::ostream& str, const SPH::SPHOpenSystem< Model >& sphSimulation )
{
   TNL::Logger logger( 100, str );

   sphSimulation.writeProlog( logger );

   return str;
}

#include "SPHOpen_impl.h"

