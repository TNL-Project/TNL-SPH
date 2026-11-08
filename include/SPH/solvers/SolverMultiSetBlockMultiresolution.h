// SPDX-FileComment: This file is part of TNL-SPH (https://gitlab.com/tnl-project/tnl-sph)
// SPDX-License-Identifier: MIT

#pragma once

#include <TNL/Meshes/Grid.h>

#include <SPH/solvers/SolverMultiSetBase.h>
#include <SPH/BoundaryGhostParticles.h>
#include <SPH/DecompositionTopology.h>
#include <SPH/MultiresolutionRectangleBuffer.h>
#include <SPH/VariablesGhostInit.h>

namespace TNL::SPH {

template< typename Model >
class SolverMultiSetBlockMultiresolution : public SolverMultiSetBase< Model >
{
public:
   using BaseType = SolverMultiSetBase< Model >;
   using SimulationType = SolverMultiSetBlockMultiresolution< Model >;

   using typename BaseType::CoordinatesType;
   using typename BaseType::DeviceType;
   using typename BaseType::GlobalIndexType;
   using typename BaseType::ModelParams;
   using typename BaseType::ModelType;
   using typename BaseType::ParticlesType;
   using typename BaseType::RealType;
   using typename BaseType::SPHConfig;
   using typename BaseType::VectorType;

   using SolverTraitsType = SPHFluidTraits< SPHConfig >;
   using IndexArrayType = typename SolverTraitsType::IndexArrayType;

   using typename BaseType::BoundaryPointer;
   using typename BaseType::FluidPointer;
   using typename BaseType::FluidVariables;
   using typename BaseType::OpenBoundaryConfigType;
   using typename BaseType::OpenBoundaryPointer;
   using typename BaseType::SimulationMonitor;

   using MultiresolutionBoundary = MultiresolutionBoundary<
      ParticlesType,
      SPHConfig,
      FluidVariables,
      typename BaseType::IntegrationSchemeVariablesType,
      OpenBoundaryConfigType,
      ModelParams >;
   using MultiresolutionBoundaryPointer = Pointers::SharedPointer< MultiresolutionBoundary, DeviceType >;

   using GridType = TNL::Meshes::Grid< SPHConfig::spaceDimension, RealType >;
   using Topology = DecompositionTopologyFlat< GridType >;

   SolverMultiSetBlockMultiresolution( std::ostream& out = std::cout )
   : BaseType( out )
   {}

   void
   init( int argc, char* argv[] );

   void
   initializeBlockBasedMultiResolutionSimulation();

   void
   initParticleSets();

   void
   initMultiResolutionBoundaryPatches();

   void
   readParticlesFiles();

   void
   multiresolutionUpdate();

   void
   interact();

   void
   save( bool writeParticleCellIndex = false );

   void
   writeProlog( bool writeSystemInformation = true, bool writeMultiresolutionBuffers = true ) noexcept;

   void
   writeMultiresolutionBoundaryPatchesProlog( const std::string& title ) noexcept;

   virtual void
   writeInfo( bool writeMultiresolutionBuffers = true ) noexcept;

   void
   writeMultiresolutionBoundaryPatchesInfo( const std::string& label ) noexcept;

#ifdef HAVE_MPI
   void
   initializeDistributedSimulation();

   void
   initDistributedParticleSets(
      TNL::Config::ParameterContainer& parameters,
      TNL::Config::ParameterContainer& parametersDistributed,
      TNL::Logger& logger );

   void
   readParticleFilesDistributed(
      TNL::Config::ParameterContainer& parameters,
      TNL::Config::ParameterContainer& parametersDistributed,
      Logger& logger );
#endif

   /* Boundary ghosts - see Documentation/multiresolution-ghost-boundaries-design.md */
   void
   initBoundaryGhosts();

   void
   refreshBoundaryGhostValues();

   void
   refreshBoundaryGhostValues( int ownIdx );

   Topology topology;
   std::vector< MultiresolutionBoundaryPointer > multiresolutionBoundaryPatches;
   std::vector< typename Topology::Interface > multiresolutionBoundaryPatchInterfaces;

   std::vector< BoundaryGhostParticles< SPHConfig > > boundaryGhostInterfaces;
   BoundaryGhostUpdate boundaryGhostUpdate = BoundaryGhostUpdate::Interpolation;

   TNL::Config::ConfigDescription configSubdomains;
   TNL::Config::ParameterContainer parametersSubdomains;
};

}  // namespace TNL::SPH

template< typename Model >
std::ostream&
operator<<( std::ostream& str, const TNL::SPH::SolverMultiSetBlockMultiresolution< Model >& sphSimulation )
{
   TNL::Logger logger( 100, str );

   sphSimulation.writeProlog( logger );

   return str;
}

#include <SPH/solvers/SolverMultiSetBlockMultiresolution.hpp>
