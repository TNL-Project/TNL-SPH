#pragma once

#include <TNL/Meshes/Grid.h>

#include "SolverMultiSetBase.h"
#include "../DecompositionTopology.h"
#include "../MultiresolutionRectangleBuffer.h"

namespace TNL {
namespace SPH {

template< typename Model >
class SolverMultiSetBlockMultiresolution : public SolverMultiSetBase< Model >
{
public:

   using BaseType = SolverMultiSetBase< Model >;
   using SimulationType = SolverMultiSetBlockMultiresolution< Model >;

   using typename BaseType::DeviceType;
   using typename BaseType::ModelType;
   using typename BaseType::ModelParams;
   using typename BaseType::ParticlesType;
   using typename BaseType::SPHConfig;
   using typename BaseType::GlobalIndexType;
   using typename BaseType::RealType;
   using typename BaseType::IndexVectorType;
   using typename BaseType::VectorType;

   using typename BaseType::FluidVariables;
   using typename BaseType::FluidPointer;
   using typename BaseType::BoundaryPointer;
   using typename BaseType::OpenBoundaryPointer;
   using typename BaseType::OpenBoundaryConfigType;
   using typename BaseType::SimulationMonitor;

   using MultiresolutionBoundary = MultiresolutionBoundary<
      ParticlesType, SPHConfig, FluidVariables, typename BaseType::IntegrationSchemeVariablesType, OpenBoundaryConfigType, ModelParams >;
   using MultiresolutionBoundaryPointer = Pointers::SharedPointer< MultiresolutionBoundary, DeviceType >;

   using GridType = TNL::Meshes::Grid< SPHConfig::spaceDimension, RealType >;
   using Topology = DecompositionTopologyFlat< GridType >;

   SolverMultiSetBlockMultiresolution( std::ostream& out = std::cout ) : BaseType( out ) {};

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
   writeProlog( bool writeSystemInformation = true ) noexcept;

   void
   writeInfo() noexcept;

#ifdef HAVE_MPI
   void
   initializeDistributedSimulation();

   void
   initDistributedParticleSets( TNL::Config::ParameterContainer& parameters,
                                TNL::Config::ParameterContainer& parametersDistributed,
                                TNL::Logger& logger );

   void
   readParticleFilesDistributed( TNL::Config::ParameterContainer& parameters,
                                 TNL::Config::ParameterContainer& parametersDistributed,
                                 TNL::Logger& logger );
#endif

   Topology topology;
   std::vector< MultiresolutionBoundaryPointer > multiresolutionBoundaryPatches;

   TNL::Config::ConfigDescription configSubdomains;
   TNL::Config::ParameterContainer parametersSubdomains;

};

} // namespace SPH
} // namespace TNL

template< typename Model >
std::ostream&
operator<<( std::ostream& str, const TNL::SPH::SolverMultiSetBlockMultiresolution< Model >& sphSimulation )
{
   TNL::Logger logger( 100, str );

   sphSimulation.writeProlog( logger );

   return str;
}

#include "SolverMultiSetBlockMultiresolution.hpp"

