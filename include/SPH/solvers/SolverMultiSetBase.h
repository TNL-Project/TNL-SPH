#pragma once

#include <memory>
#include <filesystem>
#include <iterator>
#include <ostream>
#include <string>

#include <TNL/Containers/Vector.h>
#include <TNL/Algorithms/reduce.h>
#include <TNL/Config/ConfigDelimiter.h>
#include <TNL/Config/ConfigDescription.h>
#include <TNL/Config/ParameterContainer.h>
#include <TNL/Logger.h>

#include "../Fluid.h"
#include "../Boundary.h"
#include "../OpenBoundaryBuffers.h"
#include "../OpenBoundaryConfig.h"
#include "../TimeMeasurement.h"
#include "../SimulationMonitor.h"

namespace TNL {
namespace SPH {

template< typename Model >
class SolverMultiSetBase
{
public:

   using DeviceType = typename Model::DeviceType;
   using ModelType = Model;
   using ModelParams = typename ModelType::ModelParams;
   using ParticlesType = typename ModelType::ParticlesType;
   using IntegrationSchemeType = typename ModelType::IntegrationSchemeType;
   using IntegratorPointer = typename Pointers::SharedPointer< IntegrationSchemeType, DeviceType >;
   using IntegrationSchemeVariablesType = typename Model::IntegrationSchemeVariables;
   using TimeStepping = typename Model::ModelParams::TimeStepping;

   using SPHConfig = typename Model::SPHConfig;
   using GlobalIndexType = typename ParticlesType::GlobalIndexType;
   using RealType = typename ParticlesType::RealType;
   using IndexVectorType = typename ParticlesType::IndexVectorType;
   using VectorType = typename ParticlesType::PointType;

   using FluidVariables = typename Model::FluidVariables;
   using Fluid = Fluid< ParticlesType, SPHConfig, FluidVariables, IntegrationSchemeVariablesType >;
   using FluidPointer = Pointers::SharedPointer< Fluid, DeviceType >;
   using BoundaryVariables = typename Model::BoundaryVariables;
   using Boundary = Boundary< ParticlesType, SPHConfig, BoundaryVariables, IntegrationSchemeVariablesType >;
   using BoundaryPointer = Pointers::SharedPointer< Boundary, DeviceType >;
   using OpenBoundaryVariables = typename Model::OpenBoundaryVariables;
   using OpenBoundaryConfigType = typename Model::OpenBoundaryConfig;
   using OpenBoundary = OpenBoundary<
      ParticlesType, SPHConfig, OpenBoundaryVariables, IntegrationSchemeVariablesType, OpenBoundaryConfigType >;
   using OpenBoundaryPointer = Pointers::SharedPointer< OpenBoundary, DeviceType >;
   using OpenBoundaryModel = typename Model::OpenBoundaryModel;

   using Reader = TNL::ParticleSystem::Readers::VTKReader;
   using Writer = TNL::ParticleSystem::Writers::VTKWriter< ParticlesType >;
   using SimulationReaderType = TNL::ParticleSystem::ReadParticles< typename ParticlesType::Config, Reader >;
   using ComputationTimeMeasurement = TNL::SPH::TimerMeasurement;
   using SimulationMonitor = SimulationMonitor< SolverMultiSetBase< Model > >;

   SolverMultiSetBase( std::ostream& out = std::cout ) : logger( 100, out ) {};

   void
   initOpenBoundaryPatches( TNL::Config::ParameterContainer& parameters, TNL::Logger& logger );

   void
   initPeriodicBoundaryPatches( TNL::Config::ParameterContainer& parameters, TNL::Logger& logger );

   void
   performNeighborSearch( bool performBoundarySearch = false );

   template< typename ParticleSetPointer >
   void
   performNeighborSearchForObject( ParticleSetPointer& objectPointer );

   void
   removeParticlesOutOfDomain();

   void
   removeParticlesOutOfDensityLimits();

   void
   extrapolateOpenBC();

   void
   applyOpenBC( const RealType timeStepFact = 1.f );

   void
   applyPeriodicBCEnforce();

   void
   applyPeriodicBCTransfer();

   void
   computeTimeStep();

   void
   updateTime();

   void
   measure();

   template< typename Stage = int >
   void
   integrate( const Stage integrationStage = Stage{}, const bool integrateBoundary = false );

   void
   integrateVerletStep( const bool integrateBoundary = false );

   void
   symplecticVerletPredictor();

   void
   symplecticVerletCorrector();

   void
   midpointPredictor();

   void
   midpointUpdateVariables();

   void
   midpointResidualsAndRelaxationFactor();

   void
   midpointRelax();

   void
   midpointCorrector();

   virtual void
   save( bool writeParticleCellIndex = false );

   void
   makeSnapshot();

   template< typename Func >
   void
   initUserConfig( Func&& userConfigFunction );

   void
   writeProlog( bool writeSystemInformation = true ) noexcept;

   template< typename ParameterType >
   void
   writeLog( const std::string& label, const ParameterType& value, int parameterLevel = 0 );

   void
   writeInfo() noexcept;

   void
   writeEpilog() noexcept;

#ifdef HAVE_MPI

   void
   synchronizeDistributedSimulation();

   void
   resetOverlaps();

   void
   performLoadBalancing();

   void
   writeLoadBalancingInfo( const int gridResize );

   void
   balanceSubdomains();

#endif

   FluidPointer&
   fluid()
   {
      return fluidSets[ 0 ];
   }

   BoundaryPointer&
   boundary()
   {
      return boundarySets[ 0 ];
   }

   int numberOfSubsets = 0;
   std::vector< FluidPointer > fluidSets;
   std::vector< BoundaryPointer > boundarySets;
   std::vector< OpenBoundaryPointer > openBoundaryPatches;

   Model model;
   ModelParams modelParams;
   OpenBoundaryModel openBoundaryModel;

   IntegratorPointer integrator;

   TimeStepping timeStepping;
   ComputationTimeMeasurement timeMeasurement;

   std::string caseName;
   std::string verbose = "none";
   std::string outputDirectory;
   std::string particlesFormat;
   SimulationMonitor simulationMonitor;

   TNL::Config::ParameterContainer cliParams;
   TNL::Config::ConfigDescription cliConfig;

   TNL::Config::ParameterContainer parameters;
   TNL::Config::ConfigDescription config;

   TNL::Logger logger;

   TNL::Config::ConfigDescription configOpenBoundary;
   TNL::Config::ParameterContainer parametersOpenBoundary;

   TNL::Config::ParameterContainer userParams;
   TNL::Config::ConfigDescription userConfig;

   GlobalIndexType totalNumberOfParticlesOutOfDomain = 0;
   GlobalIndexType totalNumberOfParticlesOutOfDensityLimits = 0;

#ifdef HAVE_MPI
   MPI::Comm communicator = MPI_COMM_WORLD;
   TNL::Config::ConfigDescription configDistributed;
   TNL::Config::ParameterContainer parametersDistributed;
   RealType subdomainCompTimeBackup = 0;
   std::string loadBalancingMeasure;
   int loadBalancingStepInterval = 1;
#endif

};

} // namespace SPH
} // namespace TNL

#include "SolverMultiSetBase.hpp"

