#pragma once

#include <TNL/Containers/Vector.h>
#include <TNL/Algorithms/reduce.h>
#include <memory> //shared_ptr

#include "Fluid.h"
#include "Boundary.h"
#include "OpenBoundaryBuffers.h"
#include "OpenBoundaryConfig.h"
#include "TNL/Functional.h"
#include "TNL/Logger.h"
#include "TimeMeasurement.h"

#include "SimulationMonitor.h"

// custom shared modules
#include "shared/removeParticlesOutOfDensityLimits.h"

namespace TNL {
namespace SPH {

template< typename Model >
class SPHMultiset_CFD
{
public:

   using SimulationType = SPHMultiset_CFD< Model >;
   using DeviceType = typename Model::DeviceType;
   using ModelType = Model;
   using ModelParams = typename ModelType::ModelParams;
   using ParticlesType = typename ModelType::ParticlesType;;
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

   //Reader
   using Reader = TNL::ParticleSystem::Readers::VTKReader;
   using Writer = TNL::ParticleSystem::Writers::VTKWriter< ParticlesType >;
   using SimulationReaderType = TNL::ParticleSystem::ReadParticles< typename ParticlesType::Config, Reader >;
   using ComputationTimeMeasurement = TNL::SPH::TimerMeasurement;
   using SimulationMonitor = SimulationMonitor< SimulationType >;

   //SPHMultiset_CFD() = default;

   //SPHMultiset_CFD() : logger( 100, std::cout ) {};

   SPHMultiset_CFD( std::ostream& out = std::cout ) : logger( 100, out ) {};

   void
   init( int argc, char* argv[] );

   // protected
   void
   initParticleSets( TNL::Config::ParameterContainer& parameters, TNL::Logger& logger );

   // protected
   void
   initOpenBoundaryPatches( TNL::Config::ParameterContainer& parameters, TNL::Logger& logger );

   void
   initPeriodicBoundaryPatches( TNL::Config::ParameterContainer& parameters, TNL::Logger& logger );

   // protected
   void
   readParticlesFiles( TNL::Config::ParameterContainer& parameters, TNL::Logger& logger );

   // protected
   void
   readOpenBoundaryFiles( TNL::Config::ParameterContainer& parameters, TNL::Logger& logger );


#ifdef HAVE_MPI
   // protected
   void
   initDistributedParticleSets( TNL::Config::ParameterContainer& parameters,
                                TNL::Config::ParameterContainer& parametersDistributed,
                                TNL::Logger& logger );

   // protected
   void
   readParticleFilesDistributed( TNL::Config::ParameterContainer& parameters,
                                 TNL::Config::ParameterContainer& parametersDistributed,
                                 TNL::Logger& logger );
#endif

   // protected
   void
   initOverlaps( TNL::Config::ParameterContainer& parameters,
                 TNL::Config::ParameterContainer& parametersDistributed,
                 TNL::Logger& logger );

   /**
    * Perform neighbors search and fill neighborsList in Particle system variable.
    */
   void
   performNeighborSearch( bool performBoundarySearch = false );

   /**
    *
    */
   void
   removeParticlesOutOfDomain();

   /**
    * Wrapper for module removeing particles out of density limits.
    */
   void
   removeParticlesOutOfDensityLimits();
   //requires ( requires { &Model::action; std::is_member_function_pointer_v< decltype( &Model::action ) >; } )

   //void
   //performNeighborSearch( TNL::Logger& log, bool performBoundarySearch = false );

   //TODO: Should we have log in this functions?
   template< typename ParticleSetPointer >
   void
   performNeighborSearchForObject( ParticleSetPointer& objectPointer );

   template< typename ParticleSetPointer >
   void
   performNeighborSearchForOpenBoundaryPatches( TNL::Logger& log );

   /**
    * \brief Perform interaction between all particles and all particle objects
    * in the simulation.
    */
   void
   extrapolateOpenBC();

   /**
    * \brief Apply open boundary simulations i.e. processes which add and remove
    * particles into and from simulations (inflows and outlows).
    */
   void
   applyOpenBC( const RealType timeStepFact = 1.f );

   /**
    * \brief Apply periodic boundary conditions - first part. For all periodic
    * boundary zones, this function ensures that particles are copied from the
    * connected part of the domain where the periodic OP is prescribed. After
    * copying the data into ghost zones, these particles needs to be updated
    * in terms of neighbor search, so we can find them correctly through
    * neighbrosLoops.
    */
   void
   applyPeriodicBCEnforce();

   /**
    * \brief Apply periodic boundary conditions - second part. For all periodic
    * boundary zones, this function transfers the particles that entering the
    * periodic zone to connected part of the domain where the periodic OP is
    * prescribed.
    */
   void
   applyPeriodicBCTransfer();

   /**
    * \brief Perform interaction between all particles and all particle objects
    * in the simulation.
    */
   void
   interact();

   /**
    *
    */
   void
   computeTimeStep();

   /**
    *
    */
   void
   updateTime();

   /**
    * \brief Check if is time to perform measurement and if is time to perform
    * measurement, perform measurement.
    */
   void
   measure();

   /**
    * \brief General integrator wrapper
    */
   //template< typename IntegrationStage >
   //void
   //integrate( IntegrationStage stage );

   /**
    * \brief Perform integration setp - wrapper for Verlet integration scheme step.
    */
   void
   integrateVerletStep( const bool integrateBoundary = false );

   /**
    * \brief Perform integration setp - wrapper for Verlet integration scheme step.
    */
   void
   symplecticVerletPredictor();

   void
   symplecticVerletCorrector();

   /**
    * \brief Perform integration setp - wrapper for Midpoint scheme functions.
    */
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


#ifdef HAVE_MPI

   void
   synchronizeDistributedSimulation();

   void
   resetOverlaps();

   void
   performLoadBalancing();

   //TODO: Update to general dimensions
   void
   writeLoadBalancingInfo( const int gridResize );

   /**
    * \brief Wrapper around load balancing procedure.
    */
   void
   balanceSubdomains();


#endif

   /**
    * \brief Save all particle object to vtk files. Automatically saves all
    * available fileds.
    */
   void
   save( bool writeParticleCellIndex = false  );

   void
   makeSnapshot();

   void
   writeProlog( bool writeSystemInformation = true ) noexcept;

   template< typename ParameterType >
   void
   writeLog( const std::string& label, const ParameterType& value, int parameterLevel = 0 );

   void
   writeInfo() noexcept;

   void
   writeEpilog() noexcept;

//protected:

   FluidPointer fluid;
   BoundaryPointer boundary;
   std::vector< OpenBoundaryPointer > openBoundaryPatches;

   Model model;
   ModelParams modelParams;
   OpenBoundaryModel openBoundaryModel;

   IntegratorPointer integrator; // I hate this.

   TimeStepping timeStepping;
   ComputationTimeMeasurement timeMeasurement;

   // Track removed fluid particles
   GlobalIndexType totalNumberOfParticlesOutOfDomain = 0;
   GlobalIndexType totalNumberOfParticlesOutOfDensityLimits = 0;

   std::string caseName;
   std::string verbose = "none";
   std::string outputDirectory;
   std::string particlesFormat;
   SimulationMonitor simulationMonitor;

   // Initand control parameters
   TNL::Config::ParameterContainer cliParams;
   TNL::Config::ConfigDescription cliConfig;

   TNL::Config::ParameterContainer parameters;
   TNL::Config::ConfigDescription config;

   TNL::Logger logger;

   // TEMP: And btw the names are AWFUL
#ifdef HAVE_MPI
   MPI::Comm communicator = MPI_COMM_WORLD;
   TNL::Config::ConfigDescription configDistributed;
   TNL::Config::ParameterContainer parametersDistributed;
   RealType subdomainCompTimeBackup = 0;
   std::string loadBalancingMeasure;
   int loadBalancingStepInterval = 1;
#endif

   // Configurations and parameter configs (mostly required by initialization)
   TNL::Config::ConfigDescription configOpenBoundary;
   TNL::Config::ParameterContainer parametersOpenBoundary;

   TNL::Config::ConfigDescription configPeriodicBoundary;
   TNL::Config::ParameterContainer parametersPeriodicBoundary;

};

} // SPH
} // TNL

template< typename Model >
std::ostream&
operator<<( std::ostream& str, const TNL::SPH::SPHMultiset_CFD< Model >& sphSimulation )
{
   TNL::Logger logger( 100, str );

   sphSimulation.writeProlog( logger );

   return str;
}

#include "SPHMultiset_CFD.hpp"

