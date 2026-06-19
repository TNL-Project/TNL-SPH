#include "SPH/configSetup.h"
#include "SolverMultiSet.h"
#include "SPH/OpenBoundaryConfig.h"
#include "SPH/TimeMeasurement.h"
#include <TNL/Config/ConfigDelimiter.h>
#include <TNL/Config/ConfigDescription.h>
#include <TNL/Config/ParameterContainer.h>
#include <TNL/Logger.h>
#include <iterator>
#include <ostream>
#include <string>

#include "distributedUtils.h"
#include <TNL/Particles/Writers/writeBackgroundGrid.h>

#include <SPH/configInit.h>

namespace TNL {
namespace SPH {

template< typename Model >
void
SolverMultiSet< Model >::init( int argc, char* argv[] )
{
   auto& params = this->parameters;
   auto& log = this->logger;

   try {
      initialize< SimulationType >( argc, argv, this->cliParams, this->cliConfig, params, this->config );
   }
   catch ( ... ) {
      std::cerr << std::endl;
   }

   log.writeHeader( "SPH simulation initialization." );
#ifdef HAVE_MPI
   log.writeParameter( "Configuration of distributed simulation:", "" );
   Containers::StaticVector< 2, int > numberOfSubdomains = params.getXyz< Containers::StaticVector< 2, int > >( "subdomains" );
   log.writeParameter( "Number of subdomains:", numberOfSubdomains );
   for( int x = 0; x < numberOfSubdomains[ 0 ]; x++ )
      for( int y = 0; y < numberOfSubdomains[ 1 ]; y++ )
         TNL::SPH::configSetupDistributedSubdomain( x, y, this->configDistributed );
   std::string configDistributedPath = params.template getParameter< std::string >( "distributed-config" );
   parseDistributedConfig( configDistributedPath, this->parametersDistributed, this->configDistributed, log );

   initDistributedParticleSets( params, this->parametersDistributed, log );

   this->loadBalancingMeasure = params.template getParameter< std::string >( "load-balancing-measure" );
   this->loadBalancingStepInterval = params.template getParameter< int >( "load-balancing-step-inteval" );
   this->fluid()->getDistributedParticles()->setParticlesCountResizeTrashold(
         params.template getParameter< float >( "number-of-particles-balancing-coef" ) );
   this->fluid()->getDistributedParticles()->setCompTimeResizePercetnageTrashold(
         params.template getParameter< float >( "computational-time-balancing-coef" ) );

   this->timeMeasurement.addTimer( "synchronize", false );
   this->timeMeasurement.addTimer( "rebalance", false );
#else
   initParticleSets( params, log );
#endif

   if( params.template getParameter< std::string >( "open-boundary-config" ) != "" ){
      this->initOpenBoundaryPatches( params, log );

      this->timeMeasurement.addTimer( "extrapolate-openbc" );
      this->timeMeasurement.addTimer( "apply-openbc" );
   }

   if( params.template getParameter< std::string >( "periodic-boundary-config" ) != "" ){
      this->initPeriodicBoundaryPatches( params, log );

      this->timeMeasurement.addTimer( "enforce-periodic-bc" );
      this->timeMeasurement.addTimer( "transfer-periodic-bc" );
      this->timeMeasurement.addTimer( "periodicity-fluid-updateZone", false );
      this->timeMeasurement.addTimer( "periodicity-boundary-updateZone", false );
   }

   this->modelParams.init( params );

   this->timeStepping.setTimeStep( params.template getParameter< RealType >( "initial-time-step" ) );
   this->timeStepping.setEndTime( params.template getParameter< RealType >( "final-time" ) );
   this->timeStepping.addOutputTimer( "save_results", params.template getParameter< RealType >( "snapshot-period" ) );

   this->caseName = params.template getParameter< std::string >( "case-name" );
   this->verbose = params.template getParameter< std::string >( "verbose-intensity" );
   this->outputDirectory = params.template getParameter< std::string >( "output-directory" );
   this->particlesFormat = params.template getParameter< std::string >( "particles-format" );

#ifdef HAVE_MPI
   readParticleFilesDistributed( params, this->parametersDistributed, log );
#else
   readParticlesFiles( params, log );
#endif

   log.writeSeparator();
   if( params.template getParameter< std::string >( "measuretool-config" ) != "" ) {
      log.writeParameter( "Simulation monitor initialization.", "" );
      this->simulationMonitor.init( params, this->timeStepping, log );
      log.writeParameter( "Simulation monitor initialization.", "Done." );
   }

   log.writeHeader( "SPH simulation successfully initialized." );
}

template< typename Model >
void
SolverMultiSet< Model >::initParticleSets( TNL::Config::ParameterContainer& parameters, TNL::Logger& logger )
{
   this->numberOfSubsets = 1;
   this->fluidSets.resize( 1 );
   this->boundarySets.resize( 1 );

   const VectorType domainOrigin = parameters.getXyz< VectorType >( "domainOrigin" );
   const VectorType domainSize = parameters.getXyz< VectorType >( "domainSize" );
   const RealType searchRadius = parameters.getParameter< RealType >( "searchRadius" );
   const IndexVectorType gridSize = TNL::ceil( ( domainSize - domainOrigin ) / searchRadius );

   this->fluid()->initialize( parameters.getParameter< int >( "numberOfParticles" ),
                               parameters.getParameter< int >( "numberOfAllocatedParticles" ),
                               searchRadius,
                               gridSize,
                               domainOrigin );
   if constexpr( ParticlesType::specifySearchedSetExplicitly() == true )
      this->fluid()->getParticles()->setParticleSetLabel( 0 );

   this->boundary()->initialize( parameters.getParameter< int >( "numberOfBoundaryParticles" ),
                                  parameters.getParameter< int >( "numberOfAllocatedBoundaryParticles" ),
                                  searchRadius,
                                  gridSize,
                                  domainOrigin );
   if constexpr( ParticlesType::specifySearchedSetExplicitly() == true )
      this->boundary()->getParticles()->setParticleSetLabel( 1 );
}

#ifdef HAVE_MPI
template< typename Model >
void
SolverMultiSet< Model >::initDistributedParticleSets( TNL::Config::ParameterContainer& parameters,
                                                       TNL::Config::ParameterContainer& parametersDistributed,
                                                       TNL::Logger& logger )
{
   this->numberOfSubsets = 1;
   this->fluidSets.resize( 1 );
   this->boundarySets.resize( 1 );

   logger.writeHeader( "SPH simulation initialization." );
   int rank = TNL::MPI::GetRank();
   Containers::StaticVector< 2, int > numberOfSubdomains = parameters.getXyz< Containers::StaticVector< 2, int > >( "subdomains" );
   Containers::StaticVector< 2, int > subdomainCoordinates = distributed::restoreSubdomainCoordinatesFromRank( rank, numberOfSubdomains );
   const std::string subdomainKey = distributed::getSubdomainKey( rank, numberOfSubdomains );
   logger.writeParameter( "Initializing rank: ", rank );
   logger.writeParameter( "Initializing subdomain: ", subdomainCoordinates );
   logger.writeParameter( "Subdomain key: ", subdomainKey );
   logger.writeSeparator();

   const RealType searchRadius = parameters.getParameter< RealType >( "searchRadius" );
   const VectorType domainOrigin = parameters.getXyz< VectorType >( "domainOrigin" );
   const VectorType domainSize = parameters.getXyz< VectorType >( "domainSize" );
   const IndexVectorType domainGridDimension = TNL::ceil( ( domainSize - domainOrigin ) / searchRadius );
   const int numberOfOverlapLayers = parameters.getParameter< int >( "overlapWidth" );

   const VectorType subdomainOrigin = parametersDistributed.getXyz< VectorType >( subdomainKey + "origin" );
   const IndexVectorType subdomainGridDimension = parametersDistributed.getXyz< IndexVectorType >( subdomainKey + "grid-dimensions" );
   const IndexVectorType subdomainGridOriginGlobalCoords = parametersDistributed.getXyz< IndexVectorType >( subdomainKey + "origin-global-coords" );

   logger.writeParameter( "initDistributed:", "fluid->initialize" );
   this->fluid()->initializeAsDistributed( parametersDistributed.getParameter< int >( subdomainKey + "fluid_n" ),
                                            parametersDistributed.getParameter< int >( subdomainKey + "fluid_n_allocated" ),
                                            searchRadius,
                                            domainGridDimension,
                                            domainOrigin,
                                            subdomainGridDimension,
                                            subdomainGridOriginGlobalCoords,
                                            numberOfOverlapLayers,
                                            numberOfSubdomains,
                                            subdomainOrigin,
                                            logger );
   this->fluid()->getDistributedParticles()->writeProlog( logger );
   this->fluid()->setCommunicator( this->communicator );
   this->fluid()->getDistributedParticlesSynchronizer().initialize( this->fluid()->getDistributedParticles() );

   logger.writeParameter( "initDistributed:", "boundary->initialize" );
   this->boundary()->initializeAsDistributed( parametersDistributed.getParameter< int >( subdomainKey + "boundary_n" ),
                                               parametersDistributed.getParameter< int >( subdomainKey + "boundary_n_allocated" ),
                                               searchRadius,
                                               domainGridDimension,
                                               domainOrigin,
                                               subdomainGridDimension,
                                               subdomainGridOriginGlobalCoords,
                                               numberOfOverlapLayers,
                                               numberOfSubdomains,
                                               subdomainOrigin,
                                               logger );
   this->boundary()->getDistributedParticles()->writeProlog( logger );
   this->boundary()->setCommunicator( this->communicator );
   this->boundary()->getDistributedParticlesSynchronizer().initialize( this->boundary()->getDistributedParticles() );
}
#endif

template< typename Model >
void
SolverMultiSet< Model >::readParticlesFiles( TNL::Config::ParameterContainer& parameters, TNL::Logger& logger )
{
   if( parameters.getParameter< int >( "numberOfParticles" ) != 0 ){
      logger.writeParameter( "Reading fluid particles:", parameters.getParameter< std::string >( "fluid-particles" ) );
      this->fluid()->template readParticlesAndVariables< typename BaseType::SimulationReaderType >(
         parameters.getParameter< std::string >( "fluid-particles" ) );
   }
   else{
      logger.writeParameter( "Reading fluid particles:", "NO PARTICLES TO READ" );
   }

   if( parameters.getParameter< int >( "numberOfBoundaryParticles" ) != 0 ){
      logger.writeParameter( "Reading boundary particles:", parameters.getParameter< std::string >( "boundary-particles" ) );
      this->boundary()->template readParticlesAndVariables< typename BaseType::SimulationReaderType >(
         parameters.getParameter< std::string >( "boundary-particles" ) );
   }
   else{
      logger.writeParameter( "Reading boundary particles:", "NO PARTICLES TO READ" );
   }

   const int numberOfBoundaryPatches = parameters.getParameter< int >( "openBoundaryPatches" );
   if( numberOfBoundaryPatches > 0 ) {
      for( int i = 0; i < numberOfBoundaryPatches; i++ ) {
         std::string prefix = "buffer-" + std::to_string( i + 1 ) + "-";
         if( this->parametersOpenBoundary.getParameter< int >( prefix + "numberOfParticles" ) != 0 ){
            logger.writeParameter( "Reading open boundary particles:", this->parametersOpenBoundary.getParameter< std::string >( prefix + "particles" ) );
            this->openBoundaryPatches[ i ]->template readParticlesAndVariables< typename BaseType::SimulationReaderType >(
               this->parametersOpenBoundary.getParameter< std::string >( prefix + "particles" ) );
         }
         else{
            logger.writeParameter( "Reading open boundary particles:", "NO PARTICLES TO READ" );
         }
      }
   }
}

#ifdef HAVE_MPI
template< typename Model >
void
SolverMultiSet< Model >::readParticleFilesDistributed( TNL::Config::ParameterContainer& parameters,
                                                        TNL::Config::ParameterContainer& parametersDistributed,
                                                        TNL::Logger& logger )
{
   int rank = TNL::MPI::GetRank();
   Containers::StaticVector< 2, int > numberOfSubdomains = parameters.getXyz< Containers::StaticVector< 2, int > >( "subdomains" );
   const std::string subdomainKey = distributed::getSubdomainKey( rank, numberOfSubdomains );

   logger.writeParameter( "Reading fluid particles:", parametersDistributed.getParameter< std::string >( subdomainKey + "fluid-particles" ) );
   this->fluid()->template readParticlesAndVariables< typename BaseType::SimulationReaderType >(
      parametersDistributed.getParameter< std::string >( subdomainKey + "fluid-particles" ) );
   logger.writeParameter( "Reading boundary particles:", parametersDistributed.getParameter< std::string >( subdomainKey + "boundary-particles" ) );
   this->boundary()->template readParticlesAndVariables< typename BaseType::SimulationReaderType >(
      parametersDistributed.getParameter< std::string >( subdomainKey + "boundary-particles" ) );

   const int numberOfBoundaryPatches = parameters.getParameter< int >( "openBoundaryPatches" );
   if( this->openBoundaryPatches.size() > 0 ) {
      for( int i = 0; i < numberOfBoundaryPatches; i++ ) {
         std::string prefix = subdomainKey + "buffer-" + std::to_string( i + 1 ) + "-";
         this->openBoundaryPatches[ i ]->template readParticlesAndVariables< typename BaseType::SimulationReaderType >(
            parametersDistributed.getParameter< std::string >( prefix + "particles" ) );
      }
   }
}
#endif

template< typename Model >
void
SolverMultiSet< Model >::interact()
{
   this->timeMeasurement.start( "interact" );

   this->model.updateSolidBoundary( this->fluid(), this->boundary(), this->modelParams );
   if( this->openBoundaryPatches.size() > 0 ) {
      for( long unsigned int i = 0; i < this->openBoundaryPatches.size(); i++ ) {
         this->model.updateSolidBoundaryOpenBoundary( this->boundary(), this->openBoundaryPatches[ i ], this->modelParams );
      }
   }
   this->model.finalizeBoundaryInteraction( this->fluid(), this->boundary(), this->modelParams );

   this->model.interaction( this->fluid(), this->boundary(), this->modelParams );
   if( this->openBoundaryPatches.size() > 0 ) {
      for( long unsigned int i = 0; i < this->openBoundaryPatches.size(); i++ ) {
         this->openBoundaryPatches[ i ]->zone.updateParticlesInZone( this->fluid()->getParticles() );
         this->model.interactionWithOpenBoundary( this->fluid(), this->openBoundaryPatches[ i ], this->modelParams );
      }
   }
   this->model.finalizeInteraction( this->fluid(), this->boundary(), this->modelParams );

   this->timeMeasurement.stop( "interact" );
   this->writeLog( "Interact...", "Done." );
}

} // namespace SPH
} // namespace TNL
