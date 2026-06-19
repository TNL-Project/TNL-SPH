#include "SolverMultiSetBlockMultiresolution.h"
#include "../tempFunctionsToConfigMultiresolution.h"

namespace TNL {
namespace SPH {

template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::init( int argc, char* argv[] )
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
   this->caseName = params.template getParameter< std::string >( "case-name" );
   this->verbose = params.template getParameter< std::string >( "verbose-intensity" );
   this->outputDirectory = params.template getParameter< std::string >( "output-directory" );
   this->particlesFormat = params.template getParameter< std::string >( "particles-format" );

#ifdef HAVE_MPI
   initializeDistributedSimulation();
#else
   initializeBlockBasedMultiResolutionSimulation();
#endif
   log.writeHeader( "SPH simulation successfully initialized." );
}

template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::initializeBlockBasedMultiResolutionSimulation()
{
   auto& params = this->parameters;
   auto& log = this->logger;

   this->numberOfSubsets = params.template getParameter< int >( "numberOfSubdomains" );
   const std::string configSubdomainsPath = params.template getParameter< std::string >( "subdomains-config" );
   for( int subset = 0; subset < this->numberOfSubsets; subset++ )
      TNL::SPH::configSubdomain( subset, this->configSubdomains );
   parseDistributedConfig( configSubdomainsPath, parametersSubdomains, configSubdomains, log );

   topology.loadFromConfig( params, parametersSubdomains );
   topology.finalizeLinear();
   initParticleSets();
   initMultiResolutionBoundaryPatches();
   this->timeMeasurement.addTimer( "multiresolution-update" );

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
   for( int i = 0; i < this->numberOfSubsets; i++ ){
      std::string subdomainKey = "subdomain-" + std::to_string( i ) + "-";
      const float refinementFactor = parametersSubdomains.getParameter< float >( subdomainKey + "refinement-factor" );
      multiresolutionBoundaryPatches[ i ]->initMassNodes( this->modelParams, i, refinementFactor );
   }

   this->timeStepping.setTimeStep( params.template getParameter< RealType >( "initial-time-step" ) );
   this->timeStepping.setEndTime( params.template getParameter< RealType >( "final-time" ) );
   this->timeStepping.addOutputTimer( "save_results", params.template getParameter< RealType >( "snapshot-period" ) );

   readParticlesFiles();

   log.writeSeparator();
   if( params.template getParameter< std::string >( "measuretool-config" ) != "" ) {
      log.writeParameter( "Simulation monitor initialization.", "" );
      this->simulationMonitor.init( params, this->timeStepping, log );
      log.writeParameter( "Simulation monitor initialization.", "Done." );
   }
}

#ifdef HAVE_MPI
template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::initializeDistributedSimulation()
{
   auto& params = this->parameters;
   auto& log = this->logger;

   this->numberOfSubsets = 1;
   this->fluidSets.resize( 1 );
   this->boundarySets.resize( 1 );

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

   readParticleFilesDistributed( params, this->parametersDistributed, log );

   log.writeSeparator();
   if( params.template getParameter< std::string >( "measuretool-config" ) != "" ) {
      log.writeParameter( "Simulation monitor initialization.", "" );
      this->simulationMonitor.init( params, this->timeStepping, log );
      log.writeParameter( "Simulation monitor initialization.", "Done." );
   }
}
#endif

template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::initParticleSets()
{
   const int nSubsets = topology.getNumberOfSubdomains();
   this->numberOfSubsets = nSubsets;
   this->fluidSets.resize( nSubsets );
   this->boundarySets.resize( nSubsets );
   multiresolutionBoundaryPatches.resize( nSubsets );

   for( int i = 0; i < nSubsets; i++ ){
      std::string subdomainKey = "subdomain-" + std::to_string( i ) + "-";
      this->fluidSets[ i ]->initializeAsDistributed(
            parametersSubdomains.getParameter< int >( subdomainKey + "fluid_n" ),
            parametersSubdomains.getParameter< int >( subdomainKey + "fluid_n_allocated" ),
            topology.getLocalGrid( i ),
            topology.getLocalOriginCoordinates( i ),
            topology.getGlobalGrid(),
            topology.getNumberOfOverlapLayers() );
      if constexpr( ParticlesType::specifySearchedSetExplicitly() == true )
         this->fluidSets[ i ]->getParticles()->setParticleSetLabel( 0 );

      this->boundarySets[ i ]->initializeAsDistributed(
            parametersSubdomains.getParameter< int >( subdomainKey + "boundary_n" ),
            parametersSubdomains.getParameter< int >( subdomainKey + "boundary_n_allocated" ),
            topology.getLocalGrid( i ),
            topology.getLocalOriginCoordinates( i ),
            topology.getGlobalGrid(),
            topology.getNumberOfOverlapLayers() );
      if constexpr( ParticlesType::specifySearchedSetExplicitly() == true )
         this->boundarySets[ i ]->getParticles()->setParticleSetLabel( 1 );
   }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::initMultiResolutionBoundaryPatches()
{
   const int mrbCount = topology.getNumberOfSubdomainInterfaces();
   multiresolutionBoundaryPatches.resize( mrbCount );

   int mrbIdx = 0;
   for( int i = 0; i < topology.getNumberOfSubdomains(); i++ ) {
      const std::string key = "subdomain-" + std::to_string( i ) + "-";
      const float rf = parametersSubdomains.getParameter< float >( key + "refinement-factor" );

      for( const auto& iface : topology.getInterfacesOfSubdomain( i ) ) {

         multiresolutionBoundaryPatches[ mrbIdx ]->initializeAsDistributed(
            0,
            60000,
            topology.getLocalGrid( i ),
            topology.getLocalOriginCoordinates( i ),
            topology.getGlobalGrid(),
            topology.getNumberOfOverlapLayers() );

         multiresolutionBoundaryPatches[ mrbIdx ]->initZones(
            this->fluidSets[ iface.ownIdx ]->getParticles(),
            this->fluidSets[ iface.neighborIdx ]->getParticles(),
            rf );

         mrbIdx++;
      }
   }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::readParticlesFiles()
{
   auto& params = this->parameters;
   auto& log = this->logger;
   auto& paramsOB = this->parametersOpenBoundary;

   const int nSubsets = params.template getParameter< int >( "numberOfSubdomains" );

   for( int i = 0; i < nSubsets; i++ ){
      std::string subdomainKey = "subdomain-" + std::to_string( i ) + "-";
      if( parametersSubdomains.getParameter< int >( subdomainKey + "fluid_n" ) != 0 ){
         const std::string fluidFileName = parametersSubdomains.getParameter< std::string >( subdomainKey + "fluid-particles" );
         log.writeParameter( "Reading fluid particles:", fluidFileName );
         this->fluidSets[ i ]->template readParticlesAndVariables< typename BaseType::SimulationReaderType >( fluidFileName );
      }
      if( parametersSubdomains.getParameter< int >( subdomainKey + "boundary_n" ) != 0 ){
         const std::string boundaryFileName = parametersSubdomains.getParameter< std::string >( subdomainKey + "boundary-particles" );
         log.writeParameter( "Reading boundary particles:", boundaryFileName );
         this->boundarySets[ i ]->template readParticlesAndVariables< typename BaseType::SimulationReaderType >( boundaryFileName );
      }
   }

   const int numberOfBoundaryPatches = params.template getParameter< int >( "openBoundaryPatches" );
   if( numberOfBoundaryPatches > 0 ) {
      for( int i = 0; i < numberOfBoundaryPatches; i++ ) {
         std::string prefix = "buffer-" + std::to_string( i + 1 ) + "-";
         log.writeParameter( "Reading open boundary particles:", paramsOB.template getParameter< std::string >( prefix + "particles" ) );
         this->openBoundaryPatches[ i ]->template readParticlesAndVariables< typename BaseType::SimulationReaderType >(
            paramsOB.template getParameter< std::string >( prefix + "particles" ) );
      }
   }
}

#ifdef HAVE_MPI
template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::initDistributedParticleSets( TNL::Config::ParameterContainer& parameters,
                                                                          TNL::Config::ParameterContainer& parametersDistributed,
                                                                          TNL::Logger& logger )
{
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

template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::readParticleFilesDistributed( TNL::Config::ParameterContainer& parameters,
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
SolverMultiSetBlockMultiresolution< Model >::multiresolutionUpdate()
{
   this->timeMeasurement.start( "multiresolution-update" );

   multiresolutionBoundaryPatches[ 0 ]->updateInterfaceBuffer(
         this->fluidSets[ 0 ], this->fluidSets[ 1 ], this->modelParams, this->timeStepping.getTimeStep(), 0 );
   multiresolutionBoundaryPatches[ 1 ]->updateInterfaceBuffer(
         this->fluidSets[ 1 ], this->fluidSets[ 0 ], this->modelParams, this->timeStepping.getTimeStep(), 1 );

   this->timeMeasurement.stop( "multiresolution-update" );
   this->writeLog( "Update multiresolution BC...", "Done.");
}

template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::interact()
{
   this->timeMeasurement.start( "interact" );

   for( int i = 0; i < this->numberOfSubsets; i++ ){
      this->model.updateSolidBoundary( this->fluidSets[ i ], this->boundarySets[ i ], this->modelParams );
      this->model.finalizeBoundaryInteraction( this->fluidSets[ i ], this->boundarySets[ i ], this->modelParams );

      this->model.interaction( this->fluidSets[ i ], this->boundarySets[ i ], this->modelParams );
      this->model.interactionWithOpenBoundary( this->fluidSets[ i ], multiresolutionBoundaryPatches[ i ], this->modelParams );
      this->model.finalizeInteraction( this->fluidSets[ i ], this->boundarySets[ i ], this->modelParams );
   }

   this->timeMeasurement.stop( "interact" );
   this->writeLog( "Interact...", "Done." );
}

template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::save( bool writeParticleCellIndex )
{
   if( ( this->verbose == "with-snapshot" ) || ( this->verbose == "full" ) )
      this->writeInfo();

   const RealType time = this->timeStepping.getTime();

   for( int i = 0; i < this->numberOfSubsets; i++ ){
#ifdef HAVE_MPI
      std::string outputFileNameFluid = this->outputDirectory + "/fluid_rank" + std::to_string( TNL::MPI::GetRank() ) + "_" + std::to_string( time ) + "_particles.vtk";
#else
      std::string outputFileNameFluid = this->outputDirectory + "/fluid_subdomain" + std::to_string( i ) + "_" + std::to_string( time ) + "_particles.vtk";
#endif
      this->fluidSets[ i ]->template writeParticlesAndVariables< typename BaseType::Writer >( outputFileNameFluid, writeParticleCellIndex );
      this->logger.writeParameter( "Saved:", outputFileNameFluid );

#ifdef HAVE_MPI
      std::string outputFileNameBound = this->outputDirectory + "/boundary_rank" + std::to_string( TNL::MPI::GetRank() ) + "_" + std::to_string( time ) + "_particles.vtk";
#else
      std::string outputFileNameBound = this->outputDirectory + "/boundary_subdomain" + std::to_string( i ) + "_" + std::to_string( time ) + "_particles.vtk";
#endif
      this->boundarySets[ i ]->template writeParticlesAndVariables< typename BaseType::Writer >( outputFileNameBound, writeParticleCellIndex );
      this->logger.writeParameter( "Saved:", outputFileNameBound );

      if( multiresolutionBoundaryPatches.size() ) {
         std::string outputFileNameMultiresolutionBound =
            this->outputDirectory + "/multiresolutionBoundaryPatch_subdomain" + std::to_string( i ) + "_" + std::to_string( time ) + "_particles.vtk";
         multiresolutionBoundaryPatches[ i ]->template writeParticlesAndVariables< typename BaseType::Writer >(
               outputFileNameMultiresolutionBound, writeParticleCellIndex );
         this->logger.writeParameter( "Saved:", outputFileNameMultiresolutionBound );
      }

#ifdef HAVE_MPI
      std::string outputFileNameGrid = this->outputDirectory + "/grid_rank" + std::to_string( TNL::MPI::GetRank() + 1 ) + "_" + std::to_string( time ) + ".vtk";
#else
      std::string outputFileNameGrid = this->outputDirectory + "/grid_subdomain" + std::to_string( i ) + "_" + std::to_string( time ) + ".vtk";
#endif
      TNL::Writers::writeBackgroundGrid( outputFileNameGrid, this->fluidSets[ i ]->getParticles()->getGridDimensions(), this->fluidSets[ i ]->getParticles()->getGridOrigin(), this->fluidSets[ i ]->getParticles()->getSearchRadius() );
      this->logger.writeParameter( "Saved:", outputFileNameGrid );

      this->simulationMonitor.save( this->logger );
   }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::writeProlog( bool writeSystemInformation ) noexcept
{
   auto& log = this->logger;

   log.writeHeader( "SPH simulation configuration." );
   log.writeParameter( "Case name:", this->caseName );
   log.writeSeparator();

   const bool printGPUs = std::is_same< DeviceType, TNL::Devices::Cuda >::value;

   if( TNL::MPI::isInitialized() )
      log.writeParameter( "MPI processes:", TNL::MPI::GetSize() );
   log.writeParameter( "Device type:", TNL::getType< DeviceType >() );
   if( ! printGPUs ) {
      if( TNL::Devices::Host::isOMPEnabled() ) {
         log.writeParameter( "OMP enabled:", "yes", 1 );
         log.writeParameter( "OMP threads:", TNL::Devices::Host::getMaxThreadsCount(), 1 );
      }
      else
         log.writeParameter( "OMP enabled:", "no", 1 );
   }

   log.writeParameter( "Particle system type:", Model::ParticlesType::writeModelType() );
   log.writeParameter( "SPH model:", Model::writeModelType() );
   log.writeParameter( "Verbose:", this->verbose );
   log.writeParameter( "Output directory:", this->outputDirectory );
   log.writeParameter( "Particles format", this->particlesFormat );
   writePrologModel( log, this->modelParams );

   for( int i = 0; i < this->numberOfSubsets; i++ ){
      log.writeHeader( "Fluid " + std::to_string( i ) + " object information." );
      this->fluidSets[ i ]->writeProlog( log );
      log.writeHeader( "Boundary " + std::to_string( i ) + " object information:" );
      this->boundarySets[ i ]->writeProlog( log );
   }

   if( this->openBoundaryPatches.size() > 0 ) {
      for( long unsigned int i = 0; i < this->openBoundaryPatches.size(); i++ ) {
         log.writeHeader( "Open boundary buffer" + std::to_string( i + 1 ) + "." );
         this->openBoundaryPatches[ i ]->writeProlog( log );
         log.writeSeparator();
         this->openBoundaryPatches[ i ]->config.writeProlog( log );
      }
   }

   if( multiresolutionBoundaryPatches.size() > 0 ) {
      for( long unsigned int i = 0; i < multiresolutionBoundaryPatches.size(); i++ ) {
         log.writeHeader( "Multiresolution boundary buffer" + std::to_string( i + 1 ) + "." );
         multiresolutionBoundaryPatches[ i ]->writeProlog( log, i );
      }
   }

   log.writeHeader( "System information." );
   if( writeSystemInformation ) {
      log.writeSystemInformation( printGPUs );
      log.writeSeparator();
      log.writeCurrentTime( "Started at:" );
      log.writeSeparator();
   }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::writeInfo() noexcept
{
   auto& log = this->logger;

   log.writeSeparator();
   log.writeParameter( "Simulation time: " + std::to_string( this->timeStepping.getTime() )
                                  + " s, simulation step: " + std::to_string( this->timeStepping.getStep() ),
                               "" );
   log.writeCurrentTime( "Current time:" );
   for( int i = 0; i < this->numberOfSubsets; i++ ){
      log.writeParameter( "Number of fluid particles:", this->fluidSets[ i ]->getNumberOfParticles() );
      log.writeParameter( "Number of allocated fluid particles:", this->fluidSets[ i ]->getNumberOfAllocatedParticles() );
      log.writeParameter( "Number of boundary particles:", this->boundarySets[ i ]->getNumberOfParticles() );
      log.writeParameter( "Number of allocated boundary particles:", this->boundarySets[ i ]->getNumberOfAllocatedParticles() );
   }
   if( this->openBoundaryPatches.size() > 0 ) {
      for( long unsigned int i = 0; i < this->openBoundaryPatches.size(); i++ )
         log.writeParameter( "Number of buffer" + std::to_string( i + 1 ) + " particles:",
                                this->openBoundaryPatches[ i ]->getNumberOfParticles() );
   }
   if( multiresolutionBoundaryPatches.size() > 0 ) {
      for( long unsigned int i = 0; i < multiresolutionBoundaryPatches.size(); i++ )
         log.writeParameter( "Number of mr-buffer " + std::to_string( i ) + " particles:",
               multiresolutionBoundaryPatches[ i ]->getNumberOfParticles() );
   }
   log.writeSeparator();
}

} // namespace SPH
} // namespace TNL
