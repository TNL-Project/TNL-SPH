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

#include "tempFunctionsToConfigMultiresolution.h"

#include <SPH/configInit.h>

namespace TNL {
namespace SPH {

template< typename Model >
void
SolverMultiSet< Model >::init( int argc, char* argv[] )
{
   try {
      initialize< SimulationType >( argc, argv, cliParams, cliConfig, parameters, config );
   }
   catch ( ... ) {
      std::cerr << std::endl;
   }

   logger.writeHeader( "SPH simulation initialization." );
#ifdef HAVE_MPI
   // build config for distributed domain
   logger.writeParameter( "Configuration of distributed simulation:", "" );
   Containers::StaticVector< 2, int > numberOfSubdomains = parameters.getXyz< Containers::StaticVector< 2, int > >( "subdomains" );
   logger.writeParameter( "Number of subdomains:", numberOfSubdomains );
   for( int x = 0; x < numberOfSubdomains[ 0 ]; x++ )
      for( int y = 0; y < numberOfSubdomains[ 1 ]; y++ )
         TNL::SPH::configSetupDistributedSubdomain( x, y, this->configDistributed );
   // load and parse config for distributed domain
   std::string configDistributedPath = parameters.getParameter< std::string >( "distributed-config" );
   parseDistributedConfig( configDistributedPath, parametersDistributed, configDistributed, logger );

   // initialize distributed particle sets and overlaps
   initDistributedParticleSets( parameters, this->parametersDistributed, logger );

   // set balancing tresholds
   loadBalancingMeasure = parameters.getParameter< std::string >( "load-balancing-measure" );
   loadBalancingStepInterval = parameters.getParameter< int >( "load-balancing-step-inteval" );
   fluid->getDistributedParticles()->setParticlesCountResizeTrashold(
         parameters.getParameter< float >( "number-of-particles-balancing-coef" ) );
   fluid->getDistributedParticles()->setCompTimeResizePercetnageTrashold(
         parameters.getParameter< float >( "computational-time-balancing-coef" ) );

   // add timers related to synchronization and load balancing
   timeMeasurement.addTimer( "synchronize", false );
   timeMeasurement.addTimer( "rebalance", false );
#else
   this->numberOfSubsets = parameters.getParameter< int >( "numberOfSubdomains" );
   const std::string configSubdomainsPath = parameters.getParameter< std::string >( "subdomains-config" );
   for( int subset = 0; subset < numberOfSubsets; subset++ )
      TNL::SPH::configSubdomain( subset, this->configSubdomains );
   parseDistributedConfig( configSubdomainsPath, parametersSubdomains, configSubdomains, logger ); //FIXME: Wrong function

   // initialize topology
   topology.loadFromConfig( parameters, parametersSubdomains );
   topology.finalizeLinear(); //FIXME: Do only if linear, maybe hide this into the topology
   initParticleSets();
   initMultiResolutionBundaryPatches();
   timeMeasurement.addTimer( "multiresolution-update" );
#endif

   // initialize open boundary conditions
   if( parameters.getParameter< std::string >( "open-boundary-config" ) != "" ){
      initOpenBoundaryPatches( parameters, logger );

      // add custom timers related to open boundary conditions
      timeMeasurement.addTimer( "extrapolate-openbc" );
      timeMeasurement.addTimer( "apply-openbc" );
   }

   // init periodic boundary conditions
   if( parameters.getParameter< std::string >( "periodic-boundary-config" ) != "" ){
      initPeriodicBoundaryPatches( parameters, logger );

      // add custom timers related to perioric boundary conditions
      timeMeasurement.addTimer( "enforce-periodic-bc" );
      timeMeasurement.addTimer( "transfer-periodic-bc" );
      timeMeasurement.addTimer( "periodicity-fluid-updateZone", false );
      timeMeasurement.addTimer( "periodicity-boundary-updateZone", false );
   }

   // init model parameters
   modelParams.init( parameters );
   // FIXME: See note in initializeParticleSets()
   for( int i = 0; i < numberOfSubsets; i++ ){
      std::string subdomainKey = "subdomain-" + std::to_string( i ) + "-";
      const float refinementFactor = parametersSubdomains.getParameter< float >( subdomainKey + "refinement-factor" );
      multiresolutionBoundaryPatches[ i ]->initMassNodes( modelParams, i, refinementFactor );
   }

   // init time stepping
   timeStepping.setTimeStep( parameters.getParameter< RealType >( "initial-time-step" ) );
   timeStepping.setEndTime( parameters.getParameter< RealType >( "final-time" ) );
   timeStepping.addOutputTimer( "save_results", parameters.getParameter< RealType >( "snapshot-period" ) );

   // control
   caseName = parameters.getParameter< std::string >( "case-name" );
   verbose = parameters.getParameter< std::string >( "verbose-intensity" );
   outputDirectory = parameters.getParameter< std::string >( "output-directory" );
   particlesFormat = parameters.getParameter< std::string >( "particles-format" );

#ifdef HAVE_MPI
   readParticleFilesDistributed( parameters, this->parametersDistributed, logger );
#else
   readParticlesFiles();
#endif

   // initialize the measuretool
   logger.writeSeparator();
   if( parameters.getParameter< std::string >( "measuretool-config" ) != "" ) {
      logger.writeParameter( "Simulation monitor initialization.", "" );
      simulationMonitor.init( parameters, timeStepping, logger );
      logger.writeParameter( "Simulation monitor initialization.", "Done." );
   }

   //#ifdef __CUDACC__
   //TNL::Pointers::synchronizeSmartPointersOnDevice< DeviceType >();
   //#endif

   logger.writeHeader( "SPH simulation successfully initialized." );
}


template< typename Model >
void
SolverMultiSet< Model >::initParticleSets()
{
   const int numberOfSubsets = topology.getNumberOfSubdomains();
   fluidSets.resize( numberOfSubsets );
   boundarySets.resize( numberOfSubsets );
   multiresolutionBoundaryPatches.resize( numberOfSubsets );

   for( int i = 0; i < numberOfSubsets; i++ ){
      std::string subdomainKey = "subdomain-" + std::to_string( i ) + "-";
      // init fluid
      fluidSets[ i ]->initializeAsDistributed(
            parametersSubdomains.getParameter< int >( subdomainKey + "fluid_n" ),
            parametersSubdomains.getParameter< int >( subdomainKey + "fluid_n_allocated" ),
            topology.getLocalGrid( i ),
            topology.getLocalOriginCoordinates( i ), //TODO: Breaks the structure
            topology.getGlobalGrid(),
            topology.getNumberOfOverlapLayers() );
      if constexpr( ParticlesType::specifySearchedSetExplicitly() == true )
         fluidSets[ i ]->getParticles()->setParticleSetLabel( 0 );

      // init boundary
      boundarySets[ i ]->initializeAsDistributed(
            parametersSubdomains.getParameter< int >( subdomainKey + "boundary_n" ),
            parametersSubdomains.getParameter< int >( subdomainKey + "boundary_n_allocated" ),
            topology.getLocalGrid( i ),
            topology.getLocalOriginCoordinates( i ), //TODO: Breaks the structure
            topology.getGlobalGrid(),
            topology.getNumberOfOverlapLayers() );
      if constexpr( ParticlesType::specifySearchedSetExplicitly() == true )
         boundarySets[ i ]->getParticles()->setParticleSetLabel( 1 );
   }
}

template< typename Model >
void
SolverMultiSet< Model >::initMultiResolutionBundaryPatches()
{
   const int mrbCount = topology.getNumberOfSubdomainInterfaces();
   multiresolutionBoundaryPatches.resize( mrbCount );

   int mrbIdx = 0;
   for( int i = 0; i < topology.getNumberOfSubdomains(); i++ ) {
      const std::string key = "subdomain-" + std::to_string( i ) + "-";
      const float rf = parametersSubdomains.getParameter< float >( key + "refinement-factor" );

      for( const auto& iface : topology.getInterfacesOfSubdomain( i ) ) {

         multiresolutionBoundaryPatches[ mrbIdx ]->initializeAsDistributed(
            0, //FIXME: Create some estimate for the alloc size
            60000, //FIXME: Create some estimate for the alloc size
            topology.getLocalGrid( i ),
            topology.getLocalOriginCoordinates( i ), //TODO: Breaks the structure
            topology.getGlobalGrid(),
            topology.getNumberOfOverlapLayers() );

         multiresolutionBoundaryPatches[ mrbIdx ]->initZones(
            fluidSets[ iface.ownIdx ]->getParticles(),
            fluidSets[ iface.neighborIdx ]->getParticles(),
            rf ); //FIXME: Add max number of particles per cell, compute from params, I also dont need rf

         mrbIdx++;
      }
   }
}

template< typename Model >
void
SolverMultiSet< Model >::initOpenBoundaryPatches( TNL::Config::ParameterContainer& parameters, TNL::Logger& logger  )
{
   logger.writeParameter( "Initialization of open boundary patches.", "" );
   const int numberOfBoundaryPatches = parameters.getParameter< int >( "openBoundaryPatches" );
   const std::string openBoundaryConfigPath = parameters.getParameter< std::string >( "open-boundary-config" );

   // setup and parse open boundary config
   for( int i = 0; i < numberOfBoundaryPatches; i++ ) {
      std::string prefix = "buffer-" + std::to_string( i + 1 ) + "-";
      configSetupOpenBoundaryModelPatch< SPHConfig >( configOpenBoundary, prefix );
   }
   parseOpenBoundaryConfig( openBoundaryConfigPath, parametersOpenBoundary, configOpenBoundary, logger );

   // get global domain properetis
   const VectorType domainOrigin = parameters.getXyz< VectorType >( "domainOrigin" );
   const VectorType domainSize = parameters.getXyz< VectorType >( "domainSize" );
   const RealType searchRadius = parameters.getParameter< RealType >( "searchRadius" );
   const IndexVectorType gridSize = TNL::ceil( ( domainSize - domainOrigin ) / searchRadius );

   openBoundaryPatches.resize( numberOfBoundaryPatches );
   for( int i = 0; i < numberOfBoundaryPatches; i++ ) {
      std::string prefix = "buffer-" + std::to_string( i + 1 ) + "-";
      openBoundaryPatches[ i ]->config.init( parameters, parametersOpenBoundary, prefix );
      openBoundaryPatches[ i ]->initialize( parametersOpenBoundary.getParameter< int >( prefix + "numberOfParticles" ),
                                            parametersOpenBoundary.getParameter< int >( prefix + "numberOfAllocatedParticles" ),
                                            searchRadius,
                                            gridSize,
                                            domainOrigin );
   }
   logger.writeParameter( "Initialization of open boundary patches.", "Done." );
}

template< typename Model >
void
SolverMultiSet< Model >::initPeriodicBoundaryPatches( TNL::Config::ParameterContainer& parameters, TNL::Logger& logger )
{
   //FIXME: Update to multiset: logger.writeParameter( "Initialization of open boundary patches.", "" );
   //FIXME: Update to multiset: const int numberOfBoundaryPatches = parameters.getParameter< int >( "periodicBoundaryPatches" );
   //FIXME: Update to multiset: const std::string openBoundaryConfigPath = parameters.getParameter< std::string >( "periodic-boundary-config" );

   //FIXME: Update to multiset: // setup and parse open boundary config
   //FIXME: Update to multiset: for( int i = 0; i < numberOfBoundaryPatches; i++ ) {
   //FIXME: Update to multiset:    std::string prefix = "buffer-" + std::to_string( i + 1 ) + "-";
   //FIXME: Update to multiset:    configSetupOpenBoundaryModelPatch< SPHConfig >( configOpenBoundary, prefix );
   //FIXME: Update to multiset: }
   //FIXME: Update to multiset: parseOpenBoundaryConfig( openBoundaryConfigPath, parametersOpenBoundary, configOpenBoundary, logger );

   //FIXME: Update to multiset: fluid->initializePeriodicity( parameters, parametersOpenBoundary );
   //FIXME: Update to multiset: boundary->initializePeriodicity( parameters, parametersOpenBoundary );
}

#ifdef HAVE_MPI
template< typename Model >
void
SolverMultiSet< Model >::initDistributedParticleSets( TNL::Config::ParameterContainer& parameters,
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

   // global domain properties
   const RealType searchRadius = parameters.getParameter< RealType >( "searchRadius" );
   const VectorType domainOrigin = parameters.getXyz< VectorType >( "domainOrigin" );
   const VectorType domainSize = parameters.getXyz< VectorType >( "domainSize" );
   const IndexVectorType domainGridDimension = TNL::ceil( ( domainSize - domainOrigin ) / searchRadius );
   const int numberOfOverlapLayers = parameters.getParameter< int >( "overlapWidth" );

   // subdomain + ghost properties
   const VectorType subdomainOrigin = parametersDistributed.getXyz< VectorType >( subdomainKey + "origin" ); //REMOVE
   const IndexVectorType subdomainGridDimension = parametersDistributed.getXyz< IndexVectorType >( subdomainKey + "grid-dimensions" );
   const IndexVectorType subdomainGridOriginGlobalCoords = parametersDistributed.getXyz< IndexVectorType >( subdomainKey + "origin-global-coords" );

   // init fluid
   logger.writeParameter( "initDistributed:", "fluid->initialize" );
   fluid->initializeAsDistributed( parametersDistributed.getParameter< int >( subdomainKey + "fluid_n" ),
                                   parametersDistributed.getParameter< int >( subdomainKey + "fluid_n_allocated" ),
                                   searchRadius,
                                   domainGridDimension,
                                   domainOrigin,
                                   subdomainGridDimension,
                                   subdomainGridOriginGlobalCoords,
                                   numberOfOverlapLayers,
                                   numberOfSubdomains,
                                   subdomainOrigin, //REMOVE
                                   logger );
   // since we use multiple set, we need to rewrite the default communicator with the one provided by distributed solver
   fluid->getDistributedParticles()->writeProlog( logger );
   fluid->setCommunicator( this->communicator );
   fluid->getDistributedParticlesSynchronizer().initialize( fluid->getDistributedParticles() ); //FIXME

   // init boundary
   logger.writeParameter( "initDistributed:", "boundary->initialize" );
   boundary->initializeAsDistributed( parametersDistributed.getParameter< int >( subdomainKey + "boundary_n" ),
                                      parametersDistributed.getParameter< int >( subdomainKey + "boundary_n_allocated" ),
                                      searchRadius,
                                      domainGridDimension,
                                      domainOrigin,
                                      subdomainGridDimension,
                                      subdomainGridOriginGlobalCoords,
                                      numberOfOverlapLayers,
                                      numberOfSubdomains,
                                      subdomainOrigin, //REMOVE
                                      logger );
   // since we use multiple set, we need to rewrite the default communicator with the one provided by distributed solver
   boundary->getDistributedParticles()->writeProlog( logger );
   boundary->setCommunicator( this->communicator );
   boundary->getDistributedParticlesSynchronizer().initialize( boundary->getDistributedParticles() ); //FIXME
}
#endif

template< typename Model >
void
SolverMultiSet< Model >::readParticlesFiles()
{
   const int numberOfSubsets = parameters.getParameter< int >( "numberOfSubdomains" );

   for( int i = 0; i < numberOfSubsets; i++ ){
      std::string subdomainKey = "subdomain-" + std::to_string( i ) + "-";
      if( parametersSubdomains.getParameter< int >( subdomainKey + "fluid_n" ) != 0 ){
         const std::string fluidFileName = parametersSubdomains.getParameter< std::string >( subdomainKey + "fluid-particles" );
         logger.writeParameter( "Reading fluid particles:", fluidFileName );
         fluidSets[ i ]->template readParticlesAndVariables< SimulationReaderType >( fluidFileName );
      }
      const std::string boundaryFileName = parametersSubdomains.getParameter< std::string >( subdomainKey + "boundary-particles" );
      logger.writeParameter( "Reading boundary particles:", boundaryFileName );
      boundarySets[ i ]->template readParticlesAndVariables< SimulationReaderType >( boundaryFileName );
   }

   // init open boundary patches
   const int numberOfBoundaryPatches = parameters.getParameter< int >( "openBoundaryPatches" );
   //TODO: I don't like that open boundary buffer are selected in compile time.
   if( numberOfBoundaryPatches > 0 ) {
      for( int i = 0; i < numberOfBoundaryPatches; i++ ) {
         std::string prefix = "buffer-" + std::to_string( i + 1 ) + "-";
         logger.writeParameter( "Reading open boundary particles:", parametersOpenBoundary.getParameter< std::string >( prefix + "particles" ) );
         openBoundaryPatches[ i ]->template readParticlesAndVariables< SimulationReaderType >(
            parametersOpenBoundary.getParameter< std::string >( prefix + "particles" ) );
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

   // read particle data
   logger.writeParameter( "Reading fluid particles:", parametersDistributed.getParameter< std::string >( subdomainKey + "fluid-particles" ) );
   fluid->template readParticlesAndVariables< SimulationReaderType >(
      parametersDistributed.getParameter< std::string >( subdomainKey + "fluid-particles" ) );
   logger.writeParameter( "Reading boundary particles:", parametersDistributed.getParameter< std::string >( subdomainKey + "boundary-particles" ) );
   boundary->template readParticlesAndVariables< SimulationReaderType >(
      parametersDistributed.getParameter< std::string >( subdomainKey + "boundary-particles" ) );

   const int numberOfBoundaryPatches = parameters.getParameter< int >( "openBoundaryPatches" );
   if( openBoundaryPatches.size() > 0 ) {  //TODO: I dont like this.
      for( int i = 0; i < numberOfBoundaryPatches; i++ ) {
         std::string prefix = subdomainKey + "buffer-" + std::to_string( i + 1 ) + "-";
         openBoundaryPatches[ i ]->template readParticlesAndVariables< SimulationReaderType >(
            parametersDistributed.getParameter< std::string >( prefix + "particles" ) );
      }
   }
}
#endif

template< typename Model >
void
SolverMultiSet< Model >::performNeighborSearch( bool performBoundarySearch )
{
   timeMeasurement.start( "search" );
   for( int i = 0; i < numberOfSubsets; i++ ){
      std::cout << "Particle set: " <<  i << std::endl;
      if constexpr( ParticlesType::specifySearchedSetExplicitly() == false ){
         fluidSets[ i ]->searchForNeighbors();
         writeLog( "Fluid search procedure:", "Done." );

         if( timeStepping.getStep() == 0 || performBoundarySearch == true ){
            boundarySets[ i ]->searchForNeighbors();
            writeLog( "Boundary search procedure:", "Done." );
         }
         multiresolutionBoundaryPatches[ i ]->searchForNeighbors();
         writeLog( "Multiresolution patch search procedure:", "Done." );

      }
      else if constexpr( ParticlesType::specifySearchedSetExplicitly() == true ){
         // FIXME: Support
      }
   }
   std::cout << "Fluid search - done." << std::endl;

   // search fluid patches
   if constexpr( ParticlesType::specifySearchedSetExplicitly() == false ){
         if( openBoundaryPatches.size() > 0 )
            for( auto& openBoundaryPatch : openBoundaryPatches ){
               openBoundaryPatch->searchForNeighbors();
               writeLog( "Open boundary patch search procedure:", "Done." );
            }
   }
   else if constexpr( ParticlesType::specifySearchedSetExplicitly() == true ){
      // FIXME: Support
   }
   timeMeasurement.stop( "search" );
   writeLog("Search...", "Done." );
}

template< typename Model >
void
SolverMultiSet< Model >::removeParticlesOutOfDomain()
{
   for( int i = 0; i < numberOfSubsets; i++ ){
      std::cout << "Remove particles out of domain... ";
      const int numberOfParticlesToRemove = fluidSets[ i ]->getParticles()->getNumberOfParticlesToRemove();
      fluidSets[ i ]->getParticles()->removeParitclesOutOfDomain();

      if( fluidSets[ i ]->getParticles()->getNumberOfParticlesToRemove() > numberOfParticlesToRemove ){
         const int numberOfParticlesOutOfDomain = fluidSets[ i ]->getParticles()->getNumberOfParticlesToRemove() - numberOfParticlesToRemove;
         writeLog( "Number of out of domain removed particles:", numberOfParticlesOutOfDomain  );
         // search for neighbros
         //FIXME: I can not search dist
         //timeMeasurement.start( "search" );
         //this->performNeighborSearch( log );
         //timeMeasurement.stop( "search" );
      }
      std::cout << "... done." << std::endl;
   }
}

template< typename Model >
template< typename ParticleSetPointer >
void
SolverMultiSet< Model >::performNeighborSearchForObject( ParticleSetPointer& objectPointer )
{
   objectPointer->getParticles()->resetListWithIndices();
   objectPointer->getParticles()->computeParticleCellIndices();
   objectPointer->sortParticles();
   objectPointer->getParticles()->particlesToCells();
}

template< typename Model >
void
SolverMultiSet< Model >::extrapolateOpenBC()
{
   timeMeasurement.start( "extrapolate-openbc" );
   //FIXME: for( long unsigned int i = 0; i < std::size( openBoundaryPatches ); i++ ) {
   //FIXME:    //TODO Check if open boundary buffer is really open boundary buffer
   //FIXME:    model.extrapolateOpenBoundaryData( fluid, openBoundaryPatches[ i ], modelParams, openBoundaryPatches[ i ]->config );
   //FIXME: }
   timeMeasurement.stop( "extrapolate-openbc" );
   writeLog( "Extrapolate open BC...", "Done." );
}

template< typename Model >
void
SolverMultiSet< Model >::applyOpenBC( const RealType timeStepFact )
{
   timeMeasurement.start( "apply-openbc" );
   //FIXME: for( long unsigned int i = 0; i < std::size( openBoundaryPatches ); i++ ) {
   //FIXME:    //TODO Check if open boundary buffer is really open boundary buffer
   //FIXME:    openBoundaryModel.applyOpenBoundary(
   //FIXME:       timeStepFact * timeStepping.getTimeStep(), fluid, openBoundaryPatches[ i ], openBoundaryPatches[ i ]->config );
   //FIXME: }
   timeMeasurement.stop( "apply-openbc" );
   writeLog( "Update open BC...", "Done." );
}

template< typename Model >
void
SolverMultiSet< Model >::applyPeriodicBCEnforce()
{
   //FIXME: timeMeasurement.start( "periodicity-fluid-updateZone" );
   //FIXME: fluid->enforcePeriodicPatches();
   //FIXME: timeMeasurement.stop( "periodicity-boundary-updateZone" );
   //FIXME: timeMeasurement.start( "periodicity-fluid-updateZone" );
   //FIXME: boundary->enforcePeriodicPatches();
   //FIXME: timeMeasurement.stop( "periodicity-boundary-updateZone" );
}

template< typename Model >
void
SolverMultiSet< Model >::applyPeriodicBCTransfer()
{
   //FIXME: const long unsigned int numberOfPeriodicPatches = std::size( fluid->periodicPatches );
   //FIXME: for( long unsigned int i = 0; i < numberOfPeriodicPatches; i++ ) {
   //FIXME:    openBoundaryModel.periodicityParticleTransfer( fluid, fluid->periodicPatches[ i ] );
   //FIXME: }
}

template< typename Model >
void
SolverMultiSet< Model >::multiresolutionUpdate()
{
   //FIXME: Describe by loops
   timeMeasurement.start( "multiresolution-update" );

   multiresolutionBoundaryPatches[ 0 ]->updateInterfaceBuffer(
         fluidSets[ 0 ], fluidSets[ 1 ], modelParams, timeStepping.getTimeStep(), 0 );
   multiresolutionBoundaryPatches[ 1 ]->updateInterfaceBuffer(
         fluidSets[ 1 ], fluidSets[ 0 ], modelParams, timeStepping.getTimeStep(), 1 );

   timeMeasurement.stop( "multiresolution-update" );
   writeLog( "Update multiresolution BC...", "Done.");
}

template< typename Model >
void
SolverMultiSet< Model >::interact()
{
   timeMeasurement.start( "interact" );

   for( int i = 0; i < numberOfSubsets; i++ ){
      // update solid boundary conditions
      //std::cout << "Set i: " << i << " updating boundary." << std::endl;
      model.updateSolidBoundary( fluidSets[ i ], boundarySets[ i ], modelParams );
      //FIXME: if( openBoundaryPatches.size() > 0 ) {
      //FIXME:    for( long unsigned int j = 0; j < std::size( openBoundaryPatches ); j++ ) {
      //FIXME:       //FIXME: At this point, updateSolidBoundaryOpenBoundary doesn't use ghost zones.
      //FIXME:       //       Here, we should update boundary ghost zones similarly to fluid procedure.
      //FIXME:       model.updateSolidBoundaryOpenBoundary( boundarySets[ i ], openBoundaryPatches[ j ], modelParams );
      //FIXME:    }
      //FIXME: }
      model.finalizeBoundaryInteraction( fluidSets[ i ], boundarySets[ i ], modelParams );

      // updat fluid
      //std::cout << "Set i: " << i << " updating fluid." << std::endl;
      model.interaction( fluidSets[ i ], boundarySets[ i ], modelParams );
      //FIXME: if( openBoundaryPatches.size() > 0 ) {
      //FIXME:    for( long unsigned int j = 0; j < std::size( openBoundaryPatches ); j++ ) {
      //FIXME:       openBoundaryPatches[ j ]->zone.updateParticlesInZone( fluidSets[ i ]->getParticles() );
      //FIXME:       model.interactionWithOpenBoundary( fluidSets[ i ], openBoundaryPatches[ j ], modelParams );
      //FIXME:    }
      //FIXME: }
            if( i == 1 )
            model.interactionWithOpenBoundary( fluidSets[ i ], multiresolutionBoundaryPatches[ i ], modelParams );
      model.finalizeInteraction( fluidSets[ i ], boundarySets[ i ], modelParams );
   }

   timeMeasurement.stop( "interact" );
   writeLog( "Interact...", "Done." );
}

template< typename Model >
void
SolverMultiSet< Model >::computeTimeStep()
{
   //FIXME: Individual time stepping in subdomains
   //FIXME: timeStepping.computeTimeStep( fluid, modelParams );
}

template< typename Model >
void
SolverMultiSet< Model >::updateTime()
{
   timeStepping.updateTimeStep();
}

//FIXME: Need update for multiset scheme
template< typename Model >
void
SolverMultiSet< Model >::measure()
{
   for( int i = 0; i < numberOfSubsets; i++ )
      simulationMonitor.template measure< typename ModelParams::KernelFunction, typename ModelParams::EOS >(
            fluidSets[ i ], boundarySets[ i ], modelParams, timeStepping, logger, verbose );
}

template< typename Model >
template< typename Stage >
void
SolverMultiSet< Model >::integrate( const Stage integrationStage, const bool integrateBoundary )
{
   timeMeasurement.start( "integrate" );
   for( int i = 0; i < numberOfSubsets; i++ )
      integrator->integratStepVerlet( fluidSets[ i ], boundarySets[ i ], timeStepping, integrateBoundary );
   timeMeasurement.stop( "integrate" );
   writeLog( "Integrate...", "Done." );
}

#ifdef HAVE_MPI

template< typename Model >
void
SPHMultiset_CFD< Model >::synchronizeDistributedSimulation()
{
   writeLog( "Starting synchronization.", "" );
   timeMeasurement.start( "synchronize" );
   fluid->synchronizeObject();
   boundary->synchronizeObject();
   timeMeasurement.stop( "synchronize" );
   writeLog( "Synchronize...", "Done." );
}

template< typename Model >
void
SPHMultiset_CFD< Model >::resetOverlaps()
{
   //TODO: This should be paritcles method
   fluid->getParticles()->removeParitclesOutOfDomain();
   boundary->getParticles()->removeParitclesOutOfDomain();
}

template< typename Model >
void
SPHMultiset_CFD< Model >::performLoadBalancing()
{

   //synchronize comp. time
   fluid->getDistributedParticles()->setNumberOfParticlesForLoadBalancing( fluid->getNumberOfParticles() ); //TODO: Remove ptcs duplicity
   fluid->getDistributedParticles()->setCompTimeForLoadBalancing( timeMeasurement.getTotalTime() - subdomainCompTimeBackup );
   fluid->synchronizeBalancingMeasures();

   //compare computational time / number of particles
   std::pair< IndexVectorType, VectorType > subdomainAdjustment;
   if( loadBalancingMeasure == "computationalTime" )
      subdomainAdjustment = fluid->getDistributedParticles()->loadBalancingDomainAdjustmentCompTime();
   else if( loadBalancingMeasure == "numberOfParticles" )
      subdomainAdjustment = fluid->getDistributedParticles()->loadBalancingDomainAdjustment();
   else
      std::cerr << "Invalid load balancing metrics. Load balancing metrics is: " << loadBalancingMeasure << "." << std::endl;

   const IndexVectorType gridDimensionsAdjustment = subdomainAdjustment.first;
   const VectorType gridOriginAdjustment = subdomainAdjustment.second * fluid->getParticles()->getSearchRadius();

   const IndexVectorType updatedGridDimensions = fluid->getParticles()->getGridDimensions() + gridDimensionsAdjustment;
   const VectorType updatedGridOrigin = fluid->getParticles()->getGridOrigin() + gridOriginAdjustment;

   logger.writeParameter( "Load balancing - subdomain adjustment: ", "" );
   logger.writeParameter( "Grid dimensions adjustment: ", gridDimensionsAdjustment );
   logger.writeParameter( "Grid origin adjustment: ", gridOriginAdjustment );
   logger.writeParameter( "Old grid dimensions: ", fluid->getParticles()->getGridDimensions() );
   logger.writeParameter( "Old grid origin adjustment: ", fluid->getParticles()->getGridOrigin() );
   logger.writeParameter( "Old firstLastCellParticleList size: ", fluid->getParticles()->getCellFirstLastParticleList().getSize() );

   //update size of subdomain
   fluid->getParticles()->setGridDimensions( updatedGridDimensions );
   fluid->getParticles()->setGridOrigin( updatedGridOrigin );
   boundary->getParticles()->setGridDimensions( updatedGridDimensions );
   boundary->getParticles()->setGridOrigin( updatedGridOrigin );

   const IndexVectorType updatedGridOriginGlobalCoords = fluid->getParticles()->getGridOriginGlobalCoords() + subdomainAdjustment.second;
   fluid->getParticles()->setGridOriginGlobalCoords( updatedGridOriginGlobalCoords );
   boundary->getParticles()->setGridOriginGlobalCoords( updatedGridOriginGlobalCoords );

   logger.writeParameter( "New grid dimensions: ", fluid->getParticles()->getGridDimensions() );
   logger.writeParameter( "New grid origin adjustment: ", fluid->getParticles()->getGridOrigin() );
   logger.writeParameter( "New firstLastCellParticleList size: ", fluid->getParticles()->getCellFirstLastParticleList().getSize() );

   //update distributed particles and overlaps
   //TODO: 1 stands for overlapWidth, pass as parameter
   fluid->getDistributedParticles()->updateDistriutedGridParameters( updatedGridDimensions,
                                                                     updatedGridOrigin,
                                                                     1,
                                                                     fluid->getParticles()->getSearchRadius() );
   boundary->getDistributedParticles()->updateDistriutedGridParameters( updatedGridDimensions,
                                                                        updatedGridOrigin,
                                                                        1,
                                                                        boundary->getParticles()->getSearchRadius() );
   writeLoadBalancingInfo( gridDimensionsAdjustment[ 0 ] );
   this->subdomainCompTimeBackup = timeMeasurement.getTotalTime();

}

template< typename Model >
void
SPHMultiset_CFD< Model >::writeLoadBalancingInfo( const int gridResize )
{
   const std::string outputPath = outputDirectory + "/loadBalancingMetrics_rank" + std::to_string( TNL::MPI::GetRank() ) + ".dat";
   std::ofstream outfile;
   outfile.open(outputPath, std::ios_base::app );
   outfile << timeStepping.getStep() << " "
           << timeStepping.getTime() << " "
           << fluid->getNumberOfParticles() << " "
           << fluid->getNumberOfAllocatedParticles() << " "
           << boundary->getNumberOfParticles() << " "
           << boundary->getNumberOfAllocatedParticles() << " "
           << timeMeasurement.getTotalTime() - subdomainCompTimeBackup << " "
           << fluid->getParticles()->getCellFirstLastParticleList().getSize() << " "
           << gridResize << std::endl;
}

template< typename Model >
void
SPHMultiset_CFD< Model >::balanceSubdomains()
{
   if( ( timeStepping.getStep() % loadBalancingStepInterval  == 0 ) && ( timeStepping.getStep() > 1 ) ){

      performNeighborSearch( true );

      logger.writeSeparator();
      writeLog( "Starting load balancing.", "" );
      timeMeasurement.start( "rebalance" );
      performLoadBalancing();
      timeMeasurement.stop( "rebalance" );
      writeLog( "Load balancing...", "Done." );
      logger.writeSeparator();
      TNL::MPI::Barrier( communicator );

      // reset overlaps and synchronize with new domain syze
      resetOverlaps();
      writeLog( "Reset overlaps...", "Done." );
      TNL::MPI::Barrier( communicator );

      performNeighborSearch( true );
      TNL::MPI::Barrier( communicator );

      synchronizeDistributedSimulation();
      TNL::MPI::Barrier( communicator );
   }
}


#endif

template< typename Model >
void
SolverMultiSet< Model >::save( bool writeParticleCellIndex )
{
   if( ( verbose == "with-snapshot" ) || ( verbose == "full" ) )
      writeInfo();

   const int step = timeStepping.getStep();
   const RealType time = timeStepping.getTime();

   for( int i = 0; i < numberOfSubsets; i++ ){
#ifdef HAVE_MPI
      std::string outputFileNameFluid = outputDirectory + "/fluid_rank" + std::to_string( TNL::MPI::GetRank() ) + "_" + std::to_string( time ) + "_particles.vtk";
#else
      std::string outputFileNameFluid = outputDirectory + "/fluid_subdomain" + std::to_string( i ) + "_" + std::to_string( time ) + "_particles.vtk";
#endif
      fluidSets[ i ]->template writeParticlesAndVariables< Writer >( outputFileNameFluid, writeParticleCellIndex );
      logger.writeParameter( "Saved:", outputFileNameFluid );

#ifdef HAVE_MPI
      std::string outputFileNameBound = outputDirectory + "/boundary_rank" + std::to_string( TNL::MPI::GetRank() ) + "_" + std::to_string( time ) + "_particles.vtk";
#else
      std::string outputFileNameBound = outputDirectory + "/boundary_subdomain" + std::to_string( i ) + "_" +  std::to_string( time ) + "_particles.vtk";
#endif
      boundarySets[ i ]->template writeParticlesAndVariables< Writer >( outputFileNameBound, writeParticleCellIndex );
      logger.writeParameter( "Saved:", outputFileNameBound );

      //FIXME: if( openBoundaryPatches.size() ) {
      //FIXME:    for( auto& openBoundaryPatch : openBoundaryPatches ) {
      //FIXME:       std::string outputFileNameOpenBound =
      //FIXME:          outputDirectory + "/" + openBoundaryPatch->parameters.identifier + "_" + std::to_string( time ) + "_particles.vtk";
      //FIXME:       openBoundaryPatch->template writeParticlesAndVariables< Writer >( outputFileNameOpenBound, writeParticleCellIndex );
      //FIXME:       logger.writeParameter( "Saved:", outputFileNameOpenBound );
      //FIXME:    }
      //FIXME: }

      if( multiresolutionBoundaryPatches.size() ) {
         //for( auto& multiresolutionBoundaryPatch : multiresolutionBoundaryPatches ) {
            std::string outputFileNameMultiresolutionBound =
               outputDirectory + "/multiresolutionBoundaryPatch_subdomain" + std::to_string( i ) + "_" + std::to_string( time ) + "_particles.vtk";
            multiresolutionBoundaryPatches[ i ]->template writeParticlesAndVariables< Writer >(
                  outputFileNameMultiresolutionBound, writeParticleCellIndex );
            logger.writeParameter( "Saved:", outputFileNameMultiresolutionBound );
         //}
      }

#ifdef HAVE_MPI
      std::string outputFileNameGrid = outputDirectory + "/grid_rank" + std::to_string( TNL::MPI::GetRank() + 1 ) + "_" + std::to_string( time ) + ".vtk";
#else
      std::string outputFileNameGrid = outputDirectory + "/grid_subdomain" + std::to_string( i ) + "_" + std::to_string( time ) + ".vtk";
#endif
      TNL::Writers::writeBackgroundGrid( outputFileNameGrid, fluidSets[ i ]->getParticles()->getGridDimensions(), fluidSets[ i ]->getParticles()->getGridOrigin(), fluidSets[ i ]->getParticles()->getSearchRadius() );
      logger.writeParameter( "Saved:", outputFileNameGrid );

      // output simulation sensors to files
      simulationMonitor.save( logger );
   }
}

template< typename Model >
void
SolverMultiSet< Model >::makeSnapshot()
{
   const bool savePressure = true;
   if( timeStepping.checkOutputTimer( "save_results" ) ){
     if( savePressure ){
         for( int i = 0; i < numberOfSubsets; i++ ){
            model.computePressureFromDensity( fluidSets[ i ], modelParams );
            model.computePressureFromDensity( boundarySets[ i ], modelParams );
         }
     }
     save();
     writeLog( "Save results...", "Done." );
   }
}

template< typename Model >
void
SolverMultiSet< Model >::writeProlog( bool writeSystemInformation ) noexcept
{
   logger.writeHeader( "SPH simulation configuration." );
   logger.writeParameter( "Case name:", caseName );
   logger.writeSeparator();

   const bool printGPUs = std::is_same< DeviceType, TNL::Devices::Cuda >::value;

   if( TNL::MPI::isInitialized() )
      logger.writeParameter( "MPI processes:", TNL::MPI::GetSize() );
   logger.writeParameter( "Device type:", TNL::getType< DeviceType >() );
   if( ! printGPUs ) {
      if( TNL::Devices::Host::isOMPEnabled() ) {
         logger.writeParameter( "OMP enabled:", "yes", 1 );
         logger.writeParameter( "OMP threads:", TNL::Devices::Host::getMaxThreadsCount(), 1 );
      }
      else
         logger.writeParameter( "OMP enabled:", "no", 1 );
   }

   logger.writeParameter( "Particle system type:", Model::ParticlesType::writeModelType() );
   logger.writeParameter( "SPH model:", Model::writeModelType() );
   logger.writeParameter( "Verbose:", verbose );
   logger.writeParameter( "Output directory:", outputDirectory );
   logger.writeParameter( "Particles format", particlesFormat );
   writePrologModel( logger, modelParams );

   std::cout << "NUMBER OF SETS:" << numberOfSubsets << std::endl;
   for( int i = 0; i < numberOfSubsets; i++ ){
      std::cout << "Print objects!" << std::endl;
      logger.writeHeader( "Fluid " + std::to_string( i ) + " object information." );
      fluidSets[ i ]->writeProlog( logger );
      logger.writeHeader( "Boundary " + std::to_string( i ) + " object information:" );
      boundarySets[ i ]->writeProlog( logger );
      std::cout << "Print done!" << std::endl;
   }

   if( openBoundaryPatches.size() > 0 ) {
      for( long unsigned int i = 0; i < openBoundaryPatches.size(); i++ ) {
         logger.writeHeader( "Open boundary buffer" + std::to_string( i + 1 ) + "." );
         openBoundaryPatches[ i ]->writeProlog( logger );
         logger.writeSeparator();
         openBoundaryPatches[ i ]->config.writeProlog( logger );
      }
   }
   //FIXME: if constexpr( Model::ModelConfigType::SPHConfig::numberOfPeriodicBuffers > 0 ) {
   //FIXME:    const long unsigned int numberOfPeriodicPatches = std::size( fluid->periodicPatches );
   //FIXME:    for( long unsigned int i = 0; i < numberOfPeriodicPatches; i++ ) {
   //FIXME:       logger.writeHeader( "Periodic boundary patch " + std::to_string( i + 1 ) + "." );
   //FIXME:       fluid->periodicPatches[ i ]->writeProlog( logger );
   //FIXME:       //periodicity is the same for the boundary, so we don't need to print it
   //FIXME:       //boundary->periodicPatches[ i ]->writeProlog( logger );
   //FIXME:    }
   //FIXME: }

   if( multiresolutionBoundaryPatches.size() > 0 ) {
      for( long unsigned int i = 0; i < multiresolutionBoundaryPatches.size(); i++ ) {
         logger.writeHeader( "Multiresolutin boundary buffer" + std::to_string( i + 1 ) + "." );
         multiresolutionBoundaryPatches[ i ]->writeProlog( logger, i );
      }
   }

   logger.writeHeader( "System information." );
   if( writeSystemInformation ) {
      logger.writeSystemInformation( printGPUs );
      logger.writeSeparator();
      logger.writeCurrentTime( "Started at:" );
      logger.writeSeparator();
   }
}

template< typename Model >
template< typename ParameterType >
void
SolverMultiSet< Model >::writeLog( const std::string& label, const ParameterType& value, int parameterLevel )
{
   if( verbose == "full" )
      logger.writeParameter( label, value, parameterLevel );
}


template< typename Model >
void
SolverMultiSet< Model >::writeInfo() noexcept
{
   logger.writeSeparator();
   logger.writeParameter( "Simulation time: " + std::to_string( timeStepping.getTime() )
                             + " s, simulation step: " + std::to_string( timeStepping.getStep() ),
                          "" );
   logger.writeCurrentTime( "Current time:" );
   for( int i = 0; i < numberOfSubsets; i++ ){
      logger.writeParameter( "Number of fluid particles:", fluidSets[ i ]->getNumberOfParticles() );
      logger.writeParameter( "Number of allocated fluid particles:", fluidSets[ i ]->getNumberOfAllocatedParticles() );
      logger.writeParameter( "Number of boundary particles:", boundarySets[ i ]->getNumberOfParticles() );
      logger.writeParameter( "Number of allocated fluid particles:", boundarySets[ i ]->getNumberOfAllocatedParticles() );
   }
   if( openBoundaryPatches.size() > 0 ) {
      for( long unsigned int i = 0; i < openBoundaryPatches.size(); i++ )
         logger.writeParameter( "Number of buffer" + std::to_string( i + 1 ) + " particles:",
                                openBoundaryPatches[ i ]->getNumberOfParticles() );
   }
   //FIXME: if( openBoundaryPatches.size() > 0 ) {
   //FIXME:    if( verbose == "full" ) {
   //FIXME:       for( long unsigned int i = 0; i < fluid->periodicPatches.size(); i++ )
   //FIXME:          logger.writeParameter( "Number of fluid particles in periodic patch " + std::to_string( i + 1 ) + ": ",
   //FIXME:                                 fluid->periodicPatches[ i ]->particleZone.getNumberOfParticles() );
   //FIXME:       for( long unsigned int i = 0; i < boundary->periodicPatches.size(); i++ )
   //FIXME:          logger.writeParameter( "Number of boundary particles in periodic patch " + std::to_string( i + 1 ) + " :",
   //FIXME:                                 boundary->periodicPatches[ i ]->particleZone.getNumberOfParticles() );
   //FIXME:    }
   //FIXME: }


   if( multiresolutionBoundaryPatches.size() > 0 ) {
      for( long unsigned int i = 0; i < multiresolutionBoundaryPatches.size(); i++ )
         logger.writeParameter( "Number of mr-buffer " + std::to_string( i ) + " particles:",
               multiresolutionBoundaryPatches[ i ]->getNumberOfParticles() );
   }
   logger.writeSeparator();
}

template< typename Model >
void
SolverMultiSet< Model >::writeEpilog() noexcept
{
   logger.writeHeader( "SPH simulation successfully finished." );
   logger.writeCurrentTime( "Ended at:" );
   timeMeasurement.writeInfo( logger, timeStepping.getStep() );

   std::string saveTimersOutputName = outputDirectory + "/time_measurements";
   timeMeasurement.writeInfoToJson( saveTimersOutputName, timeStepping.getStep() );
}

}  //namespace SPH
}  //namespace TNL
