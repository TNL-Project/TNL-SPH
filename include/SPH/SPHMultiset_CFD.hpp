#include "SPHMultiset_CFD.h"
#include "SPH/OpenBoundaryConfig.h"
#include "SPH/TimeMeasurement.h"
#include <iterator>
#include <ostream>
#include <string>

#include "distributedUtils.h"
#include "../Writers/writeBackgroundGrid.h"

namespace TNL {
namespace SPH {

template< typename Model >
void
SPHMultiset_CFD< Model >::init( TNL::Config::ParameterContainer& parameters, TNL::Logger& logger )
{
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
   initOverlaps( parameters, this->parametersDistributed, logger );
#else
   // initialize particle sets
   initParticleSets( parameters, logger );
#endif

   // init periodic boundary conditions
   //TODO: I don't like that open boundary buffer are selected in compile time.
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfPeriodicBuffers > 0 ) {
      fluid->initializePeriodicity( parameters );
      boundary->initializePeriodicity( parameters );

      timeMeasurement.addTimer( "periodicity-fluid-updateZone", false );
      timeMeasurement.addTimer( "periodicity-boundary-updateZone", false );
   }

   // init model parameters
   modelParams.init( parameters );

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
   readParticlesFiles( parameters, logger );
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
SPHMultiset_CFD< Model >::initParticleSets( TNL::Config::ParameterContainer& parameters, TNL::Logger& logger )
{
   // get global domain properetis
   const VectorType domainOrigin = parameters.getXyz< VectorType >( "domainOrigin" );
   const VectorType domainSize = parameters.getXyz< VectorType >( "domainSize" );
   const RealType searchRadius = parameters.getParameter< RealType >( "searchRadius" );
   const IndexVectorType gridSize = TNL::ceil( ( domainSize - domainOrigin ) / searchRadius );

   // init fluid
   fluid->initialize( parameters.getParameter< int >( "numberOfParticles" ),
                      parameters.getParameter< int >( "numberOfAllocatedParticles" ),
                      searchRadius,
                      gridSize,
                      domainOrigin );
   if constexpr( ParticlesType::specifySearchedSetExplicitly() == true )
      fluid->particles->setParticleSetLabel( 0 );

   // init boundary
   boundary->initialize( parameters.getParameter< int >( "numberOfBoundaryParticles" ),
                         parameters.getParameter< int >( "numberOfAllocatedBoundaryParticles" ),
                         searchRadius,
                         gridSize,
                         domainOrigin );
   if constexpr( ParticlesType::specifySearchedSetExplicitly() == true )
      boundary->particles->setParticleSetLabel( 1 );

   // init open boundary patches
   const int numberOfBoundaryPatches = parameters.getParameter< int >( "openBoundaryPatches" );
   //TODO: I don't like that open boundary buffer are selected in compile time.
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ) {
      openBoundaryPatches.resize( numberOfBoundaryPatches );
      for( int i = 0; i < numberOfBoundaryPatches; i++ ) {
         std::string prefix = "buffer-" + std::to_string( i + 1 ) + "-";
         openBoundaryPatches[ i ]->config.init( parameters, prefix );
         openBoundaryPatches[ i ]->initialize( parameters.getParameter< int >( prefix + "numberOfParticles" ),
                                               parameters.getParameter< int >( prefix + "numberOfAllocatedParticles" ),
                                               searchRadius,
                                               gridSize,
                                               domainOrigin );
      }
   }
}

#ifdef HAVE_MPI
template< typename Model >
void
SPHMultiset_CFD< Model >::initDistributedParticleSets( TNL::Config::ParameterContainer& parameters,
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
   const GlobalIndexType numberOfOverlapLayers = parameters.getParameter< int >( "overlapWidth" );

   // subdomain + ghost properties
   const VectorType subdomainOrigin = parametersDistributed.getXyz< VectorType >( subdomainKey + "origin" );
   //const VectorType subdomainSize = parametersDistributed.getXyz< VectorType >( subdomainKey + "size" ) ;
   //const IndexVectorType subdomainGridDimension = TNL::ceil( subdomainSize / searchRadius );
   const IndexVectorType subdomainGridDimension = parametersDistributed.getXyz< IndexVectorType >( subdomainKey + "grid-dimensions" );
   const IndexVectorType subdomainGridOriginGlobalCoords = parametersDistributed.getXyz< IndexVectorType >( subdomainKey + "origin-global-coords" );

   // init fluid
   logger.writeParameter( "initDistributed:", "fluid->initialize" );
   fluid->initializeAsDistributed( parametersDistributed.getParameter< int >( subdomainKey + "fluid_n" ),
                                   parametersDistributed.getParameter< int >( subdomainKey + "fluid_n_allocated" ),
                                   searchRadius,
                                   subdomainGridDimension,
                                   subdomainOrigin,
                                   subdomainGridOriginGlobalCoords,
                                   domainOrigin,
                                   logger );
   //fluid->particles->interiorSize = subdomainSize; //FIXME Getter, Setter

   // init boundary
   logger.writeParameter( "initDistributed:", "boundary->initialize" );
   boundary->initializeAsDistributed( parametersDistributed.getParameter< int >( subdomainKey + "boundary_n" ),
                                      parametersDistributed.getParameter< int >( subdomainKey + "boundary_n_allocated" ),
                                      searchRadius,
                                      subdomainGridDimension,
                                      subdomainOrigin,
                                      subdomainGridOriginGlobalCoords,
                                      domainOrigin,
                                      logger );
   //boundary->particles->interiorSize = subdomainSize; //FIXME Getter, Setter

   // set distributed particle system: FIXME: All this lines are ugly
  fluid->distributedParticles->setDistributedGridParameters( domainGridDimension,
                                                             domainOrigin,
                                                             subdomainGridDimension,
                                                             subdomainOrigin,
                                                             numberOfOverlapLayers,
                                                             searchRadius,
                                                             numberOfSubdomains,
                                                             this->communicator );
  fluid->distributedParticles->writeProlog( logger );
  fluid->synchronizer.initialize( fluid->distributedParticles );
  fluid->synchronizer.setCommunicator( this->communicator );

  boundary->distributedParticles->setDistributedGridParameters( domainGridDimension,
                                                                domainOrigin,
                                                                subdomainGridDimension,
                                                                subdomainOrigin,
                                                                numberOfOverlapLayers,
                                                                searchRadius,
                                                                numberOfSubdomains,
                                                                this->communicator );
  boundary->distributedParticles->writeProlog( logger );
  boundary->synchronizer.initialize( boundary->distributedParticles );
  boundary->synchronizer.setCommunicator( this->communicator );
}
#endif

#ifdef HAVE_MPI
template< typename Model >
void
SPHMultiset_CFD< Model >::initOverlaps( TNL::Config::ParameterContainer& parameters,
                                        TNL::Config::ParameterContainer& parametersDistributed,
                                        TNL::Logger& logger )
{

   //TODO: This whole header can be hidden to distributed utils
   int rank = TNL::MPI::GetRank();
   Containers::StaticVector< 2, int > numberOfSubdomains = parameters.getXyz< Containers::StaticVector< 2, int > >( "subdomains" );
   const std::string subdomainKey = distributed::getSubdomainKey( TNL::MPI::GetRank(), numberOfSubdomains );

   // global domain properties
   const RealType searchRadius = parameters.getParameter< RealType >( "searchRadius" );
   const VectorType domainOrigin = parameters.getXyz< VectorType >( "domainOrigin" );
   //const VectorType domainSize = parameters.getXyz< VectorType >( "domainSize" );
   //const IndexVectorType domainGridDimension = TNL::ceil( ( domainSize - domainOrigin ) / searchRadius );

   // subdomain properties
   const VectorType subdomainOrigin = parametersDistributed.getXyz< VectorType >( subdomainKey + "origin" );
   //const VectorType subdomainSize = parametersDistributed.getXyz< VectorType >(  subdomainKey + "size" );
   //const IndexVectorType subdomainGridSize = TNL::ceil( subdomainSize / searchRadius );
   const IndexVectorType subdomainGridSize = parametersDistributed.getXyz< IndexVectorType >( subdomainKey + "grid-dimensions" );
   const IndexVectorType subdomainGridOriginGlobalCoords = parametersDistributed.getXyz< IndexVectorType >( subdomainKey + "origin-global-coords" );

   int overlapCellsCount = 0;
   IndexVectorType resizedSubdomainGridSize = 0;
   VectorType resizedSubdomainGridOrigin = 0.f;

   //TODO: Consider whether the overlap is in all dimensions
   if constexpr( Model::SPHConfig::spaceDimension == 2 ) {
      overlapCellsCount = ( subdomainGridSize[ 0 ] + subdomainGridSize[ 1 ] + 4 ) * 2;

      //TODO: Maybe modify directly
      resizedSubdomainGridSize = { subdomainGridSize[ 0 ] + 2, subdomainGridSize[ 1 ] + 2 };
      resizedSubdomainGridOrigin = { subdomainOrigin[ 0 ] - searchRadius,  subdomainOrigin[ 1 ] - searchRadius };
   }
   else if constexpr( Model::SPHConfig::spaceDimension == 3 ) {
      const int xy = ( subdomainGridSize[ 0 ] + 2 ) * ( subdomainGridSize[ 1 ] + 2 );
      const int xz = ( subdomainGridSize[ 0 ] + 2 ) * ( subdomainGridSize[ 2 ] );
      const int yz = ( subdomainGridSize[ 1 ] ) * ( subdomainGridSize[ 2 ] );
      overlapCellsCount = 2 * xy + 2 * xz + 2 * yz;

      resizedSubdomainGridSize = { subdomainGridSize[ 0 ] + 2, subdomainGridSize[ 1 ] + 2, subdomainGridSize[ 2 ] + 2 };
      resizedSubdomainGridOrigin = { subdomainOrigin[ 0 ] - searchRadius,
                                     subdomainOrigin[ 1 ] - searchRadius,
                                     subdomainOrigin[ 2 ] - searchRadius };
   }
   int numberOfParticlesPerCell = parameters.getParameter< int >( "numberOfParticlesPerCell" );

   // FIXME: Here, the arguments are probably fcked
   fluidOverlap->initializeAsDistributed( 0,
                                          numberOfParticlesPerCell * overlapCellsCount,
                                          searchRadius,
                                          resizedSubdomainGridSize,
                                          resizedSubdomainGridOrigin,
                                          subdomainGridOriginGlobalCoords,
                                          domainOrigin,
                                          logger );

   // FIXME: Here, the arguments are probably fcked
   boundaryOverlap->initializeAsDistributed( 0,
                                             numberOfParticlesPerCell * overlapCellsCount,
                                             searchRadius,
                                             resizedSubdomainGridSize,
                                             resizedSubdomainGridOrigin,
                                             subdomainGridOriginGlobalCoords,
                                             domainOrigin,
                                             logger );
}
#endif

template< typename Model >
void
SPHMultiset_CFD< Model >::readParticlesFiles( TNL::Config::ParameterContainer& parameters, TNL::Logger& logger )
{
   if( parameters.getParameter< int >( "numberOfParticles" ) != 0 ){
      logger.writeParameter( "Reading fluid particles:", parameters.getParameter< std::string >( "fluid-particles" ) );
      fluid->template readParticlesAndVariables< SimulationReaderType >(
         parameters.getParameter< std::string >( "fluid-particles" ) );
   }
   logger.writeParameter( "Reading boundary particles:", parameters.getParameter< std::string >( "boundary-particles" ) );
   boundary->template readParticlesAndVariables< SimulationReaderType >(
      parameters.getParameter< std::string >( "boundary-particles" ) );

   // init open boundary patches
   const int numberOfBoundaryPatches = parameters.getParameter< int >( "openBoundaryPatches" );
   //TODO: I don't like that open boundary buffer are selected in compile time.
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ) {
      for( int i = 0; i < numberOfBoundaryPatches; i++ ) {
         std::string prefix = "buffer-" + std::to_string( i + 1 ) + "-";
         logger.writeParameter( "Reading open boundary particles:", parameters.getParameter< std::string >( prefix + "particles" ) );
         openBoundaryPatches[ i ]->template readParticlesAndVariables< SimulationReaderType >(
            parameters.getParameter< std::string >( prefix + "particles" ) );
      }
   }
}

#ifdef HAVE_MPI
template< typename Model >
void
SPHMultiset_CFD< Model >::readParticleFilesDistributed( TNL::Config::ParameterContainer& parameters,
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
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ) {  //TODO: I dont like this.
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
SPHMultiset_CFD< Model >::performNeighborSearch( TNL::Logger& logger, bool performBoundarySearch )
{
   if constexpr( ParticlesType::specifySearchedSetExplicitly() == false ){
      fluid->searchForNeighbors();
      if( verbose == "full" )
         logger.writeParameter( "Fluid search procedure:", "Done." );

      if( timeStepping.getStep() == 0 || performBoundarySearch == true ){
         boundary->searchForNeighbors();
         if( verbose == "full" )
            logger.writeParameter( "Boundary search procedure:", "Done." );
      }

      if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 )
         for( auto& openBoundaryPatch : openBoundaryPatches ){
            openBoundaryPatch->searchForNeighbors();
            if( verbose == "full" )
               logger.writeParameter( "Open boundary patch search procedure:", "Done." );
         }
   }
   else if constexpr( ParticlesType::specifySearchedSetExplicitly() == true ){

      fluid->makeSetSearchable();
      if( verbose == "full" )
         logger.writeParameter( "Fluid-fluid search procedure:", "Done." );

      if( timeStepping.getStep() == 0 || performBoundarySearch == true ){
         boundary->makeSetSearchable();
         if( verbose == "full" )
            logger.writeParameter( "Boundary-boundary search procedure:", "Done." );
      }

               // D:
               //const int numberOfPtcs = fluid->getNumberOfParticles();
               //const int offsetParticle = 0;
               //for( int j = 0; j < 50; j++ ){
               //   std::cout << fluid->particles->getNeighborListStorage().getElement( j * numberOfPtcs + offsetParticle ) << " ";

               //}
               //std::cout << std::endl;

      fluid->particles->addToParticleList( fluid->getParticles() );
      if( verbose == "full" )
         logger.writeParameter( "Fluid-boundary search procedure:", "Done." );

      fluid->particles->addToParticleList( boundary->getParticles() );
      if( verbose == "full" )
         logger.writeParameter( "Fluid-boundary search procedure:", "Done." );

      boundary->particles->addToParticleList( fluid->getParticles() );
      if( verbose == "full" )
         logger.writeParameter( "Boundary-fluid search procedure:", "Done." );

               // D:
               //for( int j = 0; j < 50; j++ ){
               //   std::cout << fluid->particles->getNeighborListStorage().getElement( j * numberOfPtcs + offsetParticle ) << " ";
               //}
               //std::cout << std::endl;

               // D:
               //const int numberOfPtcsBoundary = boundary->getNumberOfParticles();
               //const int offsetParticleBoundary = 0;
               //for( int j = 0; j < 50; j++ ){
               //   std::cout << boundary->particles->getNeighborListStorage().getElement( j * numberOfPtcsBoundary + offsetParticleBoundary ) << " ";
               //}
               //std::cout << std::endl;
   }
}

template< typename Model >
template< typename ParticleSetPointer >
void
SPHMultiset_CFD< Model >::performNeighborSearchForObject( ParticleSetPointer& objectPointer )
{
   objectPointer->particles->resetListWithIndices();
   objectPointer->particles->computeParticleCellIndices();
   objectPointer->sortParticles();
   objectPointer->particles->particlesToCells();
}

template< typename Model >
void
SPHMultiset_CFD< Model >::extrapolateOpenBC()
{
   for( long unsigned int i = 0; i < std::size( openBoundaryPatches ); i++ ) {
      //TODO Check if open boundary buffer is really open boundary buffer
      model.extrapolateOpenBoundaryData( fluid, openBoundaryPatches[ i ], modelParams, openBoundaryPatches[ i ]->config );
   }
}

template< typename Model >
void
SPHMultiset_CFD< Model >::applyOpenBC( const RealType timeStepFact )
{
   for( long unsigned int i = 0; i < std::size( openBoundaryPatches ); i++ ) {
      //TODO Check if open boundary buffer is really open boundary buffer
      openBoundaryModel.applyOpenBoundary(
         timeStepFact * timeStepping.getTimeStep(), fluid, openBoundaryPatches[ i ], openBoundaryPatches[ i ]->config );
   }
}

template< typename Model >
void
SPHMultiset_CFD< Model >::applyPeriodicBCEnforce()
{
   timeMeasurement.start( "periodicity-fluid-updateZone" );
   fluid->enforcePeriodicPatches();
   timeMeasurement.stop( "periodicity-boundary-updateZone" );
   timeMeasurement.start( "periodicity-fluid-updateZone" );
   boundary->enforcePeriodicPatches();
   timeMeasurement.stop( "periodicity-boundary-updateZone" );
}

template< typename Model >
void
SPHMultiset_CFD< Model >::applyPeriodicBCTransfer()
{
   const long unsigned int numberOfPeriodicPatches = std::size( fluid->periodicPatches );
   for( long unsigned int i = 0; i < numberOfPeriodicPatches; i++ ) {
      openBoundaryModel.periodicityParticleTransfer( fluid, fluid->periodicPatches[ i ] );
   }
}

template< typename Model >
void
SPHMultiset_CFD< Model >::interact()
{
   // update solid boundary conditions
   model.updateSolidBoundary( fluid, boundary, modelParams );
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ) {
      for( long unsigned int i = 0; i < std::size( openBoundaryPatches ); i++ ) {
         //FIXME: At this point, updateSolidBoundaryOpenBoundary doesn't use ghost zones.
         //       Here, we should update boundary ghost zones similarly to fluid procedure.
         model.updateSolidBoundaryOpenBoundary( boundary, openBoundaryPatches[ i ], modelParams );
      }
   }
   model.finalizeBoundaryInteraction( fluid, boundary, modelParams );

   // updat fluid
   model.interaction( fluid, boundary, modelParams );
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ) {
      for( long unsigned int i = 0; i < std::size( openBoundaryPatches ); i++ ) {
         openBoundaryPatches[ i ]->zone.updateParticlesInZone( fluid->particles );
         model.interactionWithOpenBoundary( fluid, openBoundaryPatches[ i ], modelParams );
      }
   }
   model.finalizeInteraction( fluid, boundary, modelParams );
}

template< typename Model >
void
SPHMultiset_CFD< Model >::updateTimeStep()
{
   timeStepping.computeTimeStep( fluid, modelParams );
   timeStepping.updateTimeStep();
}

template< typename Model >
template< typename SPHKernelFunction, typename EOS >
void
SPHMultiset_CFD< Model >::measure( TNL::Logger& logger )
{
   simulationMonitor.template measure< SPHKernelFunction, EOS >( fluid, boundary, modelParams, timeStepping, logger, verbose );
}

#ifdef HAVE_MPI

template< typename Model >
void
SPHMultiset_CFD< Model >::synchronizeDistributedSimulation( TNL::Logger& logger )
{
   fluid->synchronizeObject( fluidOverlap, logger );
   boundary->synchronizeObject( boundaryOverlap, logger );
}

template< typename Model >
void
SPHMultiset_CFD< Model >::resetOverlaps()
{
   //TODO: This should be paritcles method
   fluid->particles->removeParitclesOutOfDomain();
   boundary->particles->removeParitclesOutOfDomain();
}

template< typename Model >
void
SPHMultiset_CFD< Model >::performLoadBalancing( TNL::Logger& logger )
{
   //setup tresholds:
   fluid->distributedParticles->setParticlesCountResizeTrashold( 1000 );
   //fluid->distributedParticles->setCompTimeResizePercetnageTrashold( 0.05 );

   //synchronize comp. time
   fluid->distributedParticles->setNumberOfParticlesForLoadBalancing( fluid->getNumberOfParticles() );
   fluid->distributedParticles->setCompTimeForLoadBalancing( 0. ); //FIXME
   fluid->synchronizeBalancingMeasures();

   //compare computational time / number of particles
   std::pair< IndexVectorType, VectorType > subdomainAdjustment = fluid->distributedParticles->loadBalancingDomainAdjustment();
   const IndexVectorType gridDimensionsAdjustment = subdomainAdjustment.first;
   const VectorType gridOriginAdjustment = subdomainAdjustment.second * fluid->particles->getSearchRadius();

   const IndexVectorType updatedGridDimensions = fluid->particles->getGridDimensions() + gridDimensionsAdjustment;
   const VectorType updatedGridOrigin = fluid->particles->getGridOrigin() + gridOriginAdjustment;

   logger.writeParameter( "Load balancing - subdomain adjustment: ", "" );
   logger.writeParameter( "Grid dimensions adjustment: ", gridDimensionsAdjustment );
   logger.writeParameter( "Grid origin adjustment: ", gridOriginAdjustment );
   logger.writeParameter( "Old grid dimensions: ", fluid->particles->getGridDimensions() );
   logger.writeParameter( "Old grid origin adjustment: ", fluid->particles->getGridOrigin() );
   logger.writeParameter( "Old firstLastCellParticleList size: ", fluid->particles->getCellFirstLastParticleList().getSize() );

   //update size of subdomain
   fluid->particles->setGridDimensions( updatedGridDimensions );
   fluid->particles->setGridOrigin( updatedGridOrigin );
   boundary->particles->setGridDimensions( updatedGridDimensions );
   boundary->particles->setGridOrigin( updatedGridOrigin );

   const IndexVectorType updatedGridOriginGlobalCoords = fluid->particles->getGridOriginGlobalCoords() + subdomainAdjustment.second;
   fluid->particles->setGridOriginGlobalCoords( updatedGridOriginGlobalCoords );
   boundary->particles->setGridOriginGlobalCoords( updatedGridOriginGlobalCoords );

   logger.writeParameter( "New grid dimensions: ", fluid->particles->getGridDimensions() );
   logger.writeParameter( "New grid origin adjustment: ", fluid->particles->getGridOrigin() );
   logger.writeParameter( "New firstLastCellParticleList size: ", fluid->particles->getCellFirstLastParticleList().getSize() );

   //update distributed particles and overlaps
   //TODO: 1 stands for overlapWidth, pass as parameter
   fluid->distributedParticles->updateDistriutedGridParameters( updatedGridDimensions, updatedGridOrigin, 1, fluid->particles->getSearchRadius() );
   boundary->distributedParticles->updateDistriutedGridParameters( updatedGridDimensions, updatedGridOrigin, 1, boundary->particles->getSearchRadius() );

}

#endif

template< typename Model >
void
SPHMultiset_CFD< Model >::save( TNL::Logger& logger, bool writeParticleCellIndex )
{
   if( ( verbose == "with-snapshot" ) || ( verbose == "full" ) )
      writeInfo( logger );

   const int step = timeStepping.getStep();
#ifdef HAVE_MPI
   std::string outputFileNameFluid = outputDirectory + "/particles_rank" + std::to_string( TNL::MPI::GetRank() ) + "_" + std::to_string( step ) + "_fluid.vtk";
#else
   std::string outputFileNameFluid = outputDirectory + "/particles" + std::to_string( step ) + "_fluid.vtk";
#endif
   fluid->template writeParticlesAndVariables< Writer >( outputFileNameFluid, writeParticleCellIndex );
   logger.writeParameter( "Saved:", outputFileNameFluid );

#ifdef HAVE_MPI
   std::string outputFileNameBound = outputDirectory + "/particles_rank" + std::to_string( TNL::MPI::GetRank() ) + "_" + std::to_string( step ) + "_boundary.vtk";
#else
   std::string outputFileNameBound = outputDirectory + "/particles" + std::to_string( step ) + "_boundary.vtk";
#endif
   boundary->template writeParticlesAndVariables< Writer >( outputFileNameBound, writeParticleCellIndex );
   logger.writeParameter( "Saved:", outputFileNameBound );

   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ) {
      for( auto& openBoundaryPatch : openBoundaryPatches ) {
         std::string outputFileNameOpenBound =
            outputDirectory + "/particles" + std::to_string( step ) + "_" + openBoundaryPatch->parameters.identifier + ".vtk";
         openBoundaryPatch->template writeParticlesAndVariables< Writer >( outputFileNameOpenBound, writeParticleCellIndex );
         logger.writeParameter( "Saved:", outputFileNameOpenBound );
      }
   }

#ifdef HAVE_MPI
   std::string outputFileNameGrid = outputDirectory + "/grid_rank" + std::to_string( TNL::MPI::GetRank() + 1 ) + "_" + std::to_string( step ) + ".vtk";
#else
   std::string outputFileNameGrid = outputDirectory + "/grid" + std::to_string( step ) + ".vtk";
#endif
   TNL::Writers::writeBackgroundGrid( outputFileNameGrid, fluid->particles->getGridDimensions(), fluid->particles->getGridOrigin(), fluid->particles->getSearchRadius() );
   logger.writeParameter( "Saved:", outputFileNameGrid );

   // output simulation sensors to files
   simulationMonitor.save( logger );
}

template< typename Model >
void
SPHMultiset_CFD< Model >::makeSnapshot( TNL::Logger& logger )
{
   const bool savePressure = true;
   if( timeStepping.checkOutputTimer( "save_results" ) ){
     if( savePressure ){
        model.computePressureFromDensity( fluid, modelParams );
        model.computePressureFromDensity( boundary, modelParams );
     }
      save( logger );
      writeLog( logger, "Save results...", "Done." );
   }
}

template< typename Model >
void
SPHMultiset_CFD< Model >::writeProlog( TNL::Logger& logger, bool writeSystemInformation ) const noexcept
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
   logger.writeHeader( "Fluid object information." );
   fluid->writeProlog( logger );
   logger.writeHeader( "Boundary object information:" );
   boundary->writeProlog( logger );
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ) {
      for( long unsigned int i = 0; i < openBoundaryPatches.size(); i++ ) {
         logger.writeHeader( "Open boundary buffer" + std::to_string( i + 1 ) + "." );
         openBoundaryPatches[ i ]->writeProlog( logger );
         logger.writeSeparator();
         openBoundaryPatches[ i ]->config.writeProlog( logger );
      }
   }
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfPeriodicBuffers > 0 ) {
      const long unsigned int numberOfPeriodicPatches = std::size( fluid->periodicPatches );
      for( long unsigned int i = 0; i < numberOfPeriodicPatches; i++ ) {
         logger.writeHeader( "Periodic boundary patch " + std::to_string( i + 1 ) + "." );
         fluid->periodicPatches[ i ]->writeProlog( logger );
         //periodicity is the same for the boundary, so we don't need to print it
         //boundary->periodicPatches[ i ]->writeProlog( logger );
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
SPHMultiset_CFD< Model >::writeLog( TNL::Logger& logger,
                                    const std::string& label,
                                    const ParameterType& value,
                                    int parameterLevel )
{
   if( verbose == "full" )
      logger.writeParameter( label, value, parameterLevel );
}

template< typename Model >
void
SPHMultiset_CFD< Model >::writeInfo( TNL::Logger& logger ) const noexcept
{
   logger.writeSeparator();
   logger.writeParameter( "Simulation time: " + std::to_string( timeStepping.getTime() )
                             + " s, simulation step: " + std::to_string( timeStepping.getStep() ),
                          "" );
   logger.writeCurrentTime( "Current time:" );
   logger.writeParameter( "Number of fluid particles:", fluid->getNumberOfParticles() );
   logger.writeParameter( "Number of allocated fluid particles:", fluid->getNumberOfAllocatedParticles() );
   logger.writeParameter( "Number of boundary particles:", boundary->getNumberOfParticles() );
   logger.writeParameter( "Number of allocated fluid particles:", boundary->getNumberOfAllocatedParticles() );
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ) {
      for( long unsigned int i = 0; i < openBoundaryPatches.size(); i++ )
         logger.writeParameter( "Number of buffer" + std::to_string( i + 1 ) + " particles:",
                                openBoundaryPatches[ i ]->getNumberOfParticles() );
   }
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfPeriodicBuffers > 0 ) {
      if( verbose == "full" ) {
         for( long unsigned int i = 0; i < fluid->periodicPatches.size(); i++ )
            logger.writeParameter( "Number of fluid particles in periodic patch " + std::to_string( i + 1 ) + ": ",
                                   fluid->periodicPatches[ i ]->particleZone.getNumberOfParticles() );
         for( long unsigned int i = 0; i < boundary->periodicPatches.size(); i++ )
            logger.writeParameter( "Number of boundary particles in periodic patch " + std::to_string( i + 1 ) + " :",
                                   boundary->periodicPatches[ i ]->particleZone.getNumberOfParticles() );
      }
   }
   logger.writeSeparator();
}

template< typename Model >
void
SPHMultiset_CFD< Model >::writeEpilog( TNL::Logger& logger ) noexcept
{
   logger.writeHeader( "SPH simulation successfully finished." );
   logger.writeCurrentTime( "Ended at:" );
   timeMeasurement.writeInfo( logger, timeStepping.getStep() );

   std::string saveTimersOutputName = outputDirectory + "/time_measurements";
   timeMeasurement.writeInfoToJson( saveTimersOutputName, timeStepping.getStep() );
}

}  //namespace SPH
}  //namespace TNL
