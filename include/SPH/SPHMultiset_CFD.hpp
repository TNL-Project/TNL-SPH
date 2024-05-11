#include "SPHMultiset_CFD.h"
#include "SPH/OpenBoundaryConfig.h"
#include "SPH/TimeMeasurement.h"
#include <string>

#include "distributedUtils.h"

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

   // initialize particle sets and overlaps
   initDistributed( parameters, this->parametersDistributed, logger );
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
   outputDirecotry = parameters.getParameter< std::string >( "output-directory" );
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

   // init boundary
   boundary->initialize( parameters.getParameter< int >( "numberOfBoundaryParticles" ),
                         parameters.getParameter< int >( "numberOfAllocatedBoundaryParticles" ),
                         searchRadius,
                         gridSize,
                         domainOrigin );

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

template< typename Model >
void
SPHMultiset_CFD< Model >::initDistributed( TNL::Config::ParameterContainer& parameters,
                                           TNL::Config::ParameterContainer& parametersDistributed,
                                           TNL::Logger& logger )
{
   logger.writeHeader( "SPH simulation initialization." );

   int rank = TNL::MPI::GetRank();
   Containers::StaticVector< 2, int > numberOfSubdomains = parameters.getXyz< Containers::StaticVector< 2, int > >( "subdomains" );
   Containers::StaticVector< 2, int > subdomainCoordinates = distributed::restoreSubdomainCoordinatesFromRank( rank, numberOfSubdomains );
   const std::string subdomainKey = distributed::getSubdomainKey( rank, numberOfSubdomains );

   // global domain properties
   const RealType searchRadius = parameters.getParameter< RealType >( "searchRadius" );
   const VectorType domainOrigin = parameters.getXyz< VectorType >( "domainOrigin" );
   const VectorType domainSize = parameters.getXyz< VectorType >( "domainSize" );
   const IndexVectorType domainGridDimension = TNL::ceil( ( domainSize - domainOrigin ) / searchRadius );

   // subdomain + ghost properties
   const VectorType subdomainOrigin = parametersDistributed.getXyz< VectorType >( subdomainKey + "origin" );
   const VectorType subdomainSize = parametersDistributed.getXyz< VectorType >( subdomainKey + "size" ) ;
   const IndexVectorType subdomainGridDimension = TNL::ceil( subdomainSize / searchRadius );

   //debug, ugly with MPI, sync somehow
   logger.writeParameter( "Initializing rank: ", rank );
   logger.writeParameter( "Initializing subdomain: ", subdomainCoordinates );
   logger.writeParameter( "Subdomain key: ", subdomainKey );
   logger.writeParameter( "Initializing subdomain origin:", subdomainOrigin );
   logger.writeParameter( "Initializing subdomain size:", subdomainSize );
   logger.writeParameter( "Initializing subdomain size multiplied:", searchRadius * subdomainGridDimension );
   logger.writeParameter( "Initializing subdomain grid size:", subdomainGridDimension );

   // init fluid
   logger.writeParameter( "initDistributed:", "fluid->initialize" );
   fluid->initialize( parametersDistributed.getParameter< int >( subdomainKey + "fluid_n" ),
                      parametersDistributed.getParameter< int >( subdomainKey + "fluid_n_allocated" ),
                      searchRadius,
                      subdomainGridDimension,
                      subdomainOrigin,
                      domainGridDimension,
                      domainOrigin,
                      logger );
   fluid->particles->interiorSize = subdomainSize; //FIXME Getter, Setter
   logger.writeParameter( "subdomainSize:", subdomainSize );
   logger.writeParameter( "subdomainSize - multiplied:", searchRadius * subdomainGridDimension );
   logger.writeSeparator();

   // init boundary
   logger.writeParameter( "initDistributed:", "boundary->initialize" );
   boundary->initialize( parametersDistributed.getParameter< int >( subdomainKey + "boundary_n" ),
                         parametersDistributed.getParameter< int >( subdomainKey + "boundary_n_allocated" ),
                         searchRadius,
                         subdomainGridDimension,
                         subdomainOrigin,
                         domainGridDimension,
                         domainOrigin,
                         logger );
   boundary->particles->interiorSize = subdomainSize; //FIXME Getter, Setter

   // set distributed particle system: FIXME: All this lines are ugly
  fluid->distributedParticles->setDistributedGridParameters( domainGridDimension,
                                                             domainOrigin,
                                                             subdomainGridDimension,
                                                             subdomainOrigin,
                                                             searchRadius,
                                                             numberOfSubdomains,
                                                             this->communicator,
                                                             logger );
  fluid->distributedParticles->writeProlog( logger );
  fluid->synchronizer.initialize( fluid->distributedParticles );
  fluid->synchronizer.setCommunicator( this->communicator );

  boundary->distributedParticles->setDistributedGridParameters( domainGridDimension,
                                                                domainOrigin,
                                                                subdomainGridDimension,
                                                                subdomainOrigin,
                                                                searchRadius,
                                                                numberOfSubdomains,
                                                                this->communicator,
                                                                logger );
  boundary->distributedParticles->writeProlog( logger );
  boundary->synchronizer.initialize( boundary->distributedParticles );
  boundary->synchronizer.setCommunicator( this->communicator );
}

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
   const VectorType domainSize = parameters.getXyz< VectorType >( "domainSize" );
   const IndexVectorType domainGridDimension = TNL::ceil( ( domainSize - domainOrigin ) / searchRadius );

   // subdomain properties
   const VectorType subdomainOrigin = parametersDistributed.getXyz< VectorType >( subdomainKey + "origin" );
   const VectorType subdomainSize = parametersDistributed.getXyz< VectorType >(  subdomainKey + "size" );
   const IndexVectorType subdomainGridSize = TNL::ceil( subdomainSize / searchRadius );

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

   fluidOverlap->initialize( 0,
                             numberOfParticlesPerCell * overlapCellsCount,
                             searchRadius,
                             resizedSubdomainGridSize,
                             resizedSubdomainGridOrigin,
                             domainGridDimension,
                             domainOrigin,
                             logger );

   boundaryOverlap->initialize( 0,
                                numberOfParticlesPerCell * overlapCellsCount,
                                searchRadius,
                                resizedSubdomainGridSize,
                                resizedSubdomainGridOrigin,
                                domainGridDimension,
                                domainOrigin,
                                logger );
}

template< typename Model >
void
SPHMultiset_CFD< Model >::readParticlesFiles( TNL::Config::ParameterContainer& parameters, TNL::Logger& logger )
{
   logger.writeParameter( "Reading fluid particles:", parameters.getParameter< std::string >( "fluid-particles" ) );
   fluid->template readParticlesAndVariables< SimulationReaderType >(
      parameters.getParameter< std::string >( "fluid-particles" ) );
   logger.writeParameter( "Reading boundary particles:", parameters.getParameter< std::string >( "boundary-particles" ) );
   boundary->template readParticlesAndVariables< SimulationReaderType >(
      parameters.getParameter< std::string >( "boundary-particles" ) );

   // init open boundary patches
   const int numberOfBoundaryPatches = parameters.getParameter< int >( "openBoundaryPatches" );
   //TODO: I don't like that open boundary buffer are selected in compile time.
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ) {
      for( int i = 0; i < numberOfBoundaryPatches; i++ ) {
         std::string prefix = "buffer-" + std::to_string( i + 1 ) + "-";
         openBoundaryPatches[ i ]->template readParticlesAndVariables< SimulationReaderType >(
            parameters.getParameter< std::string >( prefix + "particles" ) );
      }
   }
}

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

template< typename Model >
void
SPHMultiset_CFD< Model >::performNeighborSearch( TNL::Logger& logger )
{
   //reset cell indices
   timeMeasurement.start( "search_reset" );
   fluid->particles->resetListWithIndices();
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ) {
      for( auto& openBoundaryPatch : openBoundaryPatches )
         openBoundaryPatch->particles->resetListWithIndices();
   }

#ifdef HAVE_MPI
#else
   if( timeStepping.getStep() == 0 )
#endif
      boundary->particles->resetListWithIndices();
   timeMeasurement.stop( "search_reset" );
   writeLog( logger, "Search - reset ...", "Done." );

   //compute cell indices
   timeMeasurement.start( "search_cellIndices" );
   fluid->particles->computeParticleCellIndices();
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ) {
      for( auto& openBoundaryPatch : openBoundaryPatches )
         openBoundaryPatch->particles->computeParticleCellIndices();
   }

#ifdef HAVE_MPI
#else
   if( timeStepping.getStep() == 0 )
#endif
      boundary->particles->computeParticleCellIndices();
   timeMeasurement.stop( "search_cellIndices" );
   writeLog( logger, "Search - compute cell indices ...", "Done." );

   //sort particles
   timeMeasurement.start( "search_sort" );
   fluid->sortParticles();
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ) {
      for( auto& openBoundaryPatch : openBoundaryPatches )
         openBoundaryPatch->sortParticles();
   }

#ifdef HAVE_MPI
#else
   if( timeStepping.getStep() == 0 )
#endif
      boundary->sortParticles();
   timeMeasurement.stop( "search_sort" );
   writeLog( logger, "Search - sort ...", "Done." );

   //update number of particles TODO: Do this in elegant way.
   // --- DEBUG ---
   logger.writeParameter( "Search - remove particles", "" );
   logger.writeParameter( "fluid.particles.getNumberOfParticlesToRemove()",  fluid->particles->getNumberOfParticlesToRemove(), 1 );
   logger.writeParameter( "fluid.synchronizer.getNumberOfRecvParticles()",  fluid->synchronizer.getNumberOfRecvParticles(), 1 );
   // -------------

   fluid->particles->setNumberOfParticles( fluid->particles->getNumberOfParticles()
                                           - fluid->particles->getNumberOfParticlesToRemove() );
   fluid->particles->setLastActiveParticle( fluid->particles->getLastActiveParticle()
                                            - fluid->particles->getNumberOfParticlesToRemove() );
   fluid->setLastActiveParticle( fluid->getLastActiveParticle() - fluid->particles->getNumberOfParticlesToRemove() );
   fluid->particles->setNumberOfParticlesToRemove( 0 );

   //update number of particles of boundary object
   boundary->particles->setNumberOfParticles( boundary->particles->getNumberOfParticles()
                                           - boundary->particles->getNumberOfParticlesToRemove() );
   boundary->particles->setLastActiveParticle( boundary->particles->getLastActiveParticle()
                                            - boundary->particles->getNumberOfParticlesToRemove() );
   boundary->setLastActiveParticle( boundary->getLastActiveParticle() - boundary->particles->getNumberOfParticlesToRemove() );
   boundary->particles->setNumberOfParticlesToRemove( 0 );

   //update number of particles TODO: Do this in elegant way.
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ) {
      for( auto& openBoundaryPatch : openBoundaryPatches ) {
         fluid->particles->setNumberOfParticles( fluid->particles->getNumberOfParticles()
                                                 - openBoundaryPatch->numberOfFluidParticlesToRemove );
         fluid->particles->setLastActiveParticle( fluid->particles->getLastActiveParticle()
                                                  - openBoundaryPatch->numberOfFluidParticlesToRemove );
         fluid->setLastActiveParticle( fluid->getLastActiveParticle() - openBoundaryPatch->numberOfFluidParticlesToRemove );
         openBoundaryPatch->numberOfFluidParticlesToRemove = 0;
      }
      writeLog( logger, "Search - resize ...", "Done." );
   }

   //assign particles to cells
   timeMeasurement.start( "search_toCells" );
   fluid->particles->particlesToCells();
   writeLog( logger, "Search - particles to cells - fluid ...", "Done." );
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ) {
      for( auto& openBoundaryPatch : openBoundaryPatches )
         openBoundaryPatch->particles->particlesToCells();
   }

#ifdef HAVE_MPI
#else
   if( timeStepping.getStep() == 0 )
#endif
      boundary->particles->particlesToCells();
   timeMeasurement.stop( "search_toCells" );
   writeLog( logger, "Search - particles to cells ...", "Done." );
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
SPHMultiset_CFD< Model >::applyOpenBC()
{
   for( long unsigned int i = 0; i < std::size( openBoundaryPatches ); i++ ) {
      //TODO Check if open boundary buffer is really open boundary buffer
      openBoundaryModel.applyOpenBoundary(
         timeStepping.getTimeStep(), fluid, openBoundaryPatches[ i ], openBoundaryPatches[ i ]->config );
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
   //Interacteraction between fluid and fluid and fluid and boundary
   model.interaction( fluid, boundary, modelParams );
   model.updateSolidBoundary( fluid, boundary, modelParams );

   //Interact between fluid and open boundary patches
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ) {
      for( long unsigned int i = 0; i < std::size( openBoundaryPatches ); i++ ) {
         openBoundaryPatches[ i ]->zone.updateParticlesInZone( fluid->particles );
         model.interactionWithOpenBoundary( fluid, openBoundaryPatches[ i ], modelParams );
         model.updateSolidBoundaryOpenBoundary( boundary, openBoundaryPatches[ i ], modelParams );
      }
   }

   //Finalize the interaction
   model.finalizeInteraction( fluid, boundary, modelParams );
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
SPHMultiset_CFD< Model >::synchronizeDistributedSimulation( bool writePoints )
{
   //NOTE: Much better syntax would be
   //synchronize( flud, fluidOverlap, s ynchronizer );

   fluid->synchronizeObject( fluidOverlap, writePoints );
   boundary->synchronizeObject( boundaryOverlap, writePoints );

}

template< typename Model >
void
SPHMultiset_CFD< Model >::resetOverlaps()
{
   //TODO: This should be paritcles method
   fluid->particles->removeParitclesOutOfDomain();
   boundary->particles->removeParitclesOutOfDomain();
}

#endif

template< typename Model >
void
SPHMultiset_CFD< Model >::save( TNL::Logger& logger, bool writeParticleCellIndex )
{
   if( verbose == "with-snapshot" )
      writeInfo( logger );

   const int step = timeStepping.getStep();
#ifdef HAVE_MPI
   std::string outputFileNameFluid = outputDirecotry + "/particles_rank" + std::to_string( TNL::MPI::GetRank() ) + "_" + std::to_string( step ) + "_fluid.vtk";
#else
   std::string outputFileNameFluid = outputDirecotry + "/particles" + std::to_string( step ) + "_fluid.vtk";
#endif
   fluid->template writeParticlesAndVariables< Writer >( outputFileNameFluid, writeParticleCellIndex );
   logger.writeParameter( "Saved:", outputFileNameFluid );

#ifdef HAVE_MPI
   std::string outputFileNameBound = outputDirecotry + "/particles_rank" + std::to_string( TNL::MPI::GetRank() ) + "_" + std::to_string( step ) + "_boundary.vtk";
#else
   std::string outputFileNameBound = outputDirecotry + "/particles" + std::to_string( step ) + "_boundary.vtk";
#endif
   boundary->template writeParticlesAndVariables< Writer >( outputFileNameBound, writeParticleCellIndex );
   logger.writeParameter( "Saved:", outputFileNameBound );

   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ) {
      for( auto& openBoundaryPatch : openBoundaryPatches ) {
         std::string outputFileNameOpenBound =
            outputDirecotry + "/particles" + std::to_string( step ) + "_" + openBoundaryPatch->parameters.identifier + ".vtk";
         openBoundaryPatch->template writeParticlesAndVariables< Writer >( outputFileNameOpenBound, writeParticleCellIndex );
         logger.writeParameter( "Saved:", outputFileNameOpenBound );
      }
   }

   // output simulation sensors to files
   simulationMonitor.save( logger );
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
   logger.writeParameter( "Output directory:", outputDirecotry );
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
   if( verbose == "full" ) {
      logger.writeParameter( "Fluid object - first active particle:", fluid->getFirstActiveParticle(), 1 );
      logger.writeParameter( "Fluid particles - first active particle:", fluid->particles->getFirstActiveParticle(), 1 );
      logger.writeParameter( "Fluid object - last active particle:", fluid->getLastActiveParticle(), 1 );
      logger.writeParameter( "Fluid particles - last active particle:", fluid->particles->getLastActiveParticle(), 1 );
   }
   logger.writeParameter( "Number of boundary particles:", boundary->getNumberOfParticles() );
   if( verbose == "full" ) {
      logger.writeParameter( "Boundary object - first active particle:", boundary->getFirstActiveParticle(), 1 );
      logger.writeParameter( "Boundary particles - first active particle:", boundary->particles->getFirstActiveParticle(), 1 );
      logger.writeParameter( "Boundary object - last active particle:", boundary->getLastActiveParticle(), 1 );
      logger.writeParameter( "Boundary particles - last active particle:", boundary->particles->getLastActiveParticle(), 1 );
   }
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ) {
      for( long unsigned int i = 0; i < openBoundaryPatches.size(); i++ ) {
         logger.writeParameter( "Number of buffer" + std::to_string( i + 1 ) + " particles:",
                                openBoundaryPatches[ i ]->getNumberOfParticles() );
         if( verbose == "full" ) {
            logger.writeParameter( "Patch " + openBoundaryPatches[ i ]->config.identifier + " object - first active particle:",
                                   openBoundaryPatches[ i ]->getFirstActiveParticle(),
                                   1 );
            logger.writeParameter( "Patch " + openBoundaryPatches[ i ]->config.identifier
                                      + " particles - first active particle:",
                                   openBoundaryPatches[ i ]->particles->getFirstActiveParticle(),
                                   1 );
            logger.writeParameter( "Patch " + openBoundaryPatches[ i ]->config.identifier + " object - last active particle:",
                                   openBoundaryPatches[ i ]->getLastActiveParticle(),
                                   1 );
            logger.writeParameter( "Patch " + openBoundaryPatches[ i ]->config.identifier
                                      + " particles - last active particle:",
                                   openBoundaryPatches[ i ]->particles->getLastActiveParticle(),
                                   1 );
         }
      }
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
SPHMultiset_CFD< Model >::writeEpilog( TNL::Logger& logger ) const noexcept
{
   logger.writeHeader( "SPH simulation successfully finished." );
   logger.writeCurrentTime( "Ended at:" );
   timeMeasurement.writeInfo( logger, timeStepping.getStep() );
}

}  //namespace SPH
}  //namespace TNL
