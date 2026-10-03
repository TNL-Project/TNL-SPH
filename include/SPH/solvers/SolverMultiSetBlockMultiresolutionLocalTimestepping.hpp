#include "SolverMultiSetBlockMultiresolutionLocalTimestepping.h"
#include "../tempFunctionsToConfigMultiresolution.h"

#include <cmath>

namespace TNL {
namespace SPH {

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimestepping< Model >::init( int argc, char* argv[] )
{
   auto& params = this->parameters;
   auto& log = this->logger;

   try {
      initialize< typename BaseType::SimulationType >( argc, argv, this->cliParams, this->cliConfig, params, this->config );
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
   // LTS substepping is currently implemented for the block path only
   this->initializeDistributedSimulation();
#else
   initializeBlockBasedMultiResolutionSimulation();
#endif
   log.writeHeader( "SPH simulation successfully initialized." );
}

template< typename Model >
float
SolverMultiSetBlockMultiresolutionLocalTimestepping< Model >::getLevelRefinementFactor( int idx ) const
{
   const std::string subdomainKey = "subdomain-" + std::to_string( idx ) + "-";
   return this->parametersSubdomains.template getParameter< float >( subdomainKey + "refinement-factor" );
}

template< typename Model >
int
SolverMultiSetBlockMultiresolutionLocalTimestepping< Model >::getLevelOverlapWidth( int idx ) const
{
   const float rfOwn = getLevelRefinementFactor( idx );
   for( const auto& iface : this->topology.getInterfacesOfSubdomain( idx ) )
      if( rfOwn < getLevelRefinementFactor( iface.neighborIdx ) )
         return 2;
   return 1;
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimestepping< Model >::initializeBlockBasedMultiResolutionSimulation()
{
   auto& params = this->parameters;
   auto& log = this->logger;

   this->numberOfSubsets = params.template getParameter< int >( "numberOfSubdomains" );
   const std::string configSubdomainsPath = params.template getParameter< std::string >( "subdomains-config" );
   const std::string subdomainsInline = params.template getParameter< std::string >( "subdomains" );
    for( int subset = 0; subset < this->numberOfSubsets; subset++ )
       TNL::SPH::configSubdomain( subset, this->configSubdomains );
    for( int bufferIdx = 0; bufferIdx < 2 * this->numberOfSubsets; bufferIdx++ ) {
       TNL::SPH::configMultiresolutionBuffer( bufferIdx, this->configSubdomains );
       TNL::SPH::configBoundaryGhostBuffer( bufferIdx, this->configSubdomains );
    }
   parseDistributedConfig( configSubdomainsPath, this->parametersSubdomains, this->configSubdomains, log, subdomainsInline );

   this->topology.loadFromConfig( params, this->parametersSubdomains );
   this->topology.finalizeLinear();
   initParticleSets();
   initMultiResolutionBoundaryPatches();
   this->timeMeasurement.addTimer( "multiresolution-update" );

   const bool hasOpenBcFile = params.template getParameter< std::string >( "open-boundary-config" ) != "";
   const bool hasOpenBcInline = params.template getParameter< std::string >( "open-boundary" ) != "";
   if( hasOpenBcFile || hasOpenBcInline ){
      this->initOpenBoundaryPatches( params, log );

      this->timeMeasurement.addTimer( "extrapolate-openbc" );
      this->timeMeasurement.addTimer( "apply-openbc" );
   }

    const bool hasPeriodicBcFile = params.template getParameter< std::string >( "periodic-boundary-config" ) != "";
    const bool hasPeriodicBcInline = params.template getParameter< std::string >( "periodic-boundary" ) != "";
    if( hasPeriodicBcFile || hasPeriodicBcInline ){
       this->initPeriodicBoundaryPatches( params, log );

       this->timeMeasurement.addTimer( "enforce-periodic-bc" );
       this->timeMeasurement.addTimer( "transfer-periodic-bc" );
       this->timeMeasurement.addTimer( "periodicity-fluid-updateZone", false );
       this->timeMeasurement.addTimer( "periodicity-boundary-updateZone", false );
    }

    this->modelParams.init( params );
    for( long unsigned int p = 0; p < multiresolutionBoundaryPatchesLTS.size(); p++ ){
      const auto& iface = this->multiresolutionBoundaryPatchInterfaces[ p ];
      const float refinementFactor = getLevelRefinementFactor( iface.ownIdx );
      multiresolutionBoundaryPatchesLTS[ p ]->initMassNodes( this->modelParams, iface.ownIdx, iface.neighborIdx, refinementFactor );
   }

   /* Per-level time steps: the configured time step is the FINEST level's step (as written
      by the init script from CFL*h_min/c0); coarser levels scale with the refinement-factor
      ratio. The global clock advances with the coarsest level's step. */
   const RealType dtFinest = params.template getParameter< RealType >( "initial-time-step" );
   RealType rfMin = 1.f;
   for( int i = 0; i < this->numberOfSubsets; i++ )
      rfMin = std::min( rfMin, getLevelRefinementFactor( i ) );

   RealType dtGlobal = 0.f;
   for( int i = 0; i < this->numberOfSubsets; i++ )
      dtGlobal = std::max( dtGlobal, dtFinest * getLevelRefinementFactor( i ) / rfMin );

   this->timeStepping.setTimeStep( dtGlobal );
   this->timeStepping.setEndTime( params.template getParameter< RealType >( "final-time" ) );
   this->timeStepping.addOutputTimer( "save_results", params.template getParameter< RealType >( "snapshot-period" ) );

   levelTimeStepping.resize( this->numberOfSubsets );
   levelSubsteps.resize( this->numberOfSubsets );
   for( int i = 0; i < this->numberOfSubsets; i++ ) {
      const RealType dtLevel = dtFinest * getLevelRefinementFactor( i ) / rfMin;
      levelTimeStepping[ i ] = this->timeStepping;
      levelTimeStepping[ i ].setTimeStep( dtLevel );
      levelSubsteps[ i ] = static_cast< int >( std::lround( dtGlobal / dtLevel ) );
   }

     this->boundaryGhostUpdate =
        ( params.template getParameter< std::string >( "boundary-ghost-method" ) == "direct-from-source" ) ?
           BoundaryGhostUpdate::DirectFromSource :
           BoundaryGhostUpdate::Interpolation;

     this->readParticlesFiles();

     const int numberOfMultiresolutionBuffers = multiresolutionBoundaryPatchesLTS.size();
     for( int i = 0; i < numberOfMultiresolutionBuffers; i++ ) {
        std::string bufferKey = "multiresolution-buffer-" + std::to_string( i ) + "-";
        if( this->parametersSubdomains.template getParameter< int >( bufferKey + "n" ) != 0 ) {
           const std::string mrbFileName = this->parametersSubdomains.template getParameter< std::string >( bufferKey + "particles" );
           log.writeParameter( "Reading multiresolution buffer particles:", mrbFileName );
           multiresolutionBoundaryPatchesLTS[ i ]->template readParticlesAndVariables< typename BaseType::Reader >( mrbFileName );
        }
     }

     this->initBoundaryGhosts();

    log.writeSeparator();
    const bool hasMeasuretoolFile = params.template getParameter< std::string >( "measuretool-config" ) != "";
    const bool hasMeasuretoolInline = params.template getParameter< std::string >( "measuretool" ) != "";
    if( hasMeasuretoolFile || hasMeasuretoolInline ) {
       log.writeParameter( "Simulation monitor initialization.", "" );
       this->simulationMonitor.init( params, this->timeStepping, log );
       this->simulationMonitor.setupVolumetricFlowRateZones( this->fluidSets[ 0 ] );
       log.writeParameter( "Simulation monitor initialization.", "Done." );
    }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimestepping< Model >::initParticleSets()
{
   const int nSubsets = this->topology.getNumberOfSubdomains();
   this->numberOfSubsets = nSubsets;
   this->fluidSets.resize( nSubsets );
   this->boundarySets.resize( nSubsets );

   for( int i = 0; i < nSubsets; i++ ){
      std::string subdomainKey = "subdomain-" + std::to_string( i ) + "-";
      const int overlapWidth = getLevelOverlapWidth( i );

      this->fluidSets[ i ]->initializeAsDistributed(
            this->parametersSubdomains.template getParameter< int >( subdomainKey + "fluid_n" ),
            this->parametersSubdomains.template getParameter< int >( subdomainKey + "fluid_n_allocated" ),
            this->topology.getLocalGrid( i ),
            this->topology.getLocalOriginCoordinates( i ),
            this->topology.getGlobalGrid(),
            overlapWidth );
      if constexpr( ParticlesType::specifySearchedSetExplicitly() == true )
         this->fluidSets[ i ]->getParticles()->setParticleSetLabel( 0 );

      this->boundarySets[ i ]->initializeAsDistributed(
            this->parametersSubdomains.template getParameter< int >( subdomainKey + "boundary_n" ),
            this->parametersSubdomains.template getParameter< int >( subdomainKey + "boundary_n_allocated" ),
            this->topology.getLocalGrid( i ),
            this->topology.getLocalOriginCoordinates( i ),
            this->topology.getGlobalGrid(),
            overlapWidth );
      if constexpr( ParticlesType::specifySearchedSetExplicitly() == true )
         this->boundarySets[ i ]->getParticles()->setParticleSetLabel( 1 );
   }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimestepping< Model >::initMultiResolutionBoundaryPatches()
{
   const int mrbCount = this->topology.getNumberOfSubdomainInterfaces();
   multiresolutionBoundaryPatchesLTS.resize( mrbCount );

   int mrbIdx = 0;
   for( int i = 0; i < this->topology.getNumberOfSubdomains(); i++ ) {
      const float rf = getLevelRefinementFactor( i );

      for( const auto& iface : this->topology.getInterfacesOfSubdomain( i ) ) {

         const std::string bufferKey = "multiresolution-buffer-" + std::to_string( mrbIdx ) + "-";
         const int mrbNumOfPtcs = this->parametersSubdomains.template getParameter< int >( bufferKey + "n" );

         const int overlapWidth = ( rf < getLevelRefinementFactor( iface.neighborIdx ) ) ? 2 : 1;
         const RealType resolutionScale = overlapWidth * std::pow( rf, -int( SPHConfig::spaceDimension ) );
         const int mrbNumOfAllocPtcs =
            std::max( int( 0.1 * this->getTotalFluidParticlesCount() * resolutionScale ), 2 * mrbNumOfPtcs );

         multiresolutionBoundaryPatchesLTS[ mrbIdx ]->initializeAsDistributed(
            mrbNumOfPtcs,
            mrbNumOfAllocPtcs,
            this->topology.getLocalGrid( i ),
            this->topology.getLocalOriginCoordinates( i ),
            this->topology.getGlobalGrid(),
            overlapWidth );

         multiresolutionBoundaryPatchesLTS[ mrbIdx ]->initZones(
            this->fluidSets[ iface.ownIdx ]->getParticles(),
            this->fluidSets[ iface.neighborIdx ]->getParticles(),
            rf );

         this->multiresolutionBoundaryPatchInterfaces.push_back( iface );
         mrbIdx++;
      }
   }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimestepping< Model >::interactCoarseLevel( int coarseIdx )
{
   this->timeMeasurement.start( "interact" );

   this->model.updateSolidBoundary( this->fluidSets[ coarseIdx ], this->boundarySets[ coarseIdx ], this->modelParams );
   this->model.finalizeBoundaryInteraction( this->fluidSets[ coarseIdx ], this->boundarySets[ coarseIdx ], this->modelParams );

   this->model.interaction( this->fluidSets[ coarseIdx ], this->boundarySets[ coarseIdx ], this->modelParams );
   for( long unsigned int p = 0; p < multiresolutionBoundaryPatchesLTS.size(); p++ )
      if( this->multiresolutionBoundaryPatchInterfaces[ p ].ownIdx == coarseIdx )
         this->model.interactionWithOpenBoundary(
               this->fluidSets[ coarseIdx ], multiresolutionBoundaryPatchesLTS[ p ], this->modelParams );
   this->model.finalizeInteraction( this->fluidSets[ coarseIdx ], this->boundarySets[ coarseIdx ], this->modelParams );

   this->timeMeasurement.stop( "interact" );
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimestepping< Model >::interactFineSubstep( int fineIdx, bool computeActiveLayer )
{
   this->timeMeasurement.start( "interact" );

   this->model.updateSolidBoundary( this->fluidSets[ fineIdx ], this->boundarySets[ fineIdx ], this->modelParams );
   this->model.finalizeBoundaryInteraction( this->fluidSets[ fineIdx ], this->boundarySets[ fineIdx ], this->modelParams );

   this->model.interaction( this->fluidSets[ fineIdx ], this->boundarySets[ fineIdx ], this->modelParams );
   for( long unsigned int p = 0; p < multiresolutionBoundaryPatchesLTS.size(); p++ ) {
      if( this->multiresolutionBoundaryPatchInterfaces[ p ].ownIdx != fineIdx )
         continue;
      this->model.interactionWithOpenBoundary(
            this->fluidSets[ fineIdx ], multiresolutionBoundaryPatchesLTS[ p ], this->modelParams );
      if( computeActiveLayer )
         this->model.interactionActiveBufferLayer(
               this->fluidSets[ fineIdx ], multiresolutionBoundaryPatchesLTS[ p ], this->modelParams );
   }
   this->model.finalizeInteraction( this->fluidSets[ fineIdx ], this->boundarySets[ fineIdx ], this->modelParams );
   if( computeActiveLayer )
      for( long unsigned int p = 0; p < multiresolutionBoundaryPatchesLTS.size(); p++ )
         if( this->multiresolutionBoundaryPatchInterfaces[ p ].ownIdx == fineIdx )
            this->model.finalizeInteractionActiveBufferLayer( multiresolutionBoundaryPatchesLTS[ p ], this->modelParams );

   this->timeMeasurement.stop( "interact" );
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimestepping< Model >::integrateLevel( int idx )
{
   this->timeMeasurement.start( "integrate" );
   this->integrator->integratStepVerlet( this->fluidSets[ idx ], this->boundarySets[ idx ], levelTimeStepping[ idx ], false );
   levelTimeStepping[ idx ].updateTimeStep();
   this->timeMeasurement.stop( "integrate" );
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimestepping< Model >::advectSubstepBuffers( RealType dt, bool skipLayer1,
                                                                                    int fineLevelIdx )
{
   for( long unsigned int p = 0; p < multiresolutionBoundaryPatchesLTS.size(); p++ ) {
      const auto& iface = this->multiresolutionBoundaryPatchInterfaces[ p ];
      const int finerIdx = getLevelRefinementFactor( iface.ownIdx ) < getLevelRefinementFactor( iface.neighborIdx )
                             ? iface.ownIdx : iface.neighborIdx;
      if( finerIdx != fineLevelIdx )
         continue;
      multiresolutionBoundaryPatchesLTS[ p ]->advectBufferParticles( dt, skipLayer1 );
   }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimestepping< Model >::integrateBuffersLayer1( RealType dt, bool useEulerRestart,
                                                                                      int fineLevelIdx )
{
   for( long unsigned int p = 0; p < multiresolutionBoundaryPatchesLTS.size(); p++ ) {
      const auto& iface = this->multiresolutionBoundaryPatchInterfaces[ p ];
      const int finerIdx = getLevelRefinementFactor( iface.ownIdx ) < getLevelRefinementFactor( iface.neighborIdx )
                             ? iface.ownIdx : iface.neighborIdx;
      if( finerIdx != fineLevelIdx )
         continue;
      multiresolutionBoundaryPatchesLTS[ p ]->integrateLayer1( dt, useEulerRestart );
   }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimestepping< Model >::midpointUpdateBuffers( RealType theta, int fineLevelIdx )
{
   for( long unsigned int p = 0; p < multiresolutionBoundaryPatchesLTS.size(); p++ ) {
      const auto& iface = this->multiresolutionBoundaryPatchInterfaces[ p ];
      const int finerIdx = getLevelRefinementFactor( iface.ownIdx ) < getLevelRefinementFactor( iface.neighborIdx )
                              ? iface.ownIdx : iface.neighborIdx;
      if( finerIdx != fineLevelIdx || iface.ownIdx != finerIdx )
         continue;
      this->fluidSets[ iface.neighborIdx ]->searchForNeighbors();
      multiresolutionBoundaryPatchesLTS[ p ]->updateVariablesMidpoint(
            this->fluidSets[ iface.neighborIdx ], this->modelParams, theta );
   }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimestepping< Model >::accumulateInterfaceMassesFromNeighbor( int neighborIdx,
                                                                                                      RealType dt )
{
   for( long unsigned int p = 0; p < multiresolutionBoundaryPatchesLTS.size(); p++ )
      if( this->multiresolutionBoundaryPatchInterfaces[ p ].neighborIdx == neighborIdx )
         multiresolutionBoundaryPatchesLTS[ p ]->accumulateMasses( this->fluidSets[ neighborIdx ], this->modelParams, dt );
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimestepping< Model >::refreshSubstepSearches( int fineIdx )
{
   // positions of the fine fluid and of all buffer bands changed during the last substep
   if( ! multiresolutionBoundaryPatchesLTS.empty() )
      multiresolutionBoundaryPatchesLTS[ 0 ]->removeParticlesOutOfGrid( this->fluidSets[ fineIdx ]->getParticles() );
   this->fluidSets[ fineIdx ]->searchForNeighbors();

   for( long unsigned int p = 0; p < multiresolutionBoundaryPatchesLTS.size(); p++ ) {
      multiresolutionBoundaryPatchesLTS[ p ]->removeParticlesOutOfGrid( multiresolutionBoundaryPatchesLTS[ p ]->getParticles() );
      multiresolutionBoundaryPatchesLTS[ p ]->searchForNeighbors();
   }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimestepping< Model >::syncMultiresolutionUpdate()
{
   this->timeMeasurement.start( "multiresolution-update" );

   for( long unsigned int p = 0; p < multiresolutionBoundaryPatchesLTS.size(); p++ ) {
      const auto& iface = this->multiresolutionBoundaryPatchInterfaces[ p ];
      multiresolutionBoundaryPatchesLTS[ p ]->syncUpdateInterfaceBuffer(
            this->fluidSets[ iface.ownIdx ],
            this->fluidSets[ iface.neighborIdx ],
            this->modelParams,
            this->timeStepping.getTimeStep() );
   }

   this->timeMeasurement.stop( "multiresolution-update" );
   this->writeLog( "Update multiresolution BC (LTS)...", "Done.");
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimestepping< Model >::save( bool writeParticleCellIndex )
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

#ifdef HAVE_MPI
      std::string outputFileNameGrid = this->outputDirectory + "/grid_rank" + std::to_string( TNL::MPI::GetRank() + 1 ) + "_" + std::string( time ) + ".vtk";
#else
      std::string outputFileNameGrid = this->outputDirectory + "/grid_subdomain" + std::to_string( i ) + "_" + std::to_string( time ) + ".vtk";
#endif
      TNL::Particles::Writers::writeBackgroundGrid( outputFileNameGrid, this->fluidSets[ i ]->getParticles()->getDimensions(), this->fluidSets[ i ]->getParticles()->getOrigin(), this->fluidSets[ i ]->getParticles()->getSearchRadius() );
      this->logger.writeParameter( "Saved:", outputFileNameGrid );

      this->simulationMonitor.save( this->logger );
   }

   for( long unsigned int p = 0; p < this->multiresolutionBoundaryPatchInterfaces.size(); p++ ) {
      const auto& iface = this->multiresolutionBoundaryPatchInterfaces[ p ];
      std::string outputFileNameMultiresolutionBound =
         this->outputDirectory + "/multiresolutionBoundaryPatch_subdomain" + std::to_string( iface.ownIdx ) + "_to_"
            + std::to_string( iface.neighborIdx ) + "_" + std::to_string( time ) + "_particles.vtk";
      multiresolutionBoundaryPatchesLTS[ p ]->template writeParticlesAndVariables< typename BaseType::Writer >(
            outputFileNameMultiresolutionBound, writeParticleCellIndex );
      this->logger.writeParameter( "Saved:", outputFileNameMultiresolutionBound );
   }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimestepping< Model >::writeProlog( bool writeSystemInformation ) noexcept
{
   BaseType::writeProlog( writeSystemInformation );

   auto& log = this->logger;
   log.writeHeader( "Local timesteps." );
   for( int i = 0; i < this->numberOfSubsets; i++ ) {
      log.writeParameter( "Time step (subdomain " + std::to_string( i ) + "):",
                          levelTimeStepping[ i ].getTimeStep() );
      log.writeParameter( "Substeps per global step (subdomain " + std::to_string( i ) + "):",
                          levelSubsteps[ i ] );
   }
   log.writeSeparator();

   for( long unsigned int i = 0; i < multiresolutionBoundaryPatchesLTS.size(); i++ ) {
      log.writeHeader( "Multiresolution boundary buffer (LTS) " + std::to_string( i ) + "." );
      multiresolutionBoundaryPatchesLTS[ i ]->writeProlog( log, i );
   }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimestepping< Model >::writeInfo() noexcept
{
   BaseType::writeInfo();

   auto& log = this->logger;
   for( long unsigned int i = 0; i < multiresolutionBoundaryPatchesLTS.size(); i++ )
      log.writeParameter( "Number of mr-buffer (LTS) " + std::to_string( i ) + " particles:",
                          multiresolutionBoundaryPatchesLTS[ i ]->getNumberOfParticles() );
   log.writeSeparator();
}

} // namespace SPH
} // namespace TNL
