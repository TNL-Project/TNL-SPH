// SPDX-FileComment: This file is part of TNL-SPH (https://gitlab.com/tnl-project/tnl-sph)
// SPDX-License-Identifier: MIT

#pragma once

#include <cmath>

#include <SPH/solvers/SolverMultiSetBlockMultiresolutionLocalTimeStepping.h>
#include <SPH/tempFunctionsToConfigMultiresolution.h>

namespace TNL::SPH {

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::init( int argc, char* argv[] )
{
   auto& params = this->parameters;
   auto& log = this->logger;

   try {
      initialize< typename BaseType::SimulationType >( argc, argv, this->cliParams, this->cliConfig, params, this->config );
   }
   catch( ... ) {
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
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::getLevelRefinementFactor( int idx ) const
{
   const std::string subdomainKey = "subdomain-" + std::to_string( idx ) + "-";
   return this->parametersSubdomains.template getParameter< float >( subdomainKey + "refinement-factor" );
}

template< typename Model >
int
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::getLevelOverlapWidth( int idx ) const
{
   const float rfOwn = getLevelRefinementFactor( idx );
   for( const auto& iface : this->topology.getInterfacesOfSubdomain( idx ) )
      if( rfOwn < getLevelRefinementFactor( iface.neighborIdx ) )
         return 2;
   return 1;
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::initializeBlockBasedMultiResolutionSimulation()
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
   if( hasOpenBcFile || hasOpenBcInline ) {
      this->initOpenBoundaryPatches( params, log );

      this->timeMeasurement.addTimer( "extrapolate-openbc" );
      this->timeMeasurement.addTimer( "apply-openbc" );
   }

   const bool hasPeriodicBcFile = params.template getParameter< std::string >( "periodic-boundary-config" ) != "";
   const bool hasPeriodicBcInline = params.template getParameter< std::string >( "periodic-boundary" ) != "";
   if( hasPeriodicBcFile || hasPeriodicBcInline ) {
      this->initPeriodicBoundaryPatches( params, log );

      this->timeMeasurement.addTimer( "enforce-periodic-bc" );
      this->timeMeasurement.addTimer( "transfer-periodic-bc" );
      this->timeMeasurement.addTimer( "periodicity-fluid-updateZone", false );
      this->timeMeasurement.addTimer( "periodicity-boundary-updateZone", false );
   }

   this->modelParams.init( params );

   /* Mass nodes belong only to interface faces backed by the neighbor fluid's
      support: faces lying at the global domain boundary (tank walls, bottom) have
      no fluid beyond them, and flux bookkeeping anchored there degenerates the
      interface state. The domain extent and the excluded faces of each interface
      frame are computed by the multiresolution helpers and applied by the buffers. */
   const auto domBox = computeSubdomainUnionBox< SPHConfig::spaceDimension >( this->topology );
   const RealType globalCell = this->topology.getGlobalGrid().getSpaceSteps()[ 0 ];
   for( long unsigned int p = 0; p < multiresolutionBoundaryPatchesLTS.size(); p++ ) {
      const auto& iface = this->multiresolutionBoundaryPatchInterfaces[ p ];
      multiresolutionBoundaryPatchesLTS[ p ]->initMassNodesWithDomainExclusion(
         this->modelParams,
         iface.ownIdx,
         iface.neighborIdx,
         getLevelRefinementFactor( iface.ownIdx ),
         this->topology.getSubdomain( iface.ownIdx ),
         this->topology.getSubdomain( iface.neighborIdx ),
         globalCell,
         domBox.first,
         domBox.second );
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
      ( params.template getParameter< std::string >( "boundary-ghost-method" ) == "direct-from-source" )
         ? BoundaryGhostUpdate::DirectFromSource
         : BoundaryGhostUpdate::Interpolation;

   this->readParticlesFiles();

   const int numberOfMultiresolutionBuffers = multiresolutionBoundaryPatchesLTS.size();
   for( int i = 0; i < numberOfMultiresolutionBuffers; i++ ) {
      std::string bufferKey = "multiresolution-buffer-" + std::to_string( i ) + "-";
      if( this->parametersSubdomains.template getParameter< int >( bufferKey + "n" ) != 0 ) {
         const std::string mrbFileName =
            this->parametersSubdomains.template getParameter< std::string >( bufferKey + "particles" );
         log.writeParameter( "Reading multiresolution buffer particles:", mrbFileName );
         multiresolutionBoundaryPatchesLTS[ i ]->template readParticlesAndVariables< typename BaseType::Reader >( mrbFileName );
      }
   }

   /* Seed the midpoint anchors for any pre-seeded buffer particles: the first
      midpointUpdateBuffers runs before the first sync, and would otherwise blend
      with uninitialized rho_old/v_old (TNL arrays do not construct fundamental
      types on resize). No-op for patches with zero particles. */
   for( long unsigned int p = 0; p < multiresolutionBoundaryPatchesLTS.size(); p++ )
      multiresolutionBoundaryPatchesLTS[ p ]->initializeIntegratorVariables();

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
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::initParticleSets()
{
   const int nSubsets = this->topology.getNumberOfSubdomains();
   this->numberOfSubsets = nSubsets;
   this->fluidSets.resize( nSubsets );
   this->boundarySets.resize( nSubsets );

   for( int i = 0; i < nSubsets; i++ ) {
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
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::initMultiResolutionBoundaryPatches()
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
            this->fluidSets[ iface.ownIdx ]->getParticles(), this->fluidSets[ iface.neighborIdx ]->getParticles(), rf );

         this->multiresolutionBoundaryPatchInterfaces.push_back( iface );
         mrbIdx++;
      }
   }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::interact( int subdomainIdx )
{
   this->timeMeasurement.start( "interact" );

   this->model.updateSolidBoundary( this->fluidSets[ subdomainIdx ], this->boundarySets[ subdomainIdx ], this->modelParams );
   this->model.finalizeBoundaryInteraction(
      this->fluidSets[ subdomainIdx ], this->boundarySets[ subdomainIdx ], this->modelParams );

   this->model.interaction( this->fluidSets[ subdomainIdx ], this->boundarySets[ subdomainIdx ], this->modelParams );
   for( long unsigned int p = 0; p < multiresolutionBoundaryPatchesLTS.size(); p++ )
      if( this->multiresolutionBoundaryPatchInterfaces[ p ].ownIdx == subdomainIdx )
         this->model.interactionWithOpenBoundary(
            this->fluidSets[ subdomainIdx ], multiresolutionBoundaryPatchesLTS[ p ], this->modelParams );
   this->model.finalizeInteraction( this->fluidSets[ subdomainIdx ], this->boundarySets[ subdomainIdx ], this->modelParams );

   this->timeMeasurement.stop( "interact" );
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::integrateLevel( int idx )
{
   this->timeMeasurement.start( "integrate" );
   /* One Euler rebase per 20 steps of THIS level = one per 20*dt_L of physical
      time at any refinement ratio - the uniform drivers' semantics (their
      getStep()%20 counts fine-dt steps). Keying on the global window counter
      instead under-rebases fine levels by the refinement ratio (2x at 2:1: one
      rebase per 40 fine steps vs the control's 20), leaving the swap-based
      leapfrog's v/v_old and rho/rho_old desync unpurged twice as long; on a
      stiff, strongly oscillating flow (the impeller pocket) the desync then
      converts density excursions into net kinetic energy at the integrate
      kick. Rebasing several times per window is equally wrong (also an
      energy ratchet), hence the per-level-step key: at 1:1 it degenerates to
      the window keying bit-for-bit, at 2:1 it restores the control cadence. */
   const bool doResync = ( levelTimeStepping[ idx ].getStep() % 20 == 0 );
   this->integrator->integratStepVerlet(
      this->fluidSets[ idx ], this->boundarySets[ idx ], levelTimeStepping[ idx ], false, doResync );
   levelTimeStepping[ idx ].updateTimeStep();
   this->timeMeasurement.stop( "integrate" );
}

template< typename Model >
int
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::getNumberOfSubsteps( int level ) const
{
   return levelSubsteps[ level ];
}

template< typename Model >
template< typename SubstepUpdate >
void
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::recursiveSubdomainsLoop( SubstepUpdate&& subdomainSubstep )
{
   const int nLevels = this->numberOfSubsets;
   const RealType theta = this->modelParams.ltsTheta;

   std::function< void( int, int, int ) > advanceLevel = [ & ]( int level, int count, int base )
   {
      for( int substep = 0; substep < count; substep++ ) {
         const int g = base + substep;  // window-global substep index of this level
         const RealType dtLevel = levelTimeStepping[ level ].getTimeStep();

         if( g != 0 )
            refreshSubstepSearches( level );

         subdomainSubstep( level, substep, g, dtLevel );

         if( level + 1 < nLevels ) {
            const int ratio = getNumberOfSubsteps( level + 1 ) / getNumberOfSubsteps( level );
            advanceLevel( level + 1, ratio, g * ratio );

            /* R1 - exact coincidence refresh: the child invocation just brought the
               band's clock to the same physical time as this level; re-stamp the
               anchors so the following blends stay bracketed. Skipped after the last
               substep - the sync-point lifecycle re-grounds all bands anyway. */
            if( substep + 1 < count ) {
               refreshBuffersFromNeighbor( 1.0, level + 1 );
               stampBuffersAnchor( level + 1 );
            }
         }

         /* R2 - midpoint blend of the band owned by this level (the interface toward
            its coarser neighbor), firing once per invocation at the midpoint substep. */
         if( theta > 0 && level > 0 && substep == count / 2 - 1 )
            refreshBuffersFromNeighbor( theta, level );
      }
   };
   advanceLevel( 0, 1, 0 );
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::updateLevelTimeSteps()
{
   if constexpr( std::is_same_v< TimeStepping, ConstantTimeStep< SPHConfig > > ) {
      return;
   }
   else {
      this->timeMeasurement.start( "interact" );

      RealType rfMin = 1.f, rfMax = 1.f;
      for( int i = 0; i < this->numberOfSubsets; i++ ) {
         rfMin = std::min( rfMin, getLevelRefinementFactor( i ) );
         rfMax = std::max( rfMax, getLevelRefinementFactor( i ) );
      }

      RealType dtFinest = std::numeric_limits< RealType >::max();
      for( int i = 0; i < this->numberOfSubsets; i++ ) {
         auto levelParams = this->modelParams;
         levelParams.h = this->modelParams.h * getLevelRefinementFactor( i );
         const RealType dtLevel = levelTimeStepping[ i ].getTimeStepCFL( this->fluidSets[ i ], levelParams );

         /* implied finest-level step: the most restrictive level governs the window */
         dtFinest = std::min( dtFinest, dtLevel * rfMin / getLevelRefinementFactor( i ) );
      }
      dtFinest = std::max( std::min( dtFinest, this->modelParams.dtInit ), this->modelParams.dtMin );

      this->timeStepping.setTimeStep( dtFinest * rfMax / rfMin );
      for( int i = 0; i < this->numberOfSubsets; i++ )
         levelTimeStepping[ i ].setTimeStep( dtFinest * getLevelRefinementFactor( i ) / rfMin );

      this->timeMeasurement.stop( "interact" );
   }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::advectSubstepBuffers( RealType dt, int fineLevelIdx )
{
   for( long unsigned int p = 0; p < multiresolutionBoundaryPatchesLTS.size(); p++ ) {
      const auto& iface = this->multiresolutionBoundaryPatchInterfaces[ p ];
      const int finerIdx = getLevelRefinementFactor( iface.ownIdx ) < getLevelRefinementFactor( iface.neighborIdx )
                            ? iface.ownIdx
                            : iface.neighborIdx;
      if( finerIdx != fineLevelIdx )
         continue;
      multiresolutionBoundaryPatchesLTS[ p ]->advectBufferParticles( dt );
   }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::refreshBuffersFromNeighbor( RealType theta, int fineLevelIdx )
{
   for( long unsigned int p = 0; p < multiresolutionBoundaryPatchesLTS.size(); p++ ) {
      const auto& iface = this->multiresolutionBoundaryPatchInterfaces[ p ];
      const int finerIdx = getLevelRefinementFactor( iface.ownIdx ) < getLevelRefinementFactor( iface.neighborIdx )
                            ? iface.ownIdx
                            : iface.neighborIdx;
      if( finerIdx != fineLevelIdx || iface.ownIdx != finerIdx )
         continue;
      this->fluidSets[ iface.neighborIdx ]->searchForNeighbors();
      multiresolutionBoundaryPatchesLTS[ p ]->refreshVariablesFromNeighbor(
         this->fluidSets[ iface.neighborIdx ], this->modelParams, theta );
   }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::stampBuffersAnchor( int fineLevelIdx )
{
   for( long unsigned int p = 0; p < multiresolutionBoundaryPatchesLTS.size(); p++ ) {
      const auto& iface = this->multiresolutionBoundaryPatchInterfaces[ p ];
      const int finerIdx = getLevelRefinementFactor( iface.ownIdx ) < getLevelRefinementFactor( iface.neighborIdx )
                            ? iface.ownIdx
                            : iface.neighborIdx;
      if( finerIdx != fineLevelIdx || iface.ownIdx != finerIdx )
         continue;
      multiresolutionBoundaryPatchesLTS[ p ]->initializeIntegratorVariables();
   }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::accumulateInterfaceMassesFromNeighbor(
   int neighborIdx,
   RealType dt )
{
   for( long unsigned int p = 0; p < multiresolutionBoundaryPatchesLTS.size(); p++ )
      if( this->multiresolutionBoundaryPatchInterfaces[ p ].neighborIdx == neighborIdx )
         multiresolutionBoundaryPatchesLTS[ p ]->accumulateMasses( this->fluidSets[ neighborIdx ], this->modelParams, dt );
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::refreshSubstepSearches( int fineIdx )
{
   // positions of the fine fluid and of all buffer bands changed during the last substep
   if( ! multiresolutionBoundaryPatchesLTS.empty() )
      multiresolutionBoundaryPatchesLTS[ 0 ]->removeParticlesOutOfGrid( this->fluidSets[ fineIdx ]->getParticles() );
   this->fluidSets[ fineIdx ]->searchForNeighbors();

   for( long unsigned int p = 0; p < multiresolutionBoundaryPatchesLTS.size(); p++ ) {
      multiresolutionBoundaryPatchesLTS[ p ]->removeParticlesOutOfGrid(
         multiresolutionBoundaryPatchesLTS[ p ]->getParticles() );
      multiresolutionBoundaryPatchesLTS[ p ]->searchForNeighbors();
   }

   /* The boundary sets of the refreshed level must be re-searched as well: the fine
      fluid was just re-sorted and the boundary markers may move inside substeps, so
      keeping the window-start boundary search leaves boundary-fluid walks built on
      stale structures. Without it the fine level slowly loses wall support (fluid
      leaks through walls, mostly at the bottom, and is culled as out-of-domain),
      degenerating the interface until the mass-flux bookkeeping explodes - the
      uniform-stepping drivers re-search the boundary every step for this reason. */
   if( fineIdx >= 0 && fineIdx < int( this->boundarySets.size() ) )
      this->boundarySets[ fineIdx ]->searchForNeighbors();
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::syncMultiresolutionUpdate()
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
   this->writeLog( "Update multiresolution BC (LTS)...", "Done." );
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::save( bool writeParticleCellIndex )
{
   if( ( this->verbose == "with-snapshot" ) || ( this->verbose == "full" ) )
      this->writeInfo();

   const RealType time = this->timeStepping.getTime();

   for( int i = 0; i < this->numberOfSubsets; i++ ) {
#ifdef HAVE_MPI
      std::string outputFileNameFluid = this->outputDirectory + "/fluid_rank" + std::to_string( TNL::MPI::GetRank() ) + "_"
                                      + std::to_string( time ) + "_particles.vtk";
#else
      std::string outputFileNameFluid =
         this->outputDirectory + "/fluid_subdomain" + std::to_string( i ) + "_" + std::to_string( time ) + "_particles.vtk";
#endif
      this->fluidSets[ i ]->template writeParticlesAndVariables< typename BaseType::Writer >(
         outputFileNameFluid, writeParticleCellIndex );
      this->logger.writeParameter( "Saved:", outputFileNameFluid );

#ifdef HAVE_MPI
      std::string outputFileNameBound = this->outputDirectory + "/boundary_rank" + std::to_string( TNL::MPI::GetRank() ) + "_"
                                      + std::to_string( time ) + "_particles.vtk";
#else
      std::string outputFileNameBound =
         this->outputDirectory + "/boundary_subdomain" + std::to_string( i ) + "_" + std::to_string( time ) + "_particles.vtk";
#endif
      this->boundarySets[ i ]->template writeParticlesAndVariables< typename BaseType::Writer >(
         outputFileNameBound, writeParticleCellIndex );
      this->logger.writeParameter( "Saved:", outputFileNameBound );

#ifdef HAVE_MPI
      std::string outputFileNameGrid =
         this->outputDirectory + "/grid_rank" + std::to_string( TNL::MPI::GetRank() + 1 ) + "_" + std::string( time ) + ".vtk";
#else
      std::string outputFileNameGrid =
         this->outputDirectory + "/grid_subdomain" + std::to_string( i ) + "_" + std::to_string( time ) + ".vtk";
#endif
      TNL::Particles::Writers::writeBackgroundGrid(
         outputFileNameGrid,
         this->fluidSets[ i ]->getParticles()->getDimensions(),
         this->fluidSets[ i ]->getParticles()->getOrigin(),
         this->fluidSets[ i ]->getParticles()->getSearchRadius() );
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
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::writeProlog( bool writeSystemInformation ) noexcept
{
   BaseType::writeProlog( writeSystemInformation );

   auto& log = this->logger;
   log.writeHeader( "Local timesteps." );
   for( int i = 0; i < this->numberOfSubsets; i++ ) {
      log.writeParameter( "Time step (subdomain " + std::to_string( i ) + "):", levelTimeStepping[ i ].getTimeStep() );
      log.writeParameter( "Substeps per global step (subdomain " + std::to_string( i ) + "):", levelSubsteps[ i ] );
   }
   log.writeSeparator();

   for( long unsigned int i = 0; i < multiresolutionBoundaryPatchesLTS.size(); i++ ) {
      log.writeHeader( "Multiresolution boundary buffer (LTS) " + std::to_string( i ) + "." );
      multiresolutionBoundaryPatchesLTS[ i ]->writeProlog( log, i );
   }

   // Structured initialization report (JSON) - the base-class prolog above already
   // wrote the block-multiresolution document; rewrite it with the LTS buffer bands
   // so the final file carries this solver's own buffers.
   nlohmann::json initReport = buildInitReportDocument( *this );
   for( long unsigned int i = 0; i < multiresolutionBoundaryPatchesLTS.size(); i++ )
      initReport[ "multiresolution_buffers" ].push_back( multiresolutionBoundaryPatchesLTS[ i ]->writeInitReport() );
   initReport[ "decomposition" ] = this->topology.writeInitReport();
   writeInitReportFile( initReport, this->outputDirectory );
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::writeInfo() noexcept
{
   BaseType::writeInfo();

   auto& log = this->logger;
   for( long unsigned int i = 0; i < multiresolutionBoundaryPatchesLTS.size(); i++ )
      log.writeParameter(
         "Number of mr-buffer (LTS) " + std::to_string( i ) + " particles:",
         multiresolutionBoundaryPatchesLTS[ i ]->getNumberOfParticles() );
   log.writeSeparator();
}

}  // namespace TNL::SPH
