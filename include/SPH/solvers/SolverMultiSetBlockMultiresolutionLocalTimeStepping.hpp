// SPDX-FileComment: This file is part of TNL-SPH (https://gitlab.com/tnl-project/tnl-sph)
// SPDX-License-Identifier: MIT

#pragma once

#include <cmath>

#include <SPH/solvers/SolverMultiSetBlockMultiresolutionLocalTimeStepping.h>
#include <SPH/MultiresolutionConfigHelpers.h>

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
   this->caseName = params.template getParameter< std::string >( "caseName" );
   this->verbose = params.template getParameter< std::string >( "verboseIntensity" );
   this->outputDirectory = params.template getParameter< std::string >( "outputDirectory" );
   this->particlesFormat = params.template getParameter< std::string >( "particlesFormat" );

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
   const std::string subdomainKey = "subdomain" + std::to_string( idx );
   return this->parametersSubdomains.template getParameter< float >( subdomainKey + "RefinementFactor" );
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
   const std::string configSubdomainsPath = params.template getParameter< std::string >( "subdomainsConfig" );
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

   const bool hasOpenBcFile = params.template getParameter< std::string >( "openBoundaryConfig" ) != "";
   const bool hasOpenBcInline = params.template getParameter< std::string >( "openBoundary" ) != "";
   if( hasOpenBcFile || hasOpenBcInline ) {
      this->initOpenBoundaryPatches( params, log );

      this->timeMeasurement.addTimer( "extrapolate-openbc" );
      this->timeMeasurement.addTimer( "apply-openbc" );
   }

   const bool hasPeriodicBcFile = params.template getParameter< std::string >( "periodicBoundaryConfig" ) != "";
   const bool hasPeriodicBcInline = params.template getParameter< std::string >( "periodicBoundary" ) != "";
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
   for( long unsigned int p = 0; p < this->multiresolutionBoundaryPatches.size(); p++ ) {
      const auto& iface = this->multiresolutionBoundaryPatchInterfaces[ p ];
      this->multiresolutionBoundaryPatches[ p ]->initMassNodesWithDomainExclusion(
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
   const RealType dtFinest = params.template getParameter< RealType >( "initialTimeStep" );
   RealType rfMin = 1.f;
   for( int i = 0; i < this->numberOfSubsets; i++ )
      rfMin = std::min( rfMin, getLevelRefinementFactor( i ) );

   RealType dtGlobal = 0.f;
   for( int i = 0; i < this->numberOfSubsets; i++ )
      dtGlobal = std::max( dtGlobal, dtFinest * getLevelRefinementFactor( i ) / rfMin );

   this->timeStepping.setTimeStep( dtGlobal );
   this->timeStepping.setEndTime( params.template getParameter< RealType >( "finalTime" ) );
   this->timeStepping.addOutputTimer( "save_results", params.template getParameter< RealType >( "snapshotPeriod" ) );

   levelTimeStepping.resize( this->numberOfSubsets );
   levelSubsteps.resize( this->numberOfSubsets );
   for( int i = 0; i < this->numberOfSubsets; i++ ) {
      const RealType dtLevel = dtFinest * getLevelRefinementFactor( i ) / rfMin;
      levelTimeStepping[ i ] = this->timeStepping;
      levelTimeStepping[ i ].setTimeStep( dtLevel );
      levelSubsteps[ i ] = static_cast< int >( std::lround( dtGlobal / dtLevel ) );
   }

   this->boundaryGhostUpdate = ( params.template getParameter< std::string >( "boundaryGhostMethod" ) == "direct-from-source" )
                                ? BoundaryGhostUpdate::DirectFromSource
                                : BoundaryGhostUpdate::Interpolation;

   this->readParticlesFiles();

   /* Seed the midpoint anchors for any pre-seeded buffer particles: the first
      midpointUpdateBuffers runs before the first sync, and would otherwise blend
      with uninitialized rho_old/v_old (TNL arrays do not construct fundamental
      types on resize). No-op for patches with zero particles. */
   for( long unsigned int p = 0; p < this->multiresolutionBoundaryPatches.size(); p++ )
      this->multiresolutionBoundaryPatches[ p ]->initializeIntegratorVariables();

   this->initBoundaryGhosts();

   log.writeSeparator();
   const bool hasMeasuretoolFile = params.template getParameter< std::string >( "measuretoolConfig" ) != "";
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
      std::string subdomainKey = "subdomain" + std::to_string( i );
      const int overlapWidth = getLevelOverlapWidth( i );

      this->fluidSets[ i ]->initializeAsDistributed(
         this->parametersSubdomains.template getParameter< int >( subdomainKey + "Fluid_n" ),
         this->parametersSubdomains.template getParameter< int >( subdomainKey + "Fluid_n_allocated" ),
         this->topology.getLocalGrid( i ),
         this->topology.getLocalOriginCoordinates( i ),
         this->topology.getGlobalGrid(),
         overlapWidth );
      if constexpr( ParticlesType::specifySearchedSetExplicitly() == true )
         this->fluidSets[ i ]->getParticles()->setParticleSetLabel( 0 );

      this->boundarySets[ i ]->initializeAsDistributed(
         this->parametersSubdomains.template getParameter< int >( subdomainKey + "Boundary_n" ),
         this->parametersSubdomains.template getParameter< int >( subdomainKey + "Boundary_n_allocated" ),
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
   this->multiresolutionBoundaryPatches.resize( mrbCount );

   int mrbIdx = 0;
   for( int i = 0; i < this->topology.getNumberOfSubdomains(); i++ ) {
      const float rf = getLevelRefinementFactor( i );

      for( const auto& iface : this->topology.getInterfacesOfSubdomain( i ) ) {
         const std::string bufferKey = "multiresolutionBuffer" + std::to_string( mrbIdx );
         const int mrbNumOfParticles = this->parametersSubdomains.template getParameter< int >( bufferKey + "N" );

         const int overlapWidth = ( rf < getLevelRefinementFactor( iface.neighborIdx ) ) ? 2 : 1;
         const RealType resolutionScale = overlapWidth * std::pow( rf, -int( SPHConfig::spaceDimension ) );
         const int mrbNumOfAllocatedParticles =
            std::max( int( 0.1 * this->getTotalFluidParticlesCount() * resolutionScale ), 2 * mrbNumOfParticles );

         this->multiresolutionBoundaryPatches[ mrbIdx ]->initializeAsDistributed(
            mrbNumOfParticles,
            mrbNumOfAllocatedParticles,
            this->topology.getLocalGrid( i ),
            this->topology.getLocalOriginCoordinates( i ),
            this->topology.getGlobalGrid(),
            overlapWidth );

         this->multiresolutionBoundaryPatches[ mrbIdx ]->initZones(
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

   this->refreshBoundaryGhostValues( subdomainIdx );

   this->model.interaction( this->fluidSets[ subdomainIdx ], this->boundarySets[ subdomainIdx ], this->modelParams );
   for( long unsigned int p = 0; p < this->multiresolutionBoundaryPatches.size(); p++ )
      if( this->multiresolutionBoundaryPatchInterfaces[ p ].ownIdx == subdomainIdx )
         this->model.interactionWithOpenBoundary(
            this->fluidSets[ subdomainIdx ], this->multiresolutionBoundaryPatches[ p ], this->modelParams );
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
   for( long unsigned int p = 0; p < this->multiresolutionBoundaryPatches.size(); p++ ) {
      const auto& iface = this->multiresolutionBoundaryPatchInterfaces[ p ];
      const int finerIdx = getLevelRefinementFactor( iface.ownIdx ) < getLevelRefinementFactor( iface.neighborIdx )
                            ? iface.ownIdx
                            : iface.neighborIdx;
      if( finerIdx != fineLevelIdx )
         continue;
      this->multiresolutionBoundaryPatches[ p ]->advectBufferParticles( dt );
   }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::refreshBuffersFromNeighbor( RealType theta, int fineLevelIdx )
{
   for( long unsigned int p = 0; p < this->multiresolutionBoundaryPatches.size(); p++ ) {
      const auto& iface = this->multiresolutionBoundaryPatchInterfaces[ p ];
      const int finerIdx = getLevelRefinementFactor( iface.ownIdx ) < getLevelRefinementFactor( iface.neighborIdx )
                            ? iface.ownIdx
                            : iface.neighborIdx;
      if( finerIdx != fineLevelIdx || iface.ownIdx != finerIdx )
         continue;
      this->fluidSets[ iface.neighborIdx ]->searchForNeighbors();
      this->multiresolutionBoundaryPatches[ p ]->refreshVariablesFromNeighbor(
         this->fluidSets[ iface.neighborIdx ], this->modelParams, theta );
   }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::stampBuffersAnchor( int fineLevelIdx )
{
   for( long unsigned int p = 0; p < this->multiresolutionBoundaryPatches.size(); p++ ) {
      const auto& iface = this->multiresolutionBoundaryPatchInterfaces[ p ];
      const int finerIdx = getLevelRefinementFactor( iface.ownIdx ) < getLevelRefinementFactor( iface.neighborIdx )
                            ? iface.ownIdx
                            : iface.neighborIdx;
      if( finerIdx != fineLevelIdx || iface.ownIdx != finerIdx )
         continue;
      this->multiresolutionBoundaryPatches[ p ]->initializeIntegratorVariables();
   }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::accumulateInterfaceMassesFromNeighbor(
   int neighborIdx,
   RealType dt )
{
   bool hasMatchingPatch = false;
   for( long unsigned int p = 0; p < this->multiresolutionBoundaryPatches.size(); p++ )
      if( this->multiresolutionBoundaryPatchInterfaces[ p ].neighborIdx == neighborIdx ) {
         hasMatchingPatch = true;
         break;
      }
   if( ! hasMatchingPatch )
      return;

   this->multiresolutionBoundaryPatches[ 0 ]->removeParticlesOutOfGrid( this->fluidSets[ neighborIdx ]->getParticles() );
   this->fluidSets[ neighborIdx ]->searchForNeighbors();
   for( long unsigned int p = 0; p < this->multiresolutionBoundaryPatches.size(); p++ )
      if( this->multiresolutionBoundaryPatchInterfaces[ p ].neighborIdx == neighborIdx )
         this->multiresolutionBoundaryPatches[ p ]->accumulateMasses( this->fluidSets[ neighborIdx ], this->modelParams, dt );
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::refreshSubstepSearches( int fineIdx )
{
   // positions of the fine fluid and of all buffer bands changed during the last substep
   if( ! this->multiresolutionBoundaryPatches.empty() )
      this->multiresolutionBoundaryPatches[ 0 ]->removeParticlesOutOfGrid( this->fluidSets[ fineIdx ]->getParticles() );
   this->fluidSets[ fineIdx ]->searchForNeighbors();

   for( long unsigned int p = 0; p < this->multiresolutionBoundaryPatches.size(); p++ ) {
      this->multiresolutionBoundaryPatches[ p ]->removeParticlesOutOfGrid(
         this->multiresolutionBoundaryPatches[ p ]->getParticles() );
      this->multiresolutionBoundaryPatches[ p ]->searchForNeighbors();
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

   for( long unsigned int p = 0; p < this->multiresolutionBoundaryPatches.size(); p++ ) {
      const auto& iface = this->multiresolutionBoundaryPatchInterfaces[ p ];
      this->multiresolutionBoundaryPatches[ p ]->syncUpdateInterfaceBuffer(
         this->fluidSets[ iface.ownIdx ],
         this->fluidSets[ iface.neighborIdx ],
         this->modelParams,
         this->timeStepping.getTimeStep() );
   }

   for( auto& fluid : this->fluidSets )
      fluid->searchForNeighbors();

   this->timeMeasurement.stop( "multiresolution-update" );
   this->writeLog( "Update multiresolution BC (LTS)...", "Done." );
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::writeProlog( bool writeSystemInformation ) noexcept
{
   BaseType::writeProlog( writeSystemInformation, false );

   auto& log = this->logger;
   log.writeHeader( "Local timesteps." );
   for( int i = 0; i < this->numberOfSubsets; i++ ) {
      log.writeParameter( "Time step (subdomain " + std::to_string( i ) + "):", levelTimeStepping[ i ].getTimeStep() );
      log.writeParameter( "Substeps per global step (subdomain " + std::to_string( i ) + "):", levelSubsteps[ i ] );
   }
   log.writeSeparator();

   this->writeMultiresolutionBoundaryPatchesProlog( "Multiresolution boundary buffer (LTS) " );
}

template< typename Model >
void
SolverMultiSetBlockMultiresolutionLocalTimeStepping< Model >::writeInfo( bool writeMultiresolutionBuffers ) noexcept
{
   BaseType::writeInfo( false );

   if( writeMultiresolutionBuffers ) {
      this->writeMultiresolutionBoundaryPatchesInfo( "Number of mr-buffer (LTS) " );
      this->logger.writeSeparator();
   }
}

}  // namespace TNL::SPH
