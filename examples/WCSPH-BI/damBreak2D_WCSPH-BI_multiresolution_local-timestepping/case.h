#include "template/config.h"
#include <SPH/shared/removeParticlesOutOfDensityLimits.h>

#include <functional>
#include <stdexcept>

/*
 * Local-timestepping schedules (Documentation/multiresolution-local-timestepping.md),
 * selected at runtime by the lts-mode config key. Both variants share the same
 * Berger-Colella recursion over the levels: the coarsest level (index 0) takes one
 * step of dt_global; a finer level L takes levelSubsteps[L] steps of dt_L inside the
 * coarse step, recursing into L+1 after each of its own substeps so that the finer
 * level always advances within the just-integrated window of its parent. Buffer bands
 * drift ballistically during the substeps of the level that owns them (iface(i,i+1)
 * follows level i+1's clock); full buffer lifecycle runs only at the sync point, when
 * all level clocks coincide again.
 */

/* v1: a level's first substep per global step computes the active layer-1 pass of the
 * interface above it (full SPH + Verlet on layer 1 of the band); the remaining
 * substeps are fluid-only with both layers ballistic. */
template< typename Simulation >
requires std::is_same_v<
    typename Simulation::ModelParams::IntegrationScheme,
    TNL::SPH::IntegrationSchemes::VerletScheme< typename SPHDefs::SPHConfig >
>
void execV1( Simulation& sph )
{
   MassMonitor massMonitor;
   massMonitor.init( sph.fluidSets, sph.modelParams );

   const int nLevels = sph.numberOfSubsets;

   while( sph.timeStepping.runTheSimulation() )
   {
      // sweeps and neighbor lists at the synchronized state
      for( int i = 0; i < sph.numberOfSubsets; i++ )
         TNL::SPH::customFunctions::removeParticlesOutOfDensityLimits( sph.fluidSets[ i ], sph.modelParams );
      sph.removeParticlesOutOfDomain();
      sph.performNeighborSearch();

      // coarsest level: single step spanning the whole window
      sph.interactCoarseLevel( 0 );
      BoundaryCorrection::boundaryCorrection( sph.fluidSets[ 0 ], sph.boundarySets[ 0 ], sph.modelParams,
                                              sph.levelTimeStepping[ 0 ].getTimeStep() );
      sph.integrateLevel( 0 );
      sph.accumulateInterfaceMassesFromNeighbor( 0, sph.levelTimeStepping[ 0 ].getTimeStep() );

      std::vector< int > substepsDone( nLevels, 0 );
      std::function< void( int, int ) > advanceLevel = [ & ]( int level, int count )
      {
         for( int substep = 1; substep <= count; substep++ )
         {
            substepsDone[ level ]++;
            const bool firstSubstep = ( substepsDone[ level ] == 1 );
            const typename Simulation::RealType dtFine = sph.levelTimeStepping[ level ].getTimeStep();
            const bool eulerRestart = ( sph.levelTimeStepping[ level ].getStep() % 20 == 0 );

            if( ! firstSubstep )
               sph.refreshSubstepSearches( level );

            sph.interactFineSubstep( level, firstSubstep );
            BoundaryCorrection::boundaryCorrection( sph.fluidSets[ level ], sph.boundarySets[ level ], sph.modelParams, dtFine );
            sph.integrateLevel( level );
            if( firstSubstep )
               sph.integrateBuffersLayer1( dtFine, eulerRestart, level );
            sph.advectSubstepBuffers( dtFine, firstSubstep, level );
            sph.accumulateInterfaceMassesFromNeighbor( level, dtFine );

            if( level + 1 < nLevels )
               advanceLevel( level + 1, sph.levelSubsteps[ level + 1 ] / sph.levelSubsteps[ level ] );
         }
      };
      advanceLevel( 1, sph.levelSubsteps[ 1 ] );

      // sync-point buffer lifecycle, then output and bookkeeping at the coarse clock
      sph.syncMultiresolutionUpdate();
      sph.makeSnapshot();

      massMonitor.sumTotalMass( sph.fluidSets );
      massMonitor.output( sph.outputDirectory + "/massConservation.dat", sph.timeStepping.getStep(), sph.timeStepping.getTime() );

      sph.measure();

      sph.timeStepping.updateTimeStep();
   }
}

/* v1.5 (midpoint): no active band pass at any substep - both band layers drift
 * ballistically. Instead, the fine-side band state is refreshed once per global step
 * at the global midpoint (after level 1's first substep has fully completed,
 * including the substeps of all deeper levels) as a theta-blend of the sync-time
 * anchor kept in rho_old/v_old with a fresh interpolation from the already advanced
 * coarser neighbor. */
template< typename Simulation >
requires std::is_same_v<
    typename Simulation::ModelParams::IntegrationScheme,
    TNL::SPH::IntegrationSchemes::VerletScheme< typename SPHDefs::SPHConfig >
>
void execMidpoint( Simulation& sph, const typename Simulation::RealType theta )
{
   MassMonitor massMonitor;
   massMonitor.init( sph.fluidSets, sph.modelParams );

   const int nLevels = sph.numberOfSubsets;

   while( sph.timeStepping.runTheSimulation() )
   {
      // sweeps and neighbor lists at the synchronized state
      for( int i = 0; i < sph.numberOfSubsets; i++ )
         TNL::SPH::customFunctions::removeParticlesOutOfDensityLimits( sph.fluidSets[ i ], sph.modelParams );
      sph.removeParticlesOutOfDomain();
      sph.performNeighborSearch();

      // coarsest level: single step spanning the whole window
      sph.interactCoarseLevel( 0 );
      BoundaryCorrection::boundaryCorrection( sph.fluidSets[ 0 ], sph.boundarySets[ 0 ], sph.modelParams,
                                              sph.levelTimeStepping[ 0 ].getTimeStep() );
      sph.integrateLevel( 0 );
      sph.accumulateInterfaceMassesFromNeighbor( 0, sph.levelTimeStepping[ 0 ].getTimeStep() );

      std::vector< int > substepsDone( nLevels, 0 );
      std::function< void( int, int ) > advanceLevel = [ & ]( int level, int count )
      {
         for( int substep = 1; substep <= count; substep++ )
         {
            substepsDone[ level ]++;
            const bool firstSubstep = ( substepsDone[ level ] == 1 );
            const typename Simulation::RealType dtFine = sph.levelTimeStepping[ level ].getTimeStep();

            if( ! firstSubstep )
               sph.refreshSubstepSearches( level );

            sph.interactFineSubstep( level, false );
            BoundaryCorrection::boundaryCorrection( sph.fluidSets[ level ], sph.boundarySets[ level ], sph.modelParams, dtFine );
            sph.integrateLevel( level );
            sph.advectSubstepBuffers( dtFine, false, level );
            sph.accumulateInterfaceMassesFromNeighbor( level, dtFine );

            if( level + 1 < nLevels )
               advanceLevel( level + 1, sph.levelSubsteps[ level + 1 ] / sph.levelSubsteps[ level ] );

            if( level == 1 && firstSubstep )
               for( int fineIdx = 1; fineIdx < nLevels; fineIdx++ )
                  sph.midpointUpdateBuffers( theta, fineIdx );
         }
      };
      advanceLevel( 1, sph.levelSubsteps[ 1 ] );

      // sync-point buffer lifecycle, then output and bookkeeping at the coarse clock
      sph.syncMultiresolutionUpdate();
      sph.makeSnapshot();

      massMonitor.sumTotalMass( sph.fluidSets );
      massMonitor.output( sph.outputDirectory + "/massConservation.dat", sph.timeStepping.getStep(), sph.timeStepping.getTime() );

      sph.measure();

      sph.timeStepping.updateTimeStep();
   }
}

int main( int argc, char* argv[] )
{
   Simulation sph;
   sph.init( argc, argv );

   const std::string ltsMode = sph.parameters.checkParameter( "lts-mode" )
                                 ? sph.parameters.getParameter< std::string >( "lts-mode" ) : "v1";
   const double ltsTheta = sph.parameters.checkParameter( "lts-theta" )
                             ? sph.parameters.getParameter< double >( "lts-theta" ) : 0.5;

   sph.writeProlog();
   sph.logger.writeParameter( "LTS schedule:", ltsMode );
   if( ltsMode != "v1" )
      sph.logger.writeParameter( "LTS midpoint theta:", ltsTheta );

   if( ltsMode == "v1" )
      execV1( sph );
   else if( ltsMode == "v1.5" || ltsMode == "midpoint" )
      execMidpoint( sph, static_cast< Simulation::RealType >( ltsTheta ) );
   else
      throw std::invalid_argument( "Unknown lts-mode '" + ltsMode + "'." );

   sph.writeEpilog();
}

