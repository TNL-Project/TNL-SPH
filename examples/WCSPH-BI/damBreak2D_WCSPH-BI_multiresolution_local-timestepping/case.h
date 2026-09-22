#include "template/config.h"
#include <SPH/shared/removeParticlesOutOfDensityLimits.h>

#include <functional>

/*
 * Local-timestepping schedule (Documentation/multiresolution-local-timestepping.md):
 * Berger-Colella recursion over the levels. The coarsest level (index 0) takes one
 * step of dt_global; a finer level L takes levelSubsteps[L] steps of dt_L inside the
 * coarse step, recursing into L+1 after each of its own substeps so that the finer
 * level always advances within the just-integrated window of its parent. Buffer bands
 * drift ballistically during the substeps of the level that owns them (iface(i,i+1)
 * follows level i+1's clock); full buffer lifecycle runs only at the sync point, when
 * all level clocks coincide again.
 *
 * A level's first substep per global step computes the active layer-1 pass of the
 * interface above it (full SPH + Verlet on the band); remaining substeps are fluid
 * only with both layers ballistic — the v1 "band state fixed at sync point" model.
 */
template< typename Simulation >
requires std::is_same_v<
    typename Simulation::ModelParams::IntegrationScheme,
    TNL::SPH::IntegrationSchemes::VerletScheme< typename SPHDefs::SPHConfig >
>
void exec( Simulation& sph )
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

int main( int argc, char* argv[] )
{
   Simulation sph;
   sph.init( argc, argv );
   sph.writeProlog();
   exec( sph );
   sph.writeEpilog();
}
