#include "template/config.h"
#include <SPH/shared/removeParticlesOutOfDensityLimits.h>

#include <functional>

/*
 * Local-timestepping schedule (v2, unified):
 *
 * Berger-Colella recursion over the refinement levels with the coarsest level
 * included. Level L takes getNumberOfSubsteps(L) substeps of dt_L per global window
 * and recurses into L+1 after each of them, so the finer level always advances
 * within the just-integrated window of its parent. Substeps are indexed from zero;
 * the window-global index of a substep is g = base + substep, where base is threaded
 * through the recursion (childBase = g * ratio). The genuine first substep of every
 * level (g == 0) skips the search refresh - the sync-time performNeighborSearch()
 * covers all levels; every later substep calls refreshSubstepSearches().
 *
 * The buffer band of iface(i,i+1) follows the clock of level i+1: it advects
 * ballistically during that level's substeps and its state (rho, v) is re-grounded
 * from the coarser neighbor inside the window by refreshBuffersFromNeighbor():
 *   (R1) exact (theta = 1) whenever a parent substep ends and the child invocation
 *        has just landed on the same physical time (level-clock coincidence);
 *        stampBuffersAnchor() then re-stamps the band anchors, so the following
 *        blends stay time-bracketed around the band's next substep interval,
 *   (R2) blended q <- theta * q_interp(neighbor now) + (1 - theta) * q_anchor, at
 *        the midpoint substep of each level's invocation, for the band owned by that
 *        level (its coarser neighbor may have no interior coincidence, e.g.
 *        iface(0,1) with the single-step coarsest level).
 * With theta = 0.5 every blend is the exact linear interpolation at the band's own
 * local time (coincidences and syncs refresh exactly). In a two-level run, no
 * interior coincidence exists and the schedule reduces to the previous
 * single-trigger midpoint blend.
 *
 * Membership of the buffer bands (creation, conversion, removal of particles) is
 * frozen during the window and re-evaluated only at the sync point by
 * syncMultiresolutionUpdate(), when all level clocks coincide again.
 */

template< typename Simulation >
requires std::is_same_v<
    typename Simulation::ModelParams::IntegrationScheme,
    TNL::SPH::IntegrationSchemes::VerletScheme< typename SPHDefs::SPHConfig >
>
void execLocalTimestepping( Simulation& sph )
{
   MassMonitor massMonitor;
   massMonitor.init( sph.fluidSets, sph.modelParams );

   const typename Simulation::RealType theta = sph.modelParams.ltsTheta;
   const int nLevels = sph.numberOfSubsets;

   while( sph.timeStepping.runTheSimulation() )
   {
      // sweeps and neighbor lists at the synchronized state
      for( int i = 0; i < sph.numberOfSubsets; i++ )
         TNL::SPH::customFunctions::removeParticlesOutOfDensityLimits( sph.fluidSets[ i ], sph.modelParams );
      sph.removeParticlesOutOfDomain();
      sph.performNeighborSearch();

      // boundary ghosts are part of the sync-time boundary state, frozen for the whole window
      sph.refreshBoundaryGhostValues();

      // Berger-Colella recursion with the coarsest level included
      std::function< void( int, int, int ) > advanceLevel = [ & ]( int level, int count, int base )
      {
         for( int substep = 0; substep < count; substep++ )
         {
            const int g = base + substep;  // window-global substep index of this level
            const typename Simulation::RealType dtLevel = sph.levelTimeStepping[ level ].getTimeStep();

            if( g != 0 )
               sph.refreshSubstepSearches( level );

            sph.interact( level );
            BoundaryCorrection::boundaryCorrection( sph.fluidSets[ level ], sph.boundarySets[ level ], sph.modelParams, dtLevel );
            sph.integrateLevel( level );
            sph.advectSubstepBuffers( dtLevel, level );
            sph.accumulateInterfaceMassesFromNeighbor( level, dtLevel );

            if( level + 1 < nLevels ) {
               const int ratio = sph.getNumberOfSubsteps( level + 1 ) / sph.getNumberOfSubsteps( level );
               advanceLevel( level + 1, ratio, g * ratio );

               /* R1 - exact coincidence refresh: the child invocation just brought the
                  band's clock to the same physical time as this level; re-stamp the
                  anchors so the following blends stay bracketed. Skipped after the
                  last substep - the sync-point lifecycle re-grounds all bands anyway. */
               if( substep + 1 < count ) {
                  sph.refreshBuffersFromNeighbor( 1.0, level + 1 );
                  sph.stampBuffersAnchor( level + 1 );
               }
            }

            /* R2 - midpoint blend of the band owned by this level (the interface
               toward its coarser neighbor), firing once per invocation at the
               midpoint substep. */
            if( theta > 0 && level > 0 && substep == count / 2 - 1 )
               sph.refreshBuffersFromNeighbor( theta, level );
         }
      };
      advanceLevel( 0, 1, 0 );

      // sync-point buffer lifecycle, then output and bookkeeping at the coarse clock
      sph.makeSnapshot();
      sph.syncMultiresolutionUpdate();

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
   execLocalTimestepping( sph );
   sph.writeEpilog();
}

