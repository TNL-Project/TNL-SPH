#include "template/config.h"
#include <SPH/shared/removeParticlesOutOfDensityLimits.h>

/*
 * Local-timestepping driver (schedule v2). The Berger-Colella recursion over the
 * levels, the search refreshes and the band state refreshes (exact at level-clock
 * coincidences, theta-blended at level midpoints) are driven by the solver - see
 * SolverMultiSetBlockMultiresolutionLocalTimestepping::recursiveSubdomainsLoop. The
 * case supplies only the per-substep physics via the lambda below; band membership
 * (creation, conversion, removal) is re-evaluated only at the sync point.
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

   while( sph.timeStepping.runTheSimulation() )
   {
      // sweeps and neighbor lists at the synchronized state
      for( int i = 0; i < sph.numberOfSubsets; i++ )
         TNL::SPH::customFunctions::removeParticlesOutOfDensityLimits( sph.fluidSets[ i ], sph.modelParams );
      sph.removeParticlesOutOfDomain();
      sph.performNeighborSearch();

      // boundary ghosts are part of the sync-time boundary state, frozen for the whole window
      sph.refreshBoundaryGhostValues();

      // substep physics of one level - the local-timestepping schedule drives it
      sph.recursiveSubdomainsLoop(
         [ & ]( int level, int /*substep*/, int /*globalSubstep*/, typename Simulation::RealType dt )
      {
         sph.interact( level );
         BoundaryCorrection::boundaryCorrection( sph.fluidSets[ level ], sph.boundarySets[ level ], sph.modelParams, dt );
         sph.integrateLevel( level );
         sph.advectSubstepBuffers( dt, level );
         sph.accumulateInterfaceMassesFromNeighbor( level, dt );
      } );

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
   sph.logger.writeParameter( "LTS schedule:", "v2 (unified)" );

   execLocalTimestepping( sph );

   sph.writeEpilog();
}
