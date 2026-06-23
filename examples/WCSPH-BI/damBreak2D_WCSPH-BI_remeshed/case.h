#include "template/config.h"
#include <SPH/shared/removeParticlesOutOfDensityLimits.h>

template< typename Simulation >
void exec( Simulation& sph )
{
   MassMonitor massMonitor;
   massMonitor.init( sph.fluidSets, sph.modelParams );

   while( sph.timeStepping.runTheSimulation() )
   {
      // search for neighbors
      sph.performNeighborSearch();

      // perform interaction with given model
      sph.interact();
      // custom: no-penetration bc
      BoundaryCorrection::boundaryCorrection( sph.fluidSets[ 0 ], sph.boundarySets[ 0 ], sph.modelParams, sph.timeStepping.getTimeStep() );

      // integrate
      sph.integrate();

      // output particle data
      sph.makeSnapshot();

      // remesh particles
      sph.remeshParticles();

      // CUSTOM: check mass conservation
      massMonitor.sumTotalMass( sph.fluidSets );
      massMonitor.output( sph.outputDirectory + "/massConservation.dat", sph.timeStepping.getStep(), sph.timeStepping.getTime() );

      // check timers and if measurement or interpolation should be performed, is performed
      sph.measure();

      // update time step
      sph.timeStepping.updateTimeStep();
   }
}

int main( int argc, char* argv[] )
{
   Simulation sph;
   sph.initRemeshedSimulation( argc, argv );
   sph.writeProlog();
   exec( sph );
   sph.writeEpilog();
}

