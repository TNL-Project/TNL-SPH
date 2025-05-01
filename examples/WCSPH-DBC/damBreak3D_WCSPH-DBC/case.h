#include "template/config.h"

int main( int argc, char* argv[] )
{
   Simulation sph;
   sph.init( argc, argv );
   sph.writeProlog();

   while( sph.timeStepping.runTheSimulation() )
   {
      // search for neighbros
      sph.performNeighborSearch();

      // perform interaction with given model
      sph.interact();

      // in case of variable time step, compute the step
      sph.computeTimeStep();

      //integrate
      sph.timeMeasurement.start( "integrate" );
      sph.integrator->integratStepVerlet( sph.fluid, sph.boundary, sph.timeStepping, SPHDefs::BCType::integrateInTime() );
      sph.timeMeasurement.stop( "integrate" );
      sph.writeLog( "Integrate...", "Done." );

      // check timers and if output should be performed, it is performed
      sph.makeSnapshot();

      // check timers and if measurement or interpolation should be performed, it is performed
      sph.measure();

      // update time step
      sph.updateTime();
   }

   sph.writeEpilog();
}

