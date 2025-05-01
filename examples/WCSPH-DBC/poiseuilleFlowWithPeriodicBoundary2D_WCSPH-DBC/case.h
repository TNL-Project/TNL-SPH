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

      // enforce open boundary particles
      sph.applyPeriodicBCEnforce();

      // perform interaction with given model
      sph.interact();

      //integrate
      sph.integrateVerletStep();

      // transform particles inside periodic boundary conditions
      sph.applyPeriodicBCTransfer();

      // check timers and if output should be performed, it is performed
      sph.makeSnapshot();

      //update time step
      sph.timeStepping.updateTimeStep();
   }

   sph.writeEpilog();
}

