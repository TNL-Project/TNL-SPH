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

      // extrapolate open boundary
      sph.extrapolateOpenBC();

      // perform interaction with given model
      sph.interact();

      // make integration step with Verlet scheme
      sph.integrateVerletStep();

      // apply open boundary condition
      sph.applyOpenBC();

      // check timers and if output should be performed, it is performed
      sph.makeSnapshot();

      // update time step
      sph.updateTime();
   }

   sph.writeEpilog();
}

