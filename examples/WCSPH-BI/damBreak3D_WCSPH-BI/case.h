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
      // custom: no-penetration bc
      BoundaryCorrection::boundaryCorrection( sph.fluid, sph.boundary, sph.modelParams, sph.timeStepping.getTimeStep() );

      //integrate
      sph.integrateVerletStep();

      // output particle data
      sph.makeSnapshot();

      // check timers and if measurement or interpolation should be performed, is performed
      sph.measure();

      // update time step
      sph.timeStepping.updateTimeStep();
   }

   sph.writeEpilog();
}

