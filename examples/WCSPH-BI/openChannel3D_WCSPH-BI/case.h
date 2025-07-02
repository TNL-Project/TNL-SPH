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
      // custom: no-penetration bc
      BoundaryCorrection::boundaryCorrection( sph.fluid, sph.boundary, sph.modelParams, sph.timeStepping.getTimeStep() );

      // compute new time step
      sph.computeTimeStep();

      //integrate
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

