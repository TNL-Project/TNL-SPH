#include "template/config.h"

template< typename Simulation >
void
exec( Simulation& sph )
{
   while( sph.timeStepping.runTheSimulation() )
   {
      // search for neighbros
      sph.performNeighborSearch();

      // enforce open boundary particles
      sph.applyPeriodicBCEnforce();

      // perform interaction with given model
      sph.interact();

      // in case variable time stepping is used, compute time step
      sph.computeTimeStep();

      //integrate
      sph.integrateVerletStep();

      // transform particles inside periodic boundary conditions
      sph.applyPeriodicBCTransfer();

      // output particle data
      sph.makeSnapshot();

      //update time step
      sph.updateTime();
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

