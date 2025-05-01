#include "template/config.h"

template< typename Simulation >
void
exec( Simulation& sph )
{
   // add special tools
   EnergyFields energyMonitor;
   energyMonitor.init( sph.fluid, true );
   EnergyFieldsInletOutlet energyMonitorInletOutlet;
   energyMonitorInletOutlet.addInlet( sph.openBoundaryPatches[ 0 ] );
   energyMonitorInletOutlet.addOutlet( sph.openBoundaryPatches[ 1 ] );

   while( sph.timeStepping.runTheSimulation() )
   {
      // search for neighbros
      sph.performNeighborSearch();

      // extrapolate open boundary
      sph.extrapolateOpenBC();

      // perform interaction with given model
      sph.interact();

      // compute and outpute energy levels
      energyMonitor.computeEnergyDerivatives( sph.fluid, sph.modelParams );
      energyMonitor.integrate( sph.timeStepping.getTimeStep() );

      // in case of variable time step, compute the step
      sph.computeTimeStep();

      // integrate
      sph.integrateVerletStep();

      // compute energy flow through open boudnaries
      energyMonitorInletOutlet.computeInflowEnergyLevels( sph.fluid, sph.openBoundaryPatches[ 0 ], sph.modelParams, sph.timeStepping.getTimeStep() );
      energyMonitorInletOutlet.computeOutletEnergyLevels( sph.fluid, sph.openBoundaryPatches[ 1 ], sph.modelParams, sph.timeStepping.getTimeStep() );
      energyMonitorInletOutlet.output( sph.outputDirectory + "/energyOpenBoundary.dat", sph.timeStepping.getStep(), sph.timeStepping.getTime() );

      // apply open boundary condition
      sph.applyOpenBC();

      // compute energy of particles and output the energy data
      energyMonitor.computeEnergyLevels( sph.fluid, sph.modelParams );
      energyMonitor.output( sph.outputDirectory + "/energy.dat", sph.timeStepping.getStep(), sph.timeStepping.getTime() );

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

