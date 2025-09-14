#include "template/config.h"

int main( int argc, char* argv[] )
{
   Simulation sph;
   sph.init( argc, argv );
   sph.writeProlog();

   // Initialize energy monitor:
   EnergyFields energyMonitor;
   energyMonitor.init( sph.fluid, true );
   EnergyFieldsInletOutlet energyMonitorInletOutlet;
   energyMonitorInletOutlet.addInlet( sph.openBoundaryPatches[ 0 ] );
   energyMonitorInletOutlet.addOutlet( sph.openBoundaryPatches[ 1 ] );

   // Initialize flow rate:
   FlowRateMonitor flowRateMonitor;
   flowRateMonitor.init( sph.fluid, { 0.25f, 0.f }, { 0.25f, 0.1f }, { 1.f, 0.f } );


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

      // compute and outpute energy levels
      energyMonitor.computeEnergyDerivatives( sph.fluid, sph.modelParams );
      energyMonitor.integrate( sph.timeStepping.getTimeStep() );

      // measure volumetric flow rate
      flowRateMonitor.measureVolumetricFlowRate( sph.fluid, sph.modelParams, sph.timeStepping.getTimeStep() );
      flowRateMonitor.output( sph.outputDirectory + "/volumetricFlowRate.dat", sph.timeStepping.getStep(), sph.timeStepping.getTime() );

      //integrate
      sph.integrateVerletStep();

      // compute energy flow through open boudnaries
      energyMonitorInletOutlet.computeInflowEnergyLevels( sph.fluid, sph.openBoundaryPatches[ 0 ], sph.modelParams, sph.timeStepping.getTimeStep() );
      energyMonitorInletOutlet.computeOutletEnergyLevels( sph.fluid, sph.openBoundaryPatches[ 1 ], sph.modelParams, sph.timeStepping.getTimeStep() );
      energyMonitorInletOutlet.output( sph.outputDirectory + "/energyOpenBoundary.dat", sph.timeStepping.getStep(), sph.timeStepping.getTime() );

      // apply open boundary condition
      sph.applyOpenBC();

      // compute and outpute energy levels
      energyMonitor.computeEnergyLevels( sph.fluid, sph.modelParams );
      energyMonitor.output( sph.outputDirectory + "/energy.dat", sph.timeStepping.getStep(), sph.timeStepping.getTime() );

      // output particle data
      sph.makeSnapshot();

      //update time step
      sph.timeStepping.updateTimeStep();
      sph.updateTime();
   }

   sph.writeEpilog();
}

