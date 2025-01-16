#include <TNL/Config/parseCommandLine.h>
#include <TNL/Config/parseINIConfigFile.h>
#include <TNL/Logger.h>

#include <SPH/configSetup.h>
#include <SPH/configInit.h>
#include "template/config.h"

#include <SPH/Models/WCSPH_BI/control.h>

int main( int argc, char* argv[] )
{
   // get CLI parameters
   TNL::Config::ParameterContainer cliParams;
   TNL::Config::ConfigDescription cliConfig;

   // get config parameters
   TNL::Config::ParameterContainer parameters;
   TNL::Config::ConfigDescription config;

   try {
      TNL::SPH::template initialize< Simulation >( argc, argv, cliParams, cliConfig, parameters, config );
   }
   catch ( ... ) {
      return EXIT_FAILURE;
   }

   std::string logFileName = "results/simulation.log";
   std::ofstream logFile( logFileName );
   TNL::Logger log( 100, logFile );
   Simulation sph;
   sph.init( parameters, log );
   sph.writeProlog( log );

   sph.timeMeasurement.addTimer( "extrapolate-openbc" );
   sph.timeMeasurement.addTimer( "apply-openbc" );

   // Library model:
   EnergyFields energyMonitor;
   energyMonitor.init( sph.fluid, true );
   EnergyFieldsInletOutlet energyMonitorInletOutlet;
   energyMonitorInletOutlet.addInlet( sph.openBoundaryPatches[ 0 ] );
   energyMonitorInletOutlet.addOutlet( sph.openBoundaryPatches[ 1 ] );

   while( sph.timeStepping.runTheSimulation() )
   {
      // search for neighbros
      sph.timeMeasurement.start( "search" );
      sph.performNeighborSearch( log );
      sph.timeMeasurement.stop( "search" );
      sph.writeLog( log, "Search...", "Done." );

      // extrapolate open boundary
      sph.timeMeasurement.start( "extrapolate-openbc" );
      sph.extrapolateOpenBC();
      sph.timeMeasurement.stop( "extrapolate-openbc" );
      sph.writeLog( log, "Extrapolate open BC...", "Done." );

      // perform interaction with given model
      sph.timeMeasurement.start( "interact" );
      sph.interact();
      sph.timeMeasurement.stop( "interact" );
      sph.writeLog( log, "Interact...", "Done." );

      // compute and outpute energy levels
      energyMonitor.computeEnergyDerivatives( sph.fluid, sph.modelParams );
      energyMonitor.integrate( sph.timeStepping.getTimeStep() );

      //integrate
      sph.timeMeasurement.start( "integrate" );
      sph.integrator->integratStepVerlet( sph.fluid, sph.boundary, sph.timeStepping );
      sph.timeMeasurement.stop( "integrate" );
      sph.writeLog( log, "Integrate...", "Done." );

      // compute energy flow through open boudnaries
      energyMonitorInletOutlet.computeInflowEnergyLevels( sph.fluid, sph.openBoundaryPatches[ 0 ], sph.modelParams, sph.timeStepping.getTimeStep() );
      energyMonitorInletOutlet.computeOutletEnergyLevels( sph.fluid, sph.openBoundaryPatches[ 1 ], sph.modelParams, sph.timeStepping.getTimeStep() );
      energyMonitorInletOutlet.output( sph.outputDirectory + "/energyOpenBoundary.dat", sph.timeStepping.getStep(), sph.timeStepping.getTime() );

      // apply open boundary condition
      sph.timeMeasurement.start( "apply-openbc" );
      sph.applyOpenBC();
      sph.timeMeasurement.stop( "apply-openbc" );
      sph.writeLog( log, "Update open BC...", "Done." );

      // compute and outpute energy levels
      energyMonitor.computeEnergyLevels( sph.fluid, sph.modelParams );

      // output particle data
      sph.makeSnapshot( log );
      energyMonitor.output( sph.outputDirectory + "/energy.dat", sph.timeStepping.getStep(), sph.timeStepping.getTime() );

      //update time step
      sph.timeStepping.updateTimeStep();
   }

   sph.writeEpilog( log );
}

