#include <TNL/Config/parseCommandLine.h>
#include <TNL/Config/parseINIConfigFile.h>
#include <TNL/Logger.h>

#include <SPH/configSetup.h>
#include <SPH/configInit.h>
#include "template/config.h"

#include <SPH/Models/WCSPH_DBC/control.h>

template< typename Simulation >
void
exec( Simulation& sph, TNL::Logger& log )
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

      // in case of variable time step, compute the step
      sph.computeTimeStep();

      // integrate
      sph.timeMeasurement.start( "integrate" );
      sph.integrator->integratStepVerlet( sph.fluid, sph.boundary, sph.timeStepping, SPHDefs::BCType::integrateInTime() );
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

      // compute energy of particles and output the energy data
      energyMonitor.computeEnergyLevels( sph.fluid, sph.modelParams );
      energyMonitor.output( sph.outputDirectory + "/energy.dat", sph.timeStepping.getStep(), sph.timeStepping.getTime() );

      // output particle data
      sph.makeSnapshot( log );

      //update time step
      sph.updateTime();
   }
}

int main( int argc, char* argv[] )
{
   // prepare client parameters
   TNL::Config::ParameterContainer cliParams;
   TNL::Config::ConfigDescription cliConfig;

   // prepare sph parameters
   TNL::Config::ParameterContainer parameters;
   TNL::Config::ConfigDescription config;

   try {
      TNL::SPH::template initialize< Simulation >( argc, argv, cliParams, cliConfig, parameters, config );
   }
   catch ( ... ) {
      return EXIT_FAILURE;
   }

   TNL::Logger log( 100, std::cout );
   Simulation sph;
   sph.init( parameters, log );

   // add custom timer measuremets
   sph.timeMeasurement.addTimer( "extrapolate-openbc" );
   sph.timeMeasurement.addTimer( "apply-openbc" );

   sph.writeProlog( log );
   exec( sph, log );
   sph.writeEpilog( log );
}

