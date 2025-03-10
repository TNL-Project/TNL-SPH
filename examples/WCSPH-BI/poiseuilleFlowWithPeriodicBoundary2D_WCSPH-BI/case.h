#include <TNL/Config/parseCommandLine.h>
#include <TNL/Config/parseINIConfigFile.h>
#include <TNL/Logger.h>

#include <SPH/configSetup.h>
#include <SPH/configInit.h>
#include "template/config.h"

#include <SPH/Models/WCSPH_BI/control.h>

template< typename Simulation >
void
exec( Simulation& sph, TNL::Logger& log )
{
   while( sph.timeStepping.runTheSimulation() )
   {
      // search for neighbros
      sph.timeMeasurement.start( "search" );
      sph.performNeighborSearch( log );
      sph.timeMeasurement.stop( "search" );
      sph.writeLog( log, "Search...", "Done." );

      // enforce open boundary particles
      sph.timeMeasurement.start( "enforce-periodic-bc" );
      sph.applyPeriodicBCEnforce();
      sph.timeMeasurement.stop( "enforce-periodic-bc" );
      sph.writeLog( log, "Apply periodic BC...", "Done." );

      // perform interaction with given model
      sph.timeMeasurement.start( "interact" );
      sph.interact();
      sph.timeMeasurement.stop( "interact" );
      sph.writeLog( log, "Interact...", "Done." );

      // in case variable time stepping is used, compute time step
      sph.computeTimeStep();

      //integrate
      sph.timeMeasurement.start( "integrate" );
      sph.integrator->integratStepVerlet( sph.fluid, sph.boundary, sph.timeStepping );
      sph.timeMeasurement.stop( "integrate" );
      sph.writeLog( log, "Integrate...", "Done." );

      // transform particles inside periodic boundary conditions
      sph.timeMeasurement.start( "transfer-periodic-bc" );
      sph.applyPeriodicBCTransfer();
      sph.timeMeasurement.stop( "transfer-periodic-bc" );
      sph.writeLog( log, "Transfer periodic BC...", "Done." );

      // output particle data
      sph.makeSnapshot( log );

      //update time step
      sph.updateTime();
   }
}

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

   TNL::Logger log( 100, std::cout );
   Simulation sph;
   sph.init( parameters, log );

   // add custom timers
   sph.timeMeasurement.addTimer( "enforce-periodic-bc" );
   sph.timeMeasurement.addTimer( "transfer-periodic-bc" );

   sph.writeProlog( log );
   exec( sph, log );
   sph.writeEpilog( log );
}

