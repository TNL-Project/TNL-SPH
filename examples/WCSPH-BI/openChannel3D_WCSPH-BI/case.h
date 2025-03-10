#include <TNL/Config/parseCommandLine.h>
#include <TNL/Config/parseINIConfigFile.h>
#include <TNL/Logger.h>

#include <SPH/configSetup.h>
#include <SPH/configInit.h>
#include "template/config.h"

#include <SPH/Models/WCSPH_BI/control.h>

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

   //std::string logFileName = "results/simulation.log";
   //std::ofstream logFile( logFileName );
   //TNL::Logger log( 100, logFile );
   TNL::Logger log( 100, std::cout );
   Simulation sph;
   sph.init( parameters, log );
   sph.writeProlog( log );

   // Solver model:

   //sph.init( parameters );
   //sph.writeProlog( parameters );
   //sph.exec();
   //sph.writeEpilog( parameters );
   //sph.exec();

   sph.timeMeasurement.addTimer( "extrapolate-openbc" );
   sph.timeMeasurement.addTimer( "apply-openbc" );


   // Library model:

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

      // custom: no-penetration bc
      BoundaryCorrection::boundaryCorrection( sph.fluid, sph.boundary, sph.modelParams, sph.timeStepping.getTimeStep() );

      // compute new time step
      sph.computeTimeStep();

      //integrate
      sph.timeMeasurement.start( "integrate" );
      sph.integrator->integratStepVerlet( sph.fluid, sph.boundary, sph.timeStepping );
      sph.timeMeasurement.stop( "integrate" );
      sph.writeLog( log, "Integrate...", "Done." );

      // output particle data
      sph.makeSnapshot( log );

      // apply open boundary condition
      sph.timeMeasurement.start( "apply-openbc" );
      sph.applyOpenBC();
      sph.timeMeasurement.stop( "apply-openbc" );
      sph.writeLog( log, "Update open BC...", "Done." );

      //update time step
      sph.timeStepping.updateTimeStep();
   }

   sph.writeEpilog( log );
}

