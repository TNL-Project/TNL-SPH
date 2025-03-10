#include <TNL/Config/parseCommandLine.h>
#include <TNL/Config/parseINIConfigFile.h>
#include <TNL/Logger.h>

#include <SPH/configSetup.h>
#include <SPH/configInit.h>
#include "template/config.h"

#include <SPH/Models/WCSPH_DBC/control.h>
#include <cstdlib>

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
   sph.writeProlog( log );

   // Solver model:

   //sph.init( parameters );
   //sph.writeProlog( parameters );
   //sph.exec();
   //sph.writeEpilog( parameters );

   // Library model:

   /*
   while( sph.timeStepping.runTheSimulation() )
   {
      // search for neighbros
      sph.timeMeasurement.start( "search" );
      sph.performNeighborSearch( log );
      sph.timeMeasurement.stop( "search" );
      sph.writeLog( log, "Search...", "Done." );

      // perform interaction with given model
      sph.timeMeasurement.start( "interact" );
      sph.interact();
      sph.timeMeasurement.stop( "interact" );
      sph.writeLog( log, "Interact...", "Done." );

      // in case of variable time step, compute the step
      sph.computeTimeStep();

      //integrate
      sph.timeMeasurement.start( "integrate" );
      sph.integrator->integratStepVerlet( sph.fluidSets[ 0 ], sph.boundarySets[ 0 ], sph.timeStepping, SPHDefs::BCType::integrateInTime() );
      sph.timeMeasurement.stop( "integrate" );
      sph.writeLog( log, "Integrate...", "Done." );

      // output particle data
      sph.makeSnapshot( log );
      // check timers and if measurement or interpolation should be performed, is performed
      sph.template measure< SPHDefs::KernelFunction, SPHDefs::EOS >( log );

      // update time step
      sph.updateTime();
   }
   */

   sph.writeEpilog( log );
}

