#include <TNL/Config/parseCommandLine.h>
#include <TNL/Config/parseINIConfigFile.h>
#include <TNL/Logger.h>

#include <SPH/configSetup.h>
#include <SPH/configInit.h>
#include "template/config.h"

#include <SPH/Models/WCSPH_BI/control.h>
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

   //while( sph.timeStepping.getStep() < 1 )
   while( sph.timeStepping.runTheSimulation() )
   {
      std::cout << "Step: " << sph.timeStepping.getStep() << std::endl;

      // search for neighbros
      sph.timeMeasurement.start( "search" );
      sph.performNeighborSearch( log, true );
      sph.timeMeasurement.stop( "search" );
      sph.writeLog( log, "Search...", "Done." );

      //std::cout << sph.fluidSets[0]->getParticles()->getCellFirstLastParticleList() << std::endl;
      // perform interaction with given model
      sph.timeMeasurement.start( "interact" );
      sph.interact();
      BoundaryCorrection::boundaryCorrection( sph.fluidSets[ 0 ], sph.boundarySets[ 0 ], sph.modelParams, sph.timeStepping.getTimeStep() );
      BoundaryCorrection::boundaryCorrection( sph.fluidSets[ 1 ], sph.boundarySets[ 1 ], sph.modelParams, sph.timeStepping.getTimeStep() );
      sph.timeMeasurement.stop( "interact" );
      sph.writeLog( log, "Interact...", "Done." );

      // in case of variable time step, compute the step
      sph.computeTimeStep();

      //integrate
      sph.timeMeasurement.start( "integrate" );
      sph.integrator->integratStepVerlet( sph.fluidSets[ 0 ], sph.boundarySets[ 0 ], sph.timeStepping, false );
      sph.integrator->integratStepVerlet( sph.fluidSets[ 1 ], sph.boundarySets[ 1 ], sph.timeStepping, false );
      sph.timeMeasurement.stop( "integrate" );
      sph.writeLog( log, "Integrate...", "Done." );

      // output particle data
      sph.makeSnapshot( log );
      // check timers and if measurement or interpolation should be performed, is performed
      //sph.template measure< SPHDefs::KernelFunction, SPHDefs::EOS >( log );

      // update MR buffer
      sph.applyMultiresolutionBC();

      //NOTE: Info
      if( sph.timeStepping.getStep() > 10782  )
         sph.writeInfo( log );


      // update time step
      sph.updateTime();
   }

   sph.writeEpilog( log );
}

/*

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

      // in case of variable time step, compute the step
      sph.computeTimeStep();

      // make integration step with Verlet scheme
      sph.integrateVerletStep( SPHDefs::BCType::integrateInTime() );

      // check timers and if output should be performed, it is performed
      sph.makeSnapshot();

      // check timers and if measurement or interpolation should be performed, it is performed
      sph.measure();

      // update time step
      sph.updateTime();
   }

   sph.writeEpilog();
}

*/

