#include <TNL/Config/parseCommandLine.h>
#include <TNL/Config/parseINIConfigFile.h>
#include <TNL/Logger.h>

#include <SPH/configSetup.h>
#include <SPH/configInit.h>
#include "template/config.h"

#include <SPH/Models/WCSPH_DBC/control.h>
#include <string>

int main( int argc, char* argv[] )
{

   TNL::MPI::ScopedInitializer mpi( argc, argv );

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

   //DEBUG:
   std::string logFileName = "results/simulationLog_rank" + std::to_string( TNL::MPI::GetRank() );
   std::ofstream logFile( logFileName );
   //TNL::Logger log( 100, std::cout );
   TNL::Logger log( 100, logFile );

   Simulation sph;
   TNL::MPI::Barrier( sph.communicator ); //To have clear output

   sph.init( parameters, log );
   TNL::MPI::Barrier( sph.communicator ); //To have clear output

   sph.writeProlog( log );

   // Solver model:

   //sph.init( parameters );
   //sph.writeProlog( parameters );
   //sph.exec();
   //sph.writeEpilog( parameters );
   sph.timeMeasurement.addTimer( "synchronize" );
   sph.timeMeasurement.addTimer( "rebalance" );

   while( sph.timeStepping.runTheSimulation() )
   //while( sph.timeStepping.getStep() < 2 )
   {
      if( ( sph.timeStepping.getStep() % 50  == 0 )  )
         sph.writeInfo( log );

      // SingleSet: forgot overlaps
      //sph.resetOverlaps();
      //sph.writeLog( log, "Reset overlaps...", "Done." );

      TNL::MPI::Barrier( sph.communicator ); //To have clear output

      // search for neighbros
      sph.timeMeasurement.start( "search" );
      sph.performNeighborSearch( log , true );
      sph.timeMeasurement.stop( "search" );
      sph.writeLog( log, "Search...", "Done." );

      TNL::MPI::Barrier( sph.communicator ); //To have clear output

      sph.writeLog( log, "Starting synchronization.", "" );
      sph.timeMeasurement.start( "synchronize" );
      sph.synchronizeDistributedSimulation( log );
      sph.timeMeasurement.stop( "synchronize" );
      sph.writeLog( log, "Synchronize...", "Done." );

      TNL::MPI::Barrier( sph.communicator ); //To have clear output

      if( ( sph.timeStepping.getStep() % 100  == 0 ) && ( sph.timeStepping.getStep() > 1 ) ){

         sph.timeMeasurement.start( "search" );
         sph.performNeighborSearch( log, true );
         sph.timeMeasurement.stop( "search" );
         sph.writeLog( log, "Search second...", "Done." );

         //TNL::MPI::Barrier( sph.communicator ); //To have clear output

         log.writeSeparator();
         sph.writeLog( log, "Starting load balancing.", "" );
         sph.timeMeasurement.start( "rebalance" );
         sph.performLoadBalancing( log );
         sph.timeMeasurement.stop( "rebalance" );
         sph.writeLog( log, "Load balancing...", "Done." );
         log.writeSeparator();

         TNL::MPI::Barrier( sph.communicator ); //To have clear output

         // SingleSet: forgot overlaps
         sph.resetOverlaps();
         sph.writeLog( log, "Reset overlaps...", "Done." );

         TNL::MPI::Barrier( sph.communicator ); //To have clear output


         sph.timeMeasurement.start( "search" );
         sph.performNeighborSearch( log, true );
         sph.timeMeasurement.stop( "search" );
         sph.writeLog( log, "Search second...", "Done." );

         TNL::MPI::Barrier( sph.communicator ); //To have clear output

         sph.writeLog( log, "Starting synchronization.", "" );
         sph.timeMeasurement.start( "synchronize" );
         sph.synchronizeDistributedSimulation( log );
         sph.timeMeasurement.stop( "synchronize" );
         sph.writeLog( log, "Synchronize...", "Done." );
         TNL::MPI::Barrier( sph.communicator ); //To have clear output
      }

      // SingleSet: perform second search to assign received particles
      sph.timeMeasurement.start( "search" );
      sph.performNeighborSearch( log, true );
      sph.timeMeasurement.stop( "search" );
      sph.writeLog( log, "Search second...", "Done." );

      TNL::MPI::Barrier( sph.communicator ); //To have clear output

      // output particle data
      if( sph.timeStepping.checkOutputTimer( "save_results" ) )
      {
         /**
          * Compute pressure from density.
          * This is not necessary since we do this localy, if pressure is needed.
          * It's useful for output anyway.
          */
         sph.model.computePressureFromDensity( sph.fluid, sph.modelParams );
         sph.model.computePressureFromDensity( sph.boundary, sph.modelParams );

         sph.save( log, true );
      }

      TNL::MPI::Barrier( sph.communicator ); //To have clear output

      // perform interaction with given model
      sph.timeMeasurement.start( "interact" );
      sph.interact();
      sph.timeMeasurement.stop( "interact" );
      sph.writeLog( log, "Interact...", "Done." );

      // SingleSet: forgot overlaps
      sph.resetOverlaps();

      //integrate
      sph.timeMeasurement.start( "integrate" );
      sph.integrator->integratStepVerlet( sph.fluid, sph.boundary, sph.timeStepping, SPHDefs::BCType::integrateInTime() );
      sph.timeMeasurement.stop( "integrate" );
      sph.writeLog( log, "Integrate...", "Done." );

      //// output particle data
      //if( sph.timeStepping.checkOutputTimer( "save_results" ) )
      //{
      //   /**
      //    * Compute pressure from density.
      //    * This is not necessary since we do this localy, if pressure is needed.
      //    * It's useful for output anyway.
      //    */
      //   sph.model.computePressureFromDensity( sph.fluid, sph.modelParams );
      //   sph.model.computePressureFromDensity( sph.boundary, sph.modelParams );

      //   sph.save( log );
      //}

      // check timers and if measurement or interpolation should be performed, is performed
      sph.template measure< SPHDefs::KernelFunction, SPHDefs::EOS >( log );

      // check timers and if snapshot should be done, is done
      //sph.save( log );

      // update time step
      sph.timeStepping.updateTimeStep();

      TNL::MPI::Barrier( sph.communicator ); //To have clear output
   }

   sph.writeEpilog( log );
}

