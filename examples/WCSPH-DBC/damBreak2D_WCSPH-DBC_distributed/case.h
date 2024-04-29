#include <TNL/Config/parseCommandLine.h>
#include <TNL/Config/parseINIConfigFile.h>
#include <TNL/Logger.h>

#include <SPH/configSetup.h>
#include "template/config.h"

#include <SPH/Models/WCSPH_DBC/control.h>
#include <string>

int main( int argc, char* argv[] )
{

   TNL::MPI::ScopedInitializer mpi( argc, argv );

//test;
#ifdef HAVE_MPI
   std::cout << "Running with MPI." << std::endl;
#endif

   // get CLI parameters
   TNL::Config::ParameterContainer cliParams;
   TNL::Config::ConfigDescription cliConfig;

   cliConfig.addRequiredEntry< std::string >( "config", "Path to the configuration file." );
   cliConfig.addEntry< bool >( "config-help", "Print the configuration file description and exit.", false );
   cliConfig.addEntry< bool >( "print-static-configuration", "Print the static configuration (e.g. solver and model types) and exit.", false );
   cliConfig.addEntry< std::string >( "output-directory", "Path to the output directory (overrides the corresponding option in the configuration file)." );
   cliConfig.addEntry< int >( "verbose", "Set the verbose mode. The higher number the more messages are generated.", 2 );
   cliConfig.addEntry< std::string >( "log-file", "Log file for the computation.", "log.txt" );
   cliConfig.addEntry< int >( "log-width", "Number of columns of the log table.", 80 );
   cliConfig.addEntry< bool >( "catch-exceptions",
                               "Catch C++ exceptions. Disabling it allows the program to drop into the debugger "
                               "and track the origin of the exception.",
                               true );

   if( ! TNL::Config::parseCommandLine( argc, argv, cliConfig, cliParams ) )
       return EXIT_FAILURE;

   // get config parameters
   TNL::Config::ParameterContainer parameters;
   TNL::Config::ConfigDescription config;

   // set SPH parameters
   TNL::SPH::configSetup( config );

   // set model parameters
   //Model::configSetup( config );
   TNL::SPH::template configSetupModel< SPHConfig< Device > >( config );
   //TNL::SPH::configSetupModel( config );

   if( cliParams.getParameter< bool >( "config-help" ) ) {
       // TODO: re-format the message for the config (drop the program name and "--")
       std::cout << "Priting usage." << std::endl;
       TNL::Config::printUsage( config, argv[0] );
       return EXIT_SUCCESS;
   }
   if( cliParams.getParameter< bool >( "print-static-configuration" ) ) {
       const int logWidth = cliParams.getParameter< int >( "log-width" );
       TNL::Logger consoleLogger( logWidth, std::cout );
       //TNL::MHFEM::writeProlog< Problem >( consoleLogger, false );
       TNL::SPH::writeProlog< Simulation>( consoleLogger, false );
       return EXIT_SUCCESS;
   }

   const std::string configPath = cliParams.getParameter< std::string >( "config" );
   try {
       parameters = TNL::Config::parseINIConfigFile( configPath, config );
   }
   catch ( const std::exception& e ) {
       std::cerr << "Failed to parse the configuration file " << configPath << " due to the following error:\n" << e.what() << std::endl;
       return EXIT_FAILURE;
   }
   catch (...) {
       std::cerr << "Failed to parse the configuration file " << configPath << " due to an unknown C++ exception." << std::endl;
       throw;
   }

   // --output-directory from the CLI overrides output-directory from the config
   if( cliParams.checkParameter( "output-directory" ) )
       parameters.setParameter< std::string >( "output-directory", cliParams.getParameter< std::string >( "output-directory" ) );
   if( ! parameters.checkParameter("output-directory")) {
       std::cerr << "The output-directory parameter was not found in the config and "
                    "--output-directory was not given on the command line." << std::endl;
       return EXIT_FAILURE;
   }

   //TNL::Logger log( 100, std::cout );
   //Simulation sph;
   //TNL::MPI::Barrier( sph.communicator ); //To have clear output
   //sph.init( parameters, log );
   //TNL::MPI::Barrier( sph.communicator ); //To have clear output
   //if( TNL::MPI::GetRank() == 0 )
   //   sph.writeProlog( log );

   //DEBUG:
   std::string logFileName = "results/simulationLog_rank" + std::to_string( TNL::MPI::GetRank() );
   std::ofstream logFile( logFileName );
   //TNL::Logger log( 100, std::cout );
   TNL::Logger log( 100, logFile );

   Simulation sph;
   TNL::MPI::Barrier( sph.communicator ); //To have clear output

   if( TNL::MPI::GetRank() == 0 )
      sph.init( parameters, log );
   TNL::MPI::Barrier( sph.communicator ); //To have clear output
   if( TNL::MPI::GetRank() == 1 )
      sph.init( parameters, log );
   TNL::MPI::Barrier( sph.communicator ); //To have clear output
   if( TNL::MPI::GetRank() == 2 )
      sph.init( parameters, log );
   TNL::MPI::Barrier( sph.communicator ); //To have clear output

   if( TNL::MPI::GetRank() == 0 )
      sph.writeProlog( log );

   // Solver model:

   //sph.init( parameters );
   //sph.writeProlog( parameters );
   //sph.exec();
   //sph.writeEpilog( parameters );

   // Library model:
   //return 0;

   while( sph.timeStepping.runTheSimulation() )
   //while( sph.timeStepping.getStep() < 2 )
   {
      //sph.writeInfo( log );

      // SingleSet: forgot overlaps
      //sph.resetOverlaps();
      //sph.writeLog( log, "Reset overlaps...", "Done." );

      TNL::MPI::Barrier( sph.communicator ); //To have clear output

      // search for neighbros
      sph.timeMeasurement.start( "search" );
      sph.performNeighborSearch( log );
      sph.timeMeasurement.stop( "search" );
      sph.writeLog( log, "Search...", "Done." );

      //// output particle data
      ////if( sph.timeStepping.getStep() < 950 || sph.timeStepping.getStep() > 1050 ){
      //if( sph.timeStepping.getStep() < 1500 || sph.timeStepping.getStep() > 1700 ){
      //   if( sph.timeStepping.checkOutputTimer( "save_results" ) )
      //   {
      //      /**
      //       * Compute pressure from density.
      //       * This is not necessary since we do this localy, if pressure is needed.
      //       * It's useful for output anyway.
      //       */
      //      sph.model.computePressureFromDensity( sph.fluid, sph.modelParams );
      //      sph.model.computePressureFromDensity( sph.boundary, sph.modelParams );

      //      sph.save( log, true );
      //   }
      //}
      //else{
      //      std::cout << "=========================================================== STEP " << sph.timeStepping.getStep() << " =================== " << std::endl;
      //      sph.model.computePressureFromDensity( sph.fluid, sph.modelParams );
      //      sph.model.computePressureFromDensity( sph.boundary, sph.modelParams );

      //      sph.save( log, true );
      //}
      ////if( sph.timeStepping.getStep() == 1316 )
      ////   return 0;

      TNL::MPI::Barrier( sph.communicator ); //To have clear output


      // SingleSet + Overlaps: collect particles, send them between processors
      sph.synchronizeDistributedSimulation();
      sph.writeLog( log, "Synchronize...", "Done." );
      //return 0;

      TNL::MPI::Barrier( sph.communicator ); //To have clear output

      // SingleSet: perform second search to assign received particles
      sph.timeMeasurement.start( "search" );
      sph.performNeighborSearch( log );
      sph.timeMeasurement.stop( "search" );
      sph.writeLog( log, "Search second...", "Done." );

      TNL::MPI::Barrier( sph.communicator ); //To have clear output

      // output particle data
      //if( sph.timeStepping.getStep() < 950 || sph.timeStepping.getStep() > 1050 ){
      if( sph.timeStepping.getStep() < 200 || sph.timeStepping.getStep() > 500 ){
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
      }
      else{
            std::cout << "=========================================================== STEP " << sph.timeStepping.getStep() << " =================== " << std::endl;
            sph.model.computePressureFromDensity( sph.fluid, sph.modelParams );
            sph.model.computePressureFromDensity( sph.boundary, sph.modelParams );

            sph.save( log, true );
      }
      //if( sph.timeStepping.getStep() == 1316 )
      //   return 0;

      TNL::MPI::Barrier( sph.communicator ); //To have clear output


      //TNL::MPI::Barrier( sph.communicator ); //To have clear output
      //if( TNL::MPI::GetRank() == 0 )
      //   sph.synchronizeDistributedSimulation();
      //TNL::MPI::Barrier( sph.communicator ); //To have clear output
      //if( TNL::MPI::GetRank() == 1 )
      //   sph.synchronizeDistributedSimulation();
      //TNL::MPI::Barrier( sph.communicator ); //To have clear output
      //if( TNL::MPI::GetRank() == 2 )
      //   sph.synchronizeDistributedSimulation();
      //TNL::MPI::Barrier( sph.communicator ); //To have clear output

      // perform interaction with given model
      sph.timeMeasurement.start( "interact" );
      sph.interact();
      sph.timeMeasurement.stop( "interact" );
      sph.writeLog( log, "Interact...", "Done." );

      // SingleSet: forgot overlaps
      sph.resetOverlaps();
      //sph.writeLog( log, "Reset overlaps...", "Done." );

      //integrate
      sph.timeMeasurement.start( "integrate" );
      sph.integrator->integratStepVerlet( sph.fluid, sph.boundary, sph.timeStepping );
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
      std::cout << "============================================++++++++++++++++++++++++++++++++++++============================" << std::endl;
      TNL::MPI::Barrier( sph.communicator ); //To have clear output
   }

   sph.writeEpilog( log );
}

