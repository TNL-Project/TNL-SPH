#pragma once

#include <TNL/Config/parseCommandLine.h>
#include <TNL/Config/parseINIConfigFile.h>
#include <SPH/configSetup.h>
#include <cstdlib>

namespace TNL {
namespace SPH {

template< typename Simulation >
inline void
initialize( int argc,
            char* argv[],
            TNL::Config::ParameterContainer& cliParams,
            TNL::Config::ConfigDescription& cliConfig,
            TNL::Config::ParameterContainer& parameters,
            TNL::Config::ConfigDescription& config )
{

   cliConfig.addRequiredEntry< std::string >( "config", "Path to the configuration file." );
   cliConfig.addEntry< bool >( "config-help", "Print the configuration file description and exit.", false );
   cliConfig.addEntry< bool >( "print-static-configuration", "Print the static configuration (e.g. solver and model types) and exit.", false );
   cliConfig.addEntry< std::string >( "output-directory", "Path to the output directory (overrides the corresponding option in the configuration file)." );
   cliConfig.addEntry< int >( "verbose", "Set the verbose mode. The higher number the more messages are generated.", 2 );
   cliConfig.addEntry< std::string >( "log-file", "Log file for the computation.", "log.txt" );
   cliConfig.addEntry< int >( "log-width", "Number of columns of the log table.", 100 );
   cliConfig.addEntry< bool >( "catch-exceptions",
                               "Catch C++ exceptions. Disabling it allows the program to drop into the debugger "
                               "and track the origin of the exception.",
                               true );

   if( ! TNL::Config::parseCommandLine( argc, argv, cliConfig, cliParams ) )
       std::cerr << "Failed to parse the command line arguments." << std::endl;
       //return EXIT_FAILURE;

   // set SPH parameters
   TNL::SPH::configSetup( config );

   // set model parameters
   //Model::configSetup( config );
   //TNL::SPH::template configSetupModel< SPHConfig< Device > >( config );
   //TNL::SPH::configSetupModel( config );
   Simulation::ModelType::ModelParams::configSetupModel( config );

   if( cliParams.getParameter< bool >( "config-help" ) ) {
       // TODO: re-format the message for the config (drop the program name and "--")
       std::cout << "Priting usage." << std::endl;
       TNL::Config::printUsage( config, argv[0] );
       //return EXIT_SUCCESS;
       std::exit( EXIT_SUCCESS );
   }
   if( cliParams.getParameter< bool >( "print-static-configuration" ) ) {
       const int logWidth = cliParams.getParameter< int >( "log-width" );
       TNL::Logger consoleLogger( logWidth, std::cout );
       TNL::SPH::writeProlog< Simulation>( consoleLogger, false );
       //return EXIT_SUCCESS;
       std::exit( EXIT_SUCCESS );
   }

   const std::string configPath = cliParams.getParameter< std::string >( "config" );
   try {
       parameters = TNL::Config::parseINIConfigFile( configPath, config );
   }
   catch ( const std::exception& e ) {
       std::cerr << "Failed to parse the configuration file " << configPath << " due to the following error:\n" << e.what() << std::endl;
       //return EXIT_FAILURE;
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
       //return EXIT_FAILURE;
   }

}

} // SPH
} // TNL

