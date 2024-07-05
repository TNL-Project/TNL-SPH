#pragma once
#include <TNL/Devices/Cuda.h>
#include <TNL/Devices/Host.h>
#include <TNL/MPI.h>
#include <TNL/Config/ConfigDescription.h>

namespace TNL {
namespace SPH {

inline void
configSetup( TNL::Config::ConfigDescription& config,
             std::string sectionPrefix = "SPH" )
{
    config.addEntry< std::string >( "case-name", "Name of solved case." );
    config.addEntry< std::string >( "output-directory", "Path to the output directory." );

    config.addDelimiter( sectionPrefix + " space discretisation" );
    config.addRequiredEntry< std::string >( "fluid-particles", "Input fluid particles file path." );
    config.addRequiredEntry< std::string >( "boundary-particles", "Input boundary particles file path." );

    config.addEntry< std::string >( "particles-format", "Input pa file format.", "auto" );
        config.addEntryEnum( "auto" );
        config.addEntryEnum( "vtk" );
        config.addEntryEnum( "vtu" );
    config.addEntry< std::string >( "boundary-conditions-file", "Path to the boundary conditions file." );
    config.addEntry< std::string >( "verbose-intensity", "Information write during simulation runtime." ,"with-snapshot" );
        config.addEntryEnum( "full" );
        config.addEntryEnum( "with-snapshot" );
        config.addEntryEnum( "none" );

    config.addDelimiter( sectionPrefix + " time discretisation" );
    config.addEntry< std::string >( "initial-condition", "File name with the initial condition." );
    config.addRequiredEntry< double >( "final-time", "Stop time of the time dependent problem." );
    config.addEntry< double >( "initial-time-step", "Initial time of the time dependent problem.", 0 );
    config.addRequiredEntry< double >( "snapshot-period", "Time period for writing the problem status.");
    config.addEntry< double >( "time-step", "The time step for the time discretisation.", 1.0 );

    config.addEntry< double >( "searchRadius", "The cut of radius of particle interactions.", 0. );
    config.addEntry< int >( "numberOfParticles", "The initial number of fluid particles.", 0. );
    config.addEntry< int >( "numberOfAllocatedParticles", "The allocated number of fluid particles.", 0. );
    config.addEntry< int >( "numberOfBoundaryParticles", "The initial number of fluid particles.", 0. );
    config.addEntry< int >( "numberOfAllocatedBoundaryParticles", "The allocated number of fluid particles.", 0. );

    config.addEntry< double >( "domainOrigin-x", "The origin of domain in x direction.", 0. );
    config.addEntry< double >( "domainOrigin-y", "The origin of domain in y direction.", 0. );
    config.addEntry< double >( "domainOrigin-z", "The origin of domain in z direction.", 0. );
    config.addEntry< double >( "domainSize-x", "The size of domain in x direction.", 0. );
    config.addEntry< double >( "domainSize-y", "The size of domain in y direction.", 0. );
    config.addEntry< double >( "domainSize-z", "The size of domain in y direction.", 0. );

    config.addEntry< int >( "openBoundaryPatches", "Number of open boundary patches.", 0 );
    config.addEntry< int >( "periodicBoundaryPatches", "Number of periodic boundary patces.", 0 );

    // distributed simulation parameters
    config.addEntry< int >( "subdomains-x", "Number of subdomains in the x direstion.", 0 );
    config.addEntry< int >( "subdomains-y", "Number of subdomains in the y direstion.", 0 );
    config.addEntry< std::string >( "distributed-config", "Path to the config with distributed simulation data.", "" );
    config.addEntry< int >( "overlapWidth", "Width in cells around every domain", 1 );

    // simulation monitor parameters
    config.addEntry< std::string >( "measuretool-config", "Configuration file for the measuretool config.", "" );
    config.addEntry< int >( "interpolation-planes-count", "Input boundary particles file path.", 0 );
    config.addEntry< int >( "pressure-sensors-count", "Input boundary particles file path.", 0 );
    config.addEntry< int >( "water-level-sensors-count", "Input boundary particles file path.", 0 );

    //TODO: Move this to suiteble place, it is used also for open zones
    config.addEntry< int >( "numberOfParticlesPerCell", "Max allowed number of particles per cell", 15 );
}

void
configSetupDistributedSubdomain( int subdomain_x, int subdomain_y, TNL::Config::ConfigDescription& config )
{
   std::string subdomainKey = "subdomain-x-" + std::to_string( subdomain_x ) + "-y-" + std::to_string( subdomain_y ) + "-";
   config.addRequiredEntry< std::string >( subdomainKey + "fluid-particles", "Input fluid particles file path." );
   config.addRequiredEntry< std::string >( subdomainKey + "boundary-particles", "Input boundary particles file path." );
   config.addEntry< int >( subdomainKey + "fluid_n", "The initial number of fluid particles.", 0 );
   config.addEntry< int >( subdomainKey + "fluid_n_allocated", "The allocated number of fluid particles.", 0 );
   config.addEntry< int >( subdomainKey + "boundary_n", "The initial number of fluid particles.", 0 );
   config.addEntry< int >( subdomainKey + "boundary_n_allocated", "The allocated number of fluid particles.", 0 );

   config.addEntry< double >( subdomainKey + "origin-x", "The origin of domain in x direction.", 0. );
   config.addEntry< double >( subdomainKey + "origin-y", "The origin of domain in y direction.", 0. );
   config.addEntry< double >( subdomainKey + "origin-z", "The origin of domain in z direction.", 0. );
   config.addEntry< int >( subdomainKey + "origin-global-coords-x", "The origin of domain in global cell coords. in x direction.", 0. );
   config.addEntry< int >( subdomainKey + "origin-global-coords-y", "The origin of domain in global cell coords. in y direction.", 0. );
   config.addEntry< int >( subdomainKey + "origin-global-coords-z", "The origin of domain in global cell coords. in z direction.", 0. );

   config.addEntry< double >( subdomainKey + "size-x", "The size of domain in x direction.", 0. );
   config.addEntry< double >( subdomainKey + "size-y", "The size of domain in y direction.", 0. );
   config.addEntry< double >( subdomainKey + "size-z", "The size of domain in y direction.", 0. );
   config.addEntry< int >( subdomainKey + "grid-dimensions-x", "The size of domain in cells in x direction.", 0. );
   config.addEntry< int >( subdomainKey + "grid-dimensions-y", "The size of domain in cells in y direction.", 0. );
   config.addEntry< int >( subdomainKey + "grid-dimensions-z", "The size of domain in cells in z direction.", 0. );
}


template< typename Simulation >
void writeProlog( TNL::Logger& logger, bool writeSystemInformation = true )
{
    const bool printGPUs = std::is_same< typename Simulation::DeviceType, TNL::Devices::Cuda >::value;

    //logger.writeHeader( Problem::getPrologHeader() );
    if( TNL::MPI::isInitialized() )
        logger.writeParameter( "MPI processes:", TNL::MPI::GetSize() );
    logger.writeParameter( "Device type:", TNL::getType< typename Simulation::DeviceType >() );

    if( ! printGPUs ) {
        if( TNL::Devices::Host::isOMPEnabled() ) {
            logger.writeParameter( "OMP enabled:", "yes", 1 );
            logger.writeParameter( "OMP threads:", TNL::Devices::Host::getMaxThreadsCount(), 1 );
        }
        else
            logger.writeParameter( "OMP enabled:", "no", 1 );
    }
}

void
parseDistributedConfig( const std::string& configDistributedPath,
                        TNL::Config::ParameterContainer& parametersDistributed,
                        TNL::Config::ConfigDescription& configDistributed,
                        TNL::Logger& logger )
{
   if( configDistributedPath != "" ) {
      logger.writeParameter( "Parsing distributed simulation config.", "" );
      try {
          parametersDistributed = TNL::Config::parseINIConfigFile( configDistributedPath, configDistributed );
      }
      catch ( const std::exception& e ) {
          std::cerr << "Failed to parse the measuretool configuration file " << configDistributedPath << " due to the following error:\n" << e.what() << std::endl;
      }
      catch (...) {
          std::cerr << "Failed to parse the measuretool configuration file " << configDistributedPath << " due to an unknown C++ exception." << std::endl;
          throw;
      }
      logger.writeParameter( "Parsing distributed simulation config.", "Done." );
   }
}


} // SPH
} // TNL

