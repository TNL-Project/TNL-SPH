// SPDX-FileComment: This file is part of TNL-SPH (https://gitlab.com/tnl-project/tnl-sph)
// SPDX-License-Identifier: MIT

#pragma once
#include <TNL/Devices/Cuda.h>
#include <TNL/Devices/Host.h>
#include <TNL/MPI.h>
#include <TNL/Config/ConfigDescription.h>

#include <SPH/parseConfigFile.h>

namespace TNL::SPH {

void
configSetup( TNL::Config::ConfigDescription& config, std::string sectionPrefix = "SPH" )
{
   config.addEntry< std::string >( "caseName", "Name of solved case." );
   config.addEntry< std::string >( "outputDirectory", "Path to the output directory." );

   config.addDelimiter( sectionPrefix + " space discretisation" );
   config.addRequiredEntry< std::string >( "fluidParticles", "Input fluid particles file path." );
   config.addRequiredEntry< std::string >( "boundaryParticles", "Input boundary particles file path." );

   config.addEntry< std::string >( "particlesFormat", "Input pa file format.", "auto" );
   config.addEntryEnum( "auto" );
   config.addEntryEnum( "vtk" );
   config.addEntryEnum( "vtu" );
   config.addEntry< std::string >( "boundaryConditionsFile", "Path to the boundary conditions file." );
   config.addEntry< std::string >( "verboseIntensity", "Information write during simulation runtime.", "with-snapshot" );
   config.addEntryEnum( "full" );
   config.addEntryEnum( "with-snapshot" );
   config.addEntryEnum( "none" );
   config.addEntry< bool >( "saveFluidParticles", "Save fluid particles (bool).", true );
   config.addEntry< bool >( "saveBoundaryParticles", "Save fluid particles (bool).", true );

   config.addDelimiter( sectionPrefix + " time discretisation" );
   config.addEntry< std::string >( "initialCondition", "File name with the initial condition." );
   config.addRequiredEntry< double >( "finalTime", "Stop time of the time dependent problem." );
   config.addEntry< double >( "initialTimeStep", "Initial time of the time dependent problem.", 0 );
   config.addEntry< double >( "startTime", "Initial simulation time for restart (0 = classical mode).", 0. );
   config.addRequiredEntry< double >( "snapshotPeriod", "Time period for writing the problem status." );
   config.addEntry< double >( "timeStep", "The time step for the time discretisation.", 1.0 );

   config.addEntry< double >( "searchRadius", "The cut of radius of particle interactions.", 0. );
   config.addEntry< int >( "numberOfParticles", "The initial number of fluid particles.", 0. );
   config.addEntry< int >( "numberOfAllocatedParticles", "The allocated number of fluid particles.", 0. );
   config.addEntry< int >( "numberOfBoundaryParticles", "The initial number of fluid particles.", 0. );
   config.addEntry< int >( "numberOfAllocatedBoundaryParticles", "The allocated number of fluid particles.", 0. );

   config.addEntry< double >( "domainOriginX", "The origin of domain in x direction.", 0. );
   config.addEntry< double >( "domainOriginY", "The origin of domain in y direction.", 0. );
   config.addEntry< double >( "domainOriginZ", "The origin of domain in z direction.", 0. );
   config.addEntry< double >( "domainSizeX", "The size of domain in x direction.", 0. );
   config.addEntry< double >( "domainSizeY", "The size of domain in y direction.", 0. );
   config.addEntry< double >( "domainSizeZ", "The size of domain in y direction.", 0. );

   config.addEntry< int >( "openBoundaryPatches", "Number of open boundary patches.", 0 );
   config.addEntry< std::string >( "openBoundaryConfig", "Configuration file for open boundary.", "" );
   config.addEntry< std::string >(
      "openBoundary", "Inline open boundary config as a JSON object. Used when openBoundaryConfig is empty.", "" );
   config.addEntry< int >( "periodicBoundaryPatches", "Number of periodic boundary patces.", 0 );
   config.addEntry< std::string >( "periodicBoundaryConfig", "Configuration file for periodic boundary.", "" );
   config.addEntry< std::string >(
      "periodicBoundary", "Inline periodic boundary config as a JSON object. Used when periodicBoundaryConfig is empty.", "" );

   config.addEntry< int >( "numberOfSubdomains", "Number of particles subdomains.", 0 );
   config.addEntry< std::string >( "subdomainsConfig", "Configuration file for subdomains setup.", "" );
   config.addEntry< std::string >(
      "subdomains", "Inline subdomains config as a JSON object. Used when subdomainsConfig is empty.", "" );
   config.addEntry< std::string >(
      "boundaryGhostMethod", "Update method for boundary ghost particles at multiresolution interfaces.", "interpolation" );
   config.addEntryEnum( "interpolation" );
   config.addEntryEnum( "direct-from-source" );

   // distributed simulation parameters
   config.addEntry< int >( "subdomainsX", "Number of subdomains in the x direstion.", 0 );
   config.addEntry< int >( "subdomainsY", "Number of subdomains in the y direstion.", 0 );
   config.addEntry< std::string >( "distributedConfig", "Path to the config with distributed simulation data.", "" );
   config.addEntry< std::string >(
      "distributedDomain", "Inline distributed domain config as a JSON object. Used when distributed-config is empty.", "" );
   config.addEntry< int >( "overlapWidth", "Width in cells around every domain", 1 );
   config.addEntry< std::string >( "loadBalancingMeasure", "Measure for subdomains sizes rebalancing.", "computationalTime" );
   config.addEntryEnum( "numberOfParticles" );
   config.addEntryEnum( "computationalTime" );
   config.addEntry< float >( "numberOfParticlesBalancingCoef", "Particles count treshold for rebalancing [%].", 0.05 );
   config.addEntry< float >( "computationalTimeBalancingCoef", "Comp. time treshold for rebalancing.", 0.001 );
   config.addEntry< int >( "loadBalancingStepInterval", "Number of steps after load balancing is performed", 100 );

   // simulation monitor parameters
   config.addEntry< std::string >( "measuretoolConfig", "Configuration file for the measuretool config.", "" );
   config.addEntry< std::string >(
      "measuretool", "Inline measuretool config as a JSON object. Used when measuretoolConfig is empty.", "" );
   config.addEntry< int >( "interpolationPlanesCount", "Number of interpolation planes.", 0 );
   config.addEntry< int >( "pressureSensorsCount", "Number of pressure evalutaion sensors.", 0 );
   config.addEntry< int >( "waterLevelSensorsCount", "Number of water level evaluation sensors.", 0 );
   config.addEntry< int >( "volumetricFlowRatePlanesCount", "Number of volumetric flow rate evaluation planes", 0 );

   //TODO: Move this to suiteble place, it is used also for open zones
   config.addEntry< int >( "numberOfParticlesPerCell", "Max allowed number of particles per cell", 15 );

   // user defined parameters
   config.addEntry< std::string >( "userDefinedConfig", "Configuration file for user's custom parameters.", "" );
   config.addEntry< std::string >(
      "userDefined", "Inline user-defined config as a JSON object. Used when user-defined-config is empty.", "" );
}

void
configSetupDistributedSubdomain( int subdomain_x, int subdomain_y, TNL::Config::ConfigDescription& config )
{
   std::string subdomainKey = "subdomainX" + std::to_string( subdomain_x ) + "Y" + std::to_string( subdomain_y );
   config.addRequiredEntry< std::string >( subdomainKey + "FluidParticles", "Input fluid particles file path." );
   config.addRequiredEntry< std::string >( subdomainKey + "BoundaryParticles", "Input boundary particles file path." );
   config.addEntry< int >( subdomainKey + "Fluid_n", "The initial number of fluid particles.", 0 );
   config.addEntry< int >( subdomainKey + "Fluid_n_allocated", "The allocated number of fluid particles.", 0 );
   config.addEntry< int >( subdomainKey + "Boundary_n", "The initial number of fluid particles.", 0 );
   config.addEntry< int >( subdomainKey + "Boundary_n_allocated", "The allocated number of fluid particles.", 0 );

   config.addEntry< double >( subdomainKey + "RefinementFactor", "Refinement factor of current subdomain.", 0. );

   config.addEntry< double >( subdomainKey + "OriginX", "The origin of domain in x direction.", 0. );
   config.addEntry< double >( subdomainKey + "OriginY", "The origin of domain in y direction.", 0. );
   config.addEntry< double >( subdomainKey + "OriginZ", "The origin of domain in z direction.", 0. );
   config.addEntry< int >(
      subdomainKey + "OriginGlobalCoordsX", "The origin of domain in global cell coords. in x direction.", 0. );
   config.addEntry< int >(
      subdomainKey + "OriginGlobalCoordsY", "The origin of domain in global cell coords. in y direction.", 0. );
   config.addEntry< int >(
      subdomainKey + "OriginGlobalCoordsZ", "The origin of domain in global cell coords. in z direction.", 0. );

   config.addEntry< double >( subdomainKey + "SizeX", "The size of domain in x direction.", 0. );
   config.addEntry< double >( subdomainKey + "SizeY", "The size of domain in y direction.", 0. );
   config.addEntry< double >( subdomainKey + "SizeZ", "The size of domain in y direction.", 0. );
   config.addEntry< int >( subdomainKey + "GridDimensionsX", "The size of domain in cells in x direction.", 0. );
   config.addEntry< int >( subdomainKey + "GridDimensionsY", "The size of domain in cells in y direction.", 0. );
   config.addEntry< int >( subdomainKey + "GridDimensionsZ", "The size of domain in cells in z direction.", 0. );
}

/**
 * Config entries for multiresolution buffer particles.
 *
 * Each multiresolution interface has a buffer patch on each side.
 * The buffer particles are pre-generated by the init script and read
 * from VTK files at solver startup.
 *
 * \param bufferIdx 0-based buffer index (0, 1, ...)
 * \param config    ConfigDescription to add entries to
 */
void
configMultiresolutionBuffer( int bufferIdx, TNL::Config::ConfigDescription& config )
{
   std::string key = "multiresolutionBuffer" + std::to_string( bufferIdx );
   config.addEntry< std::string >( key + "Particles", "Input multiresolution buffer particles file path.", "" );
   config.addEntry< int >( key + "N", "The initial number of multiresolution buffer particles.", 0 );
}

/**
 * Config entries for boundary ghost particles of one multiresolution interface.
 *
 * The ghost geometry (own-resolution wall lattice clipped to the interface
 * band) is pre-generated by the init script and read from a VTK file at
 * solver startup. Interfaces without a ghost file get no ghost particles.
 *
 * \param bufferIdx 0-based interface index (matches multiresolution-buffer indexing)
 * \param config    ConfigDescription to add entries to
 */
void
configBoundaryGhostBuffer( int bufferIdx, TNL::Config::ConfigDescription& config )
{
   std::string key = "boundaryGhostBuffer" + std::to_string( bufferIdx );
   config.addEntry< std::string >( key + "Particles", "Input boundary ghost particles file path.", "" );
   config.addEntry< int >( key + "N", "The initial number of boundary ghost particles.", 0 );
}

template< typename Simulation >
void
writeProlog( TNL::Logger& logger, bool writeSystemInformation = true )
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
parseOpenBoundaryConfig(
   const std::string& configOpenBoundaryPath,
   TNL::Config::ParameterContainer& parametersOpenBoundary,
   TNL::Config::ConfigDescription& configOpenBoundary,
   TNL::Logger& logger,
   const std::string& inlineConfig = "" )
{
   if( configOpenBoundaryPath != "" ) {
      logger.writeParameter( "Parsing open boundaries simulation config.", "" );
      try {
         parametersOpenBoundary = TNL::SPH::parseConfigFile( configOpenBoundaryPath, configOpenBoundary );
      }
      catch( const std::exception& e ) {
         std::cerr << "Failed to parse the open boundary configuration file " << configOpenBoundaryPath
                   << " due to the following error:\n"
                   << e.what() << std::endl;
      }
      catch( ... ) {
         std::cerr << "Failed to parse the open boundary configuration file " << configOpenBoundaryPath
                   << " due to an unknown C++ exception." << std::endl;
         throw;
      }
      logger.writeParameter( "Parsing open boundaries config.", "Done." );
   }
   else if( ! inlineConfig.empty() ) {
      logger.writeParameter( "Parsing open boundaries config.", "inline (config.jsonc subsection)" );
      try {
         parametersOpenBoundary = TNL::SPH::parseJSONConfigString( inlineConfig, "inline open-boundary", configOpenBoundary );
      }
      catch( const std::exception& e ) {
         std::cerr << "Failed to parse the inline open boundary configuration due to the following error:\n"
                   << e.what() << std::endl;
      }
      logger.writeParameter( "Parsing open boundaries config.", "Done." );
   }
}

void
parseDistributedConfig(
   const std::string& configPath,
   TNL::Config::ParameterContainer& params,
   TNL::Config::ConfigDescription& config,
   TNL::Logger& logger,
   const std::string& inlineConfig = "" )
{
   if( configPath != "" ) {
      logger.writeParameter( "Parsing distributed simulation config.", "" );
      try {
         params = TNL::SPH::parseConfigFile( configPath, config );
      }
      catch( const std::exception& e ) {
         std::cerr << "Failed to parse the distributed configuration file " << configPath << " due to the following error:\n"
                   << e.what() << std::endl;
      }
      catch( ... ) {
         std::cerr << "Failed to parse the distributed configuration file " << configPath << " due to an unknown C++ exception."
                   << std::endl;
         throw;
      }
      logger.writeParameter( "Parsing distributed simulation config.", "Done." );
   }
   else if( ! inlineConfig.empty() ) {
      logger.writeParameter( "Parsing distributed simulation config.", "inline (config.jsonc subsection)" );
      try {
         params = TNL::SPH::parseJSONConfigString( inlineConfig, "inline distributed-domain", config );
      }
      catch( const std::exception& e ) {
         std::cerr << "Failed to parse the inline distributed configuration due to the following error:\n"
                   << e.what() << std::endl;
      }
      logger.writeParameter( "Parsing distributed simulation config.", "Done." );
   }
}

void
parseUserDefinedConfig(
   const std::string& configPath,
   TNL::Config::ParameterContainer& params,
   TNL::Config::ConfigDescription& config,
   TNL::Logger& logger,
   const std::string& inlineConfig = "" )
{
   if( configPath != "" ) {
      logger.writeParameter( "Parsing user defined params config.", "" );
      try {
         params = TNL::SPH::parseConfigFile( configPath, config );
      }
      catch( const std::exception& e ) {
         std::cerr << "Failed to parse the user defined params config file " << configPath << " due to the following error:\n"
                   << e.what() << std::endl;
      }
      catch( ... ) {
         std::cerr << "Failed to parse the user defined params config file " << configPath
                   << " due to an unknown C++ exception." << std::endl;
         throw;
      }
      logger.writeParameter( "Parsing user defined params config.", "Done." );
   }
   else if( ! inlineConfig.empty() ) {
      logger.writeParameter( "Parsing user defined params config.", "inline (config.jsonc subsection)" );
      try {
         params = TNL::SPH::parseJSONConfigString( inlineConfig, "inline user-defined", config );
      }
      catch( const std::exception& e ) {
         std::cerr << "Failed to parse the inline user defined params config due to the following error:\n"
                   << e.what() << std::endl;
      }
      logger.writeParameter( "Parsing user defined params config.", "Done." );
   }
}

}  // namespace TNL::SPH
