#include "SPHMultised_CFDDomain.h"

namespace TNL {
namespace SPH {

template< typename Model >
void
SPHMultiset_CFDDomain< Model >::init( const std::string& subdomainKey,
                                      const TNL::Config::ParameterContainer& parameters,
                                      TNL::Logger& logger )
{
   logger.writeHeader( "SPH domain initialization" );
   // initialize particle sets in the domain
   initializeParticleSets( subdomainKey, parameters, logger );
#ifdef HAVE_MPI
   //initialize distributed staff
#endif
   //initialize overlaps (if overlaps)
   readParticleSets( subdomainKey, domainParameters, logger );
   //initialize model parameters
   modelParams.init( parameters );
   // initialize time stepping
   timeStepping.setTimeStep( parameters.getParameter< RealType >( "initial-time-step" ) );
   timeStepping.setEndTime( parameters.getParameter< RealType >( "final-time" ) );
   timeStepping.addOutputTimer( "save_results", parameters.getParameter< RealType >( "snapshot-period" ) );
   // initialize the measuretool
   logger.writeSeparator();
   if( parameters.getParameter< std::string >( "measuretool-config" ) != "" ) {
      logger.writeParameter( "Simulation monitor initialization.", "" );
      simulationMonitor.init( parameters, timeStepping, logger );
      logger.writeParameter( "Simulation monitor initialization.", "Done." );
   }
   logger.writeHeader( "SPH domain " + subdomainKey + " successfully initialized." );
}

template< typename Model >
void
SPHMultiset_CFDDomain< Model >::initParticleSets( const std::string& subdomainKey,
                                                  const TNL::Config::ParameterContainer& domainParameters,
                                                  TNL::Logger& logger )
{
   // compute domain properetis
   const VectorType domainOrigin = parameters.getXyz< VectorType >( subdomainKey + "origin" );
   const VectorType domainSize = parameters.getXyz< VectorType >( subdomainKey + "size" );
   const RealType searchRadius = parameters.getParameter< RealType >( subdomainKey + "searchRadius" );
   const IndexVectorType gridDimensions = TNL::ceil( ( domainSize - domainOrigin ) / searchRadius );
   // initialize main particle sets
   fluid->initialize( parameters.getParameter< int >( subdomainKey + "numberOfParticles" ),
                      parameters.getParameter< int >( subdomainKey + "numberOfAllocatedParticles" ),
                      searchRadius,
                      gridDimensions,
                      domainOrigin );
   boundary->initialize( parameters.getParameter< int >( subdomainKey + "numberOfBoundaryParticles" ),
                         parameters.getParameter< int >( subdomainKey + "numberOfAllocatedBoundaryParticles" ),
                         searchRadius,
                         gridDimensions,
                         domainOrigin );
   // init open boundary patches
   const int numberOfBoundaryPatches = parameters.getParameter< int >( "openBoundaryPatches" );
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ) {  //TODO: I dont like this.
      openBoundaryPatches.resize( numberOfBoundaryPatches );
      for( int i = 0; i < numberOfBoundaryPatches; i++ ) {
         std::string prefix = subdomainKey + "buffer-" + std::to_string( i + 1 ) + "-";
         openBoundaryPatches[ i ]->config.init( parameters, prefix );
         openBoundaryPatches[ i ]->initialize( parameters.getParameter< int >( prefix + "numberOfParticles" ),
                                               parameters.getParameter< int >( prefix + "numberOfAllocatedParticles" ),
                                               searchRadius,
                                               gridDimensions,
                                               domainOrigin );
      }
   }
}

template< typename Model >
void
SPHMultiset_CFDDomain< Model >::readParticleSets( const std::string& subdomainKey,
                                                  const TNL::Config::ParameterContainer& domainParameters,
                                                  TNL::Logger& logger )
{
   // read particle data
   logger.writeParameter( "Reading fluid particles:", domainParameters.getParameter< std::string >( subdomainKey + "fluid-particles" ) );
   fluid->template readParticlesAndVariables< SimulationReaderType >(
      domainParameters.getParameter< std::string >( subdomainKey + "fluid-particles" ) );
   logger.writeParameter( "Reading boundary particles:", domainParameters.getParameter< std::string >( subdomainKey + "boundary-particles" ) );
   boundary->template readParticlesAndVariables< SimulationReaderType >(
      domainParameters.getParameter< std::string >( subdomainKey + "boundary-particles" ) );

   const int numberOfBoundaryPatches = domainParameters.getParameter< int >( "openBoundaryPatches" );
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ) {  //TODO: I dont like this.
      for( int i = 0; i < numberOfBoundaryPatches; i++ ) {
         std::string prefix = subdomainKey + "buffer-" + std::to_string( i + 1 ) + "-";
         openBoundaryPatches[ i ]->template readParticlesAndVariables< SimulationReaderType >(
            domainParameters.getParameter< std::string >( prefix + "particles" ) );
      }
   }
}

}  //namespace SPH
}  //namespace TNL
