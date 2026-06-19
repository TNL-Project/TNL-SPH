#include "SolverMultiSetRemeshed.h"

namespace TNL {
namespace SPH {

template< typename Model >
void
SolverMultiSetRemeshed< Model >::initRemeshedSimulation( int argc, char* argv[] )
{
   auto& params = this->parameters;
   auto& log = this->logger;

   try {
      initialize< SimulationType >( argc, argv, this->cliParams, this->cliConfig, params, this->config );
   }
   catch ( ... ) {
      std::cerr << std::endl;
   }

   log.writeHeader( "Remeshed SPH simulation initialization." );

   const VectorType domainOrigin = params.template getXyz< VectorType >( "domainOrigin" );
   const VectorType domainSize = params.template getXyz< VectorType >( "domainSize" );
   const RealType searchRadius = params.template getParameter< RealType >( "searchRadius" );
   const IndexVectorType gridSize = TNL::ceil( ( domainSize - domainOrigin ) / searchRadius );

   const int numberOfParticles = params.template getParameter< int >( "numberOfParticles" );
   const int numberOfAllocatedParticles = params.template getParameter< int >( "numberOfAllocatedParticles" );
   const int numberOfBoundaryParticles = params.template getParameter< int >( "numberOfBoundaryParticles" );
   const int numberOfAllocatedBoundaryParticles = params.template getParameter< int >( "numberOfAllocatedBoundaryParticles" );

   this->numberOfSubsets = 2;
   this->fluidSets.resize( this->numberOfSubsets );
   this->boundarySets.resize( this->numberOfSubsets );

   this->fluidSets[ 0 ]->initialize( numberOfParticles,
                                      numberOfAllocatedParticles,
                                      searchRadius,
                                      gridSize,
                                      domainOrigin );
   if constexpr( ParticlesType::specifySearchedSetExplicitly() == true )
      this->fluidSets[ 0 ]->getParticles()->setParticleSetLabel( 0 );

   this->fluidSets[ 1 ]->initialize( 0,
                                      numberOfAllocatedParticles,
                                      searchRadius,
                                      gridSize,
                                      domainOrigin );
   if constexpr( ParticlesType::specifySearchedSetExplicitly() == true )
      this->fluidSets[ 1 ]->getParticles()->setParticleSetLabel( 0 );

   this->boundarySets[ 0 ]->initialize( numberOfBoundaryParticles,
                                         numberOfAllocatedBoundaryParticles,
                                         searchRadius,
                                         gridSize,
                                         domainOrigin );
   if constexpr( ParticlesType::specifySearchedSetExplicitly() == true )
      this->boundarySets[ 0 ]->getParticles()->setParticleSetLabel( 1 );

   this->boundarySets[ 1 ]->initialize( numberOfBoundaryParticles,
                                         numberOfAllocatedBoundaryParticles,
                                         searchRadius,
                                         gridSize,
                                         domainOrigin );
   if constexpr( ParticlesType::specifySearchedSetExplicitly() == true )
      this->boundarySets[ 1 ]->getParticles()->setParticleSetLabel( 1 );

   const std::string fluidFileName = params.template getParameter< std::string >( "fluid-particles" );
   const std::string boundaryFileName = params.template getParameter< std::string >( "boundary-particles" );
   log.writeParameter( "Reading fluid particles:", fluidFileName );
   this->fluidSets[ 0 ]->template readParticlesAndVariables< typename BaseType::SimulationReaderType >( fluidFileName );
   log.writeParameter( "Reading boundary particles:", boundaryFileName );
   this->boundarySets[ 0 ]->template readParticlesAndVariables< typename BaseType::SimulationReaderType >( boundaryFileName );
   this->boundarySets[ 1 ]->template readParticlesAndVariables< typename BaseType::SimulationReaderType >( boundaryFileName );

   this->modelParams.init( params );

   this->timeStepping.setTimeStep( params.template getParameter< RealType >( "initial-time-step" ) );
   this->timeStepping.setEndTime( params.template getParameter< RealType >( "final-time" ) );
   this->timeStepping.addOutputTimer( "save_results", params.template getParameter< RealType >( "snapshot-period" ) );

   this->caseName = params.template getParameter< std::string >( "case-name" );
   this->verbose = params.template getParameter< std::string >( "verbose-intensity" );
   this->outputDirectory = params.template getParameter< std::string >( "output-directory" );
   this->particlesFormat = params.template getParameter< std::string >( "particles-format" );

   this->timeMeasurement.addTimer( "remesh" );

   log.writeSeparator();
   if( params.template getParameter< std::string >( "measuretool-config" ) != "" ) {
      log.writeParameter( "Simulation monitor initialization.", "" );
      this->simulationMonitor.init( params, this->timeStepping, log );
      log.writeParameter( "Simulation monitor initialization.", "Done." );
   }

   log.writeHeader( "Remeshed SPH simulation successfully initialized." );
}

template< typename Model >
void
SolverMultiSetRemeshed< Model >::interact()
{
   this->timeMeasurement.start( "interact" );

   this->model.updateSolidBoundary( this->fluidSets[ 0 ], this->boundarySets[ 0 ], this->modelParams );
   this->model.finalizeBoundaryInteraction( this->fluidSets[ 0 ], this->boundarySets[ 0 ], this->modelParams );

   this->model.interaction( this->fluidSets[ 0 ], this->boundarySets[ 0 ], this->modelParams );
   this->model.finalizeInteraction( this->fluidSets[ 0 ], this->boundarySets[ 0 ], this->modelParams );

   this->timeMeasurement.stop( "interact" );
   this->writeLog( "Interact...", "Done." );
}

template< typename Model >
void
SolverMultiSetRemeshed< Model >::remeshParticles()
{
   this->timeMeasurement.start( "remesh" );
   remeshing::remeshParticles< ParticlesType, SPHConfig >( this->fluidSets[ 0 ], this->fluidSets[ 1 ], this->modelParams );
   this->timeMeasurement.stop( "remesh" );
   this->writeLog( "Remesh...", "Done." );
}

} // namespace SPH
} // namespace TNL
