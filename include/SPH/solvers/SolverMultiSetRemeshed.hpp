// SPDX-FileComment: This file is part of TNL-SPH (https://gitlab.com/tnl-project/tnl-sph)
// SPDX-License-Identifier: MIT

#pragma once

#include <SPH/solvers/SolverMultiSetRemeshed.h>
#include <SPH/ParseConfigFile.h>

namespace TNL::SPH {

template< typename Model >
void
SolverMultiSetRemeshed< Model >::initRemeshedSimulation( int argc, char* argv[] )
{
   auto& params = this->parameters;
   auto& log = this->logger;

   try {
      initialize< SimulationType >( argc, argv, this->cliParams, this->cliConfig, params, this->config );
   }
   catch( ... ) {
      std::cerr << std::endl;
   }

   log.writeHeader( "Remeshed SPH simulation initialization." );

   const VectorType domainOrigin = getXyzCamelCase< VectorType >( params, "domainOrigin" );
   const VectorType domainSize = getXyzCamelCase< VectorType >( params, "domainSize" );
   const RealType searchRadius = params.template getParameter< RealType >( "searchRadius" );
   const CoordinatesType dimensions = TNL::ceil( ( domainSize - domainOrigin ) / searchRadius );

   const int numberOfParticles = params.template getParameter< int >( "numberOfParticles" );
   const int numberOfAllocatedParticles = params.template getParameter< int >( "numberOfAllocatedParticles" );
   const int numberOfBoundaryParticles = params.template getParameter< int >( "numberOfBoundaryParticles" );
   const int numberOfAllocatedBoundaryParticles = params.template getParameter< int >( "numberOfAllocatedBoundaryParticles" );

   this->numberOfSubsets = 2;
   this->fluidSets.resize( this->numberOfSubsets );
   this->boundarySets.resize( this->numberOfSubsets );

   this->fluidSets[ 0 ]->initialize( numberOfParticles, numberOfAllocatedParticles, searchRadius, dimensions, domainOrigin );
   if constexpr( ParticlesType::specifySearchedSetExplicitly() == true )
      this->fluidSets[ 0 ]->getParticles()->setParticleSetLabel( 0 );

   this->fluidSets[ 1 ]->initialize( 0, numberOfAllocatedParticles, searchRadius, dimensions, domainOrigin );
   if constexpr( ParticlesType::specifySearchedSetExplicitly() == true )
      this->fluidSets[ 1 ]->getParticles()->setParticleSetLabel( 0 );

   this->boundarySets[ 0 ]->initialize(
      numberOfBoundaryParticles, numberOfAllocatedBoundaryParticles, searchRadius, dimensions, domainOrigin );
   if constexpr( ParticlesType::specifySearchedSetExplicitly() == true )
      this->boundarySets[ 0 ]->getParticles()->setParticleSetLabel( 1 );

   this->boundarySets[ 1 ]->initialize(
      numberOfBoundaryParticles, numberOfAllocatedBoundaryParticles, searchRadius, dimensions, domainOrigin );
   if constexpr( ParticlesType::specifySearchedSetExplicitly() == true )
      this->boundarySets[ 1 ]->getParticles()->setParticleSetLabel( 1 );

   const std::string fluidFileName = params.template getParameter< std::string >( "fluidParticles" );
   const std::string boundaryFileName = params.template getParameter< std::string >( "boundaryParticles" );
   log.writeParameter( "Reading fluid particles:", fluidFileName );
   this->fluidSets[ 0 ]->template readParticlesAndVariables< typename BaseType::Reader >( fluidFileName );
   log.writeParameter( "Reading boundary particles:", boundaryFileName );
   this->boundarySets[ 0 ]->template readParticlesAndVariables< typename BaseType::Reader >( boundaryFileName );
   this->boundarySets[ 1 ]->template readParticlesAndVariables< typename BaseType::Reader >( boundaryFileName );

   this->modelParams.init( params );

   this->timeStepping.setTimeStep( params.template getParameter< RealType >( "initialTimeStep" ) );
   this->timeStepping.setEndTime( params.template getParameter< RealType >( "finalTime" ) );
   this->timeStepping.addOutputTimer( "save_results", params.template getParameter< RealType >( "snapshotPeriod" ) );

   this->caseName = params.template getParameter< std::string >( "caseName" );
   this->verbose = params.template getParameter< std::string >( "verboseIntensity" );
   this->outputDirectory = params.template getParameter< std::string >( "outputDirectory" );
   this->particlesFormat = params.template getParameter< std::string >( "particlesFormat" );

   this->timeMeasurement.addTimer( "remesh" );

   log.writeSeparator();
   const bool hasMeasuretoolFile = params.template getParameter< std::string >( "measuretoolConfig" ) != "";
   const bool hasMeasuretoolInline = params.template getParameter< std::string >( "measuretool" ) != "";
   if( hasMeasuretoolFile || hasMeasuretoolInline ) {
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

}  // namespace TNL::SPH
