#include "SolverMultiSetBase.h"

#include <SPH/configInit.h>
#include <SPH/configSetup.h>
#include "../distributedUtils.h"
#include "../shared/removeParticlesOutOfDensityLimits.h"
#include <TNL/Particles/Writers/writeBackgroundGrid.h>

namespace TNL {
namespace SPH {

template< typename Model >
void
SolverMultiSetBase< Model >::initOpenBoundaryPatches( TNL::Config::ParameterContainer& parameters, TNL::Logger& logger )
{
   logger.writeParameter( "Initialization of open boundary patches.", "" );
   const int numberOfBoundaryPatches = parameters.getParameter< int >( "openBoundaryPatches" );
   const std::string openBoundaryConfigPath = parameters.getParameter< std::string >( "open-boundary-config" );

   for( int i = 0; i < numberOfBoundaryPatches; i++ ) {
      std::string prefix = "buffer-" + std::to_string( i + 1 ) + "-";
      configSetupOpenBoundaryModelPatch< SPHConfig >( configOpenBoundary, prefix );
   }
   parseOpenBoundaryConfig( openBoundaryConfigPath, parametersOpenBoundary, configOpenBoundary, logger );

   const VectorType domainOrigin = parameters.getXyz< VectorType >( "domainOrigin" );
   const VectorType domainSize = parameters.getXyz< VectorType >( "domainSize" );
   const RealType searchRadius = parameters.getParameter< RealType >( "searchRadius" );
   const IndexVectorType gridSize = TNL::ceil( ( domainSize - domainOrigin ) / searchRadius );

   openBoundaryPatches.resize( numberOfBoundaryPatches );
   for( int i = 0; i < numberOfBoundaryPatches; i++ ) {
      std::string prefix = "buffer-" + std::to_string( i + 1 ) + "-";
      openBoundaryPatches[ i ]->config.init( parameters, parametersOpenBoundary, prefix );
      openBoundaryPatches[ i ]->initialize( parametersOpenBoundary.getParameter< int >( prefix + "numberOfParticles" ),
                                            parametersOpenBoundary.getParameter< int >( prefix + "numberOfAllocatedParticles" ),
                                            searchRadius,
                                            gridSize,
                                            domainOrigin );
   }
   logger.writeParameter( "Initialization of open boundary patches.", "Done." );
}

template< typename Model >
void
SolverMultiSetBase< Model >::initPeriodicBoundaryPatches( TNL::Config::ParameterContainer& parameters, TNL::Logger& logger )
{
   logger.writeParameter( "Initialization of periodic boundary patches.", "" );
   const int numberOfBoundaryPatches = parameters.getParameter< int >( "periodicBoundaryPatches" );
   const std::string periodicBoundaryConfigPath = parameters.getParameter< std::string >( "periodic-boundary-config" );

   TNL::Config::ConfigDescription configPeriodicBoundary;
   TNL::Config::ParameterContainer parametersPeriodicBoundary;

   for( int i = 0; i < numberOfBoundaryPatches; i++ ) {
      std::string prefix = "buffer-" + std::to_string( i + 1 ) + "-";
      configSetupOpenBoundaryModelPatch< SPHConfig >( configPeriodicBoundary, prefix );
   }
   parseOpenBoundaryConfig( periodicBoundaryConfigPath, parametersPeriodicBoundary, configPeriodicBoundary, logger );

   fluidSets[ 0 ]->initializePeriodicity( parameters, parametersPeriodicBoundary );
   boundarySets[ 0 ]->initializePeriodicity( parameters, parametersPeriodicBoundary );
}

template< typename Model >
void
SolverMultiSetBase< Model >::performNeighborSearch( bool performBoundarySearch )
{
   timeMeasurement.start( "search" );

   if constexpr( ParticlesType::specifySearchedSetExplicitly() == false ){
      for( int i = 0; i < numberOfSubsets; i++ ){
         fluidSets[ i ]->searchForNeighbors();
         if( verbose == "full" )
            logger.writeParameter( "Fluid " + std::to_string( i ) + " search procedure:", "Done." );

         if( timeStepping.getStep() == 0 || performBoundarySearch == true ){
            boundarySets[ i ]->searchForNeighbors();
            if( verbose == "full" )
               logger.writeParameter( "Boundary " + std::to_string( i ) + " search procedure:", "Done." );
         }
      }

      if( openBoundaryPatches.size() > 0 )
         for( auto& openBoundaryPatch : openBoundaryPatches ){
            openBoundaryPatch->searchForNeighbors();
            if( verbose == "full" )
               logger.writeParameter( "Open boundary patch search procedure:", "Done." );
         }
   }
   else if constexpr( ParticlesType::specifySearchedSetExplicitly() == true ){

      for( int i = 0; i < numberOfSubsets; i++ ){
         fluidSets[ i ]->makeSetSearchable();
         if( verbose == "full" )
            logger.writeParameter( "Fluid " + std::to_string( i ) + "-fluid search procedure:", "Done." );

         if( timeStepping.getStep() == 0 || performBoundarySearch == true ){
            boundarySets[ i ]->makeSetSearchable();
            if( verbose == "full" )
               logger.writeParameter( "Boundary " + std::to_string( i ) + "-boundary search procedure:", "Done." );
         }

         fluidSets[ i ]->getParticles()->addToParticleList( fluidSets[ i ]->getParticles() );
         if( verbose == "full" )
            logger.writeParameter( "Fluid " + std::to_string( i ) + "-fluid search procedure:", "Done." );

         fluidSets[ i ]->getParticles()->addToParticleList( boundarySets[ i ]->getParticles() );
         if( verbose == "full" )
            logger.writeParameter( "Fluid " + std::to_string( i ) + "-boundary search procedure:", "Done." );

         boundarySets[ i ]->getParticles()->addToParticleList( fluidSets[ i ]->getParticles() );
         if( verbose == "full" )
            logger.writeParameter( "Boundary " + std::to_string( i ) + "-fluid search procedure:", "Done." );
      }
   }

   timeMeasurement.stop( "search" );
   writeLog( "Search...", "Done." );
}

template< typename Model >
template< typename ParticleSetPointer >
void
SolverMultiSetBase< Model >::performNeighborSearchForObject( ParticleSetPointer& objectPointer )
{
   objectPointer->getParticles()->resetListWithIndices();
   objectPointer->getParticles()->computeParticleCellIndices();
   objectPointer->sortParticles();
   objectPointer->getParticles()->particlesToCells();
}

template< typename Model >
void
SolverMultiSetBase< Model >::removeParticlesOutOfDomain()
{
   for( int i = 0; i < numberOfSubsets; i++ ){
      const int numberOfParticlesToRemove = fluidSets[ i ]->getParticles()->getNumberOfParticlesToRemove();
      fluidSets[ i ]->getParticles()->removeParitclesOutOfDomain();

      if( fluidSets[ i ]->getParticles()->getNumberOfParticlesToRemove() > numberOfParticlesToRemove ){
         const int numberOfParticlesOutOfDomain = fluidSets[ i ]->getParticles()->getNumberOfParticlesToRemove() - numberOfParticlesToRemove;
         this->totalNumberOfParticlesOutOfDomain += numberOfParticlesOutOfDomain;
         logger.writeParameter( "Particles out of domain removed (set " + std::to_string( i ) + "):", numberOfParticlesOutOfDomain );
         logger.writeParameter( "Total particles out of domain removed:", this->totalNumberOfParticlesOutOfDomain );
      }
   }
}

template< typename Model >
void
SolverMultiSetBase< Model >::removeParticlesOutOfDensityLimits()
{
   const int numberOfParticlesOutOfDensityLimits = customFunctions::removeParticlesOutOfDensityLimits( fluidSets[ 0 ], modelParams );
   if( numberOfParticlesOutOfDensityLimits > 0 ){
      this->totalNumberOfParticlesOutOfDensityLimits += numberOfParticlesOutOfDensityLimits;
      logger.writeParameter( "Particles out of density limits:", numberOfParticlesOutOfDensityLimits );
      logger.writeParameter( "Total particles out of density limits:", this->totalNumberOfParticlesOutOfDensityLimits );
   }
}

template< typename Model >
void
SolverMultiSetBase< Model >::extrapolateOpenBC()
{
   timeMeasurement.start( "extrapolate-openbc" );
   for( long unsigned int i = 0; i < openBoundaryPatches.size(); i++ ) {
      model.extrapolateOpenBoundaryData( fluidSets[ 0 ], openBoundaryPatches[ i ], modelParams, openBoundaryPatches[ i ]->config );
   }
   timeMeasurement.stop( "extrapolate-openbc" );
   writeLog( "Extrapolate open BC...", "Done." );
}

template< typename Model >
void
SolverMultiSetBase< Model >::applyOpenBC( const RealType timeStepFact )
{
   timeMeasurement.start( "apply-openbc" );
   for( long unsigned int i = 0; i < openBoundaryPatches.size(); i++ ) {
      openBoundaryModel.applyOpenBoundary(
         timeStepFact * timeStepping.getTimeStep(), fluidSets[ 0 ], openBoundaryPatches[ i ], openBoundaryPatches[ i ]->config );
   }
   timeMeasurement.stop( "apply-openbc" );
   writeLog( "Update open BC...", "Done." );
}

template< typename Model >
void
SolverMultiSetBase< Model >::applyPeriodicBCEnforce()
{
   timeMeasurement.start( "enforce-periodic-bc" );
   timeMeasurement.start( "periodicity-fluid-updateZone" );
   fluidSets[ 0 ]->enforcePeriodicPatches();
   timeMeasurement.stop( "periodicity-boundary-updateZone" );
   timeMeasurement.start( "periodicity-fluid-updateZone" );
   boundarySets[ 0 ]->enforcePeriodicPatches();
   timeMeasurement.stop( "periodicity-boundary-updateZone" );
   timeMeasurement.stop( "enforce-periodic-bc" );
   writeLog( "Apply periodic BC...", "Done." );
}

template< typename Model >
void
SolverMultiSetBase< Model >::applyPeriodicBCTransfer()
{
   timeMeasurement.start( "transfer-periodic-bc" );
   const long unsigned int numberOfPeriodicPatches = fluidSets[ 0 ]->periodicPatches.size();
   for( long unsigned int i = 0; i < numberOfPeriodicPatches; i++ ) {
      openBoundaryModel.periodicityParticleTransfer( fluidSets[ 0 ], fluidSets[ 0 ]->periodicPatches[ i ] );
   }
   timeMeasurement.stop( "transfer-periodic-bc" );
   writeLog( "Transfer periodic BC...", "Done." );
}

template< typename Model >
void
SolverMultiSetBase< Model >::computeTimeStep()
{
   timeStepping.computeTimeStep( fluidSets[ 0 ], modelParams );
}

template< typename Model >
void
SolverMultiSetBase< Model >::updateTime()
{
   timeStepping.updateTimeStep();
}

template< typename Model >
void
SolverMultiSetBase< Model >::measure()
{
   //FIXME: Needs generalization for multi-set simulations
   simulationMonitor.template measure< typename ModelParams::KernelFunction, typename ModelParams::EOS >(
         fluidSets[ 0 ], boundarySets[ 0 ], modelParams, timeStepping, logger, verbose );
}

template< typename Model >
template< typename Stage >
void
SolverMultiSetBase< Model >::integrate( const Stage integrationStage, const bool integrateBoundary )
{
   timeMeasurement.start( "integrate" );
   for( int i = 0; i < numberOfSubsets; i++ )
      integrator->integratStepVerlet( fluidSets[ i ], boundarySets[ i ], timeStepping, integrateBoundary );
   timeMeasurement.stop( "integrate" );
   writeLog( "Integrate...", "Done." );
}

template< typename Model >
void
SolverMultiSetBase< Model >::integrateVerletStep( const bool integrateBoundary )
{
   timeMeasurement.start( "integrate" );
   for( int i = 0; i < numberOfSubsets; i++ )
      integrator->integratStepVerlet( fluidSets[ i ], boundarySets[ i ], timeStepping, integrateBoundary );
   timeMeasurement.stop( "integrate" );
   writeLog( "Integrate...", "Done." );
}

template< typename Model >
void
SolverMultiSetBase< Model >::symplecticVerletPredictor()
{
   timeMeasurement.start( "integrate" );
   for( int i = 0; i < numberOfSubsets; i++ )
      integrator->integratePredictorStep( fluidSets[ i ], boundarySets[ i ], timeStepping );
   timeMeasurement.stop( "integrate" );
   writeLog( "Integrate - predictor step...", "Done." );
}

template< typename Model >
void
SolverMultiSetBase< Model >::symplecticVerletCorrector()
{
   timeMeasurement.start( "integrate" );
   for( int i = 0; i < numberOfSubsets; i++ )
      integrator->integrateCorrectorStep( fluidSets[ i ], boundarySets[ i ], timeStepping );
   timeMeasurement.stop( "integrate" );
   writeLog( "Integrate - corrector step...", "Done." );
}

template< typename Model >
void
SolverMultiSetBase< Model >::midpointPredictor()
{
   timeMeasurement.start( "integrate" );
   for( int i = 0; i < numberOfSubsets; i++ )
      integrator->predictor( timeStepping.getTimeStep(), fluidSets[ i ] );
   timeMeasurement.stop( "integrate" );
   writeLog( "Integrate: predictor...", "Done." );
}

template< typename Model >
void
SolverMultiSetBase< Model >::midpointUpdateVariables()
{
   timeMeasurement.start( "integrate" );
   for( int i = 0; i < numberOfSubsets; i++ )
      integrator->midpointUpdateVariables( timeStepping.getTimeStep(), fluidSets[ i ] );
   timeMeasurement.stop( "integrate" );
   writeLog( "Integrate: midpoint update...", "Done." );
}

template< typename Model >
void
SolverMultiSetBase< Model >::midpointResidualsAndRelaxationFactor()
{
   for( int i = 0; i < numberOfSubsets; i++ ){
      integrator->computeResiduals( fluidSets[ i ], modelParams );
      integrator->updateRelaxationFactor( modelParams );
   }
}

template< typename Model >
void
SolverMultiSetBase< Model >::midpointRelax()
{
   timeMeasurement.start( "integrate" );
   for( int i = 0; i < numberOfSubsets; i++ )
      integrator->relax( fluidSets[ i ], modelParams );
   timeMeasurement.stop( "integrate" );
   writeLog( "Integrate: relax...", "Done." );
}

template< typename Model >
void
SolverMultiSetBase< Model >::midpointCorrector()
{
   timeMeasurement.start( "integrate" );
   for( int i = 0; i < numberOfSubsets; i++ ){
      if( timeStepping.getStep() == 0 )
         integrator->corrector( 0, fluidSets[ i ] );
      else
         integrator->corrector( timeStepping.getTimeStep(), fluidSets[ i ] );
   }
   timeMeasurement.stop( "integrate" );
   writeLog( "Integrate: corrector...", "Done." );
}

template< typename Model >
void
SolverMultiSetBase< Model >::save( bool writeParticleCellIndex )
{
   if( ( verbose == "with-snapshot" ) || ( verbose == "full" ) )
      writeInfo();

   const int step = timeStepping.getStep();
   const RealType time = timeStepping.getTime();

   for( int i = 0; i < numberOfSubsets; i++ ){
#ifdef HAVE_MPI
      std::string outputFileNameFluid = outputDirectory + "/fluid_rank" + std::to_string( TNL::MPI::GetRank() ) + "_" + std::to_string( time ) + "_particles.vtk";
#else
      std::string outputFileNameFluid;
      if( numberOfSubsets == 1 )
         outputFileNameFluid = outputDirectory + "/fluid_" + std::to_string( time ) + "_particles.vtk";
      else
         outputFileNameFluid = outputDirectory + "/fluid_subdomain" + std::to_string( i ) + "_" + std::to_string( time ) + "_particles.vtk";
#endif
      std::string tmpOutputFileNameFluid = outputFileNameFluid + ".tmp";
      fluidSets[ i ]->template writeParticlesAndVariables< Writer >( tmpOutputFileNameFluid, writeParticleCellIndex );
      std::filesystem::rename( tmpOutputFileNameFluid, outputFileNameFluid );
      logger.writeParameter( "Saved:", outputFileNameFluid );

#ifdef HAVE_MPI
      std::string outputFileNameBound = outputDirectory + "/boundary_rank" + std::to_string( TNL::MPI::GetRank() ) + "_" + std::to_string( time ) + "_particles.vtk";
#else
      std::string outputFileNameBound;
      if( numberOfSubsets == 1 )
         outputFileNameBound = outputDirectory + "/boundary_" + std::to_string( time ) + "_particles.vtk";
      else
         outputFileNameBound = outputDirectory + "/boundary_subdomain" + std::to_string( i ) + "_" + std::to_string( time ) + "_particles.vtk";
#endif
      std::string tmpOutputFileNameBound = outputFileNameBound + ".tmp";
      boundarySets[ i ]->template writeParticlesAndVariables< Writer >( tmpOutputFileNameBound, writeParticleCellIndex );
      std::filesystem::rename( tmpOutputFileNameBound, outputFileNameBound );
      logger.writeParameter( "Saved:", outputFileNameBound );
   }

   if( openBoundaryPatches.size() ) {
      for( auto& openBoundaryPatch : openBoundaryPatches ) {
         std::string outputFileNameOpenBound =
            outputDirectory + "/" + openBoundaryPatch->parameters.identifier + "_" + std::to_string( time ) + "_particles.vtk";
         std::string tmpOutputFileNameOpenBound = outputFileNameOpenBound + ".tmp";
         openBoundaryPatch->template writeParticlesAndVariables< Writer >( tmpOutputFileNameOpenBound, writeParticleCellIndex );
         std::filesystem::rename( tmpOutputFileNameOpenBound, outputFileNameOpenBound );
         logger.writeParameter( "Saved:", outputFileNameOpenBound );
      }
   }

#ifdef HAVE_MPI
   std::string outputFileNameGrid = outputDirectory + "/grid_rank" + std::to_string( TNL::MPI::GetRank() + 1 ) + "_" + std::to_string( time ) + ".vtk";
#else
   std::string outputFileNameGrid;
   if( numberOfSubsets == 1 )
      outputFileNameGrid = outputDirectory + "/grid_" + std::to_string( time ) + ".vtk";
   else
      outputFileNameGrid = outputDirectory + "/grid_subdomain" + std::to_string( 0 ) + "_" + std::to_string( time ) + ".vtk";
#endif
   TNL::Writers::writeBackgroundGrid( outputFileNameGrid, fluidSets[ 0 ]->getParticles()->getGridDimensions(), fluidSets[ 0 ]->getParticles()->getGridOrigin(), fluidSets[ 0 ]->getParticles()->getSearchRadius() );
   logger.writeParameter( "Saved:", outputFileNameGrid );

   simulationMonitor.save( logger );
}

template< typename Model >
void
SolverMultiSetBase< Model >::makeSnapshot()
{
   const bool savePressure = true;
   if( timeStepping.checkOutputTimer( "save_results" ) ){
      if( savePressure ){
         for( int i = 0; i < numberOfSubsets; i++ ){
            model.computePressureFromDensity( fluidSets[ i ], modelParams );
            model.computePressureFromDensity( boundarySets[ i ], modelParams );
         }
      }
      save();
      writeLog( "Save results...", "Done." );
   }
}

template< typename Model >
template< typename Func >
void
SolverMultiSetBase< Model >::initUserConfig( Func&& userConfigFunction )
{
   const std::string userConfigPath = parameters.getParameter< std::string >( "user-defined-config" );
   if( userConfigPath == "" )
      return;

   userConfigFunction( userConfig );
   parseUserDefinedConfig( userConfigPath, userParams, userConfig, logger );
}

template< typename Model >
void
SolverMultiSetBase< Model >::writeProlog( bool writeSystemInformation ) noexcept
{
   logger.writeHeader( "SPH simulation configuration." );
   logger.writeParameter( "Case name:", caseName );
   logger.writeSeparator();

   const bool printGPUs = std::is_same< DeviceType, TNL::Devices::Cuda >::value;

   if( TNL::MPI::isInitialized() )
      logger.writeParameter( "MPI processes:", TNL::MPI::GetSize() );
   logger.writeParameter( "Device type:", TNL::getType< DeviceType >() );
   if( ! printGPUs ) {
      if( TNL::Devices::Host::isOMPEnabled() ) {
         logger.writeParameter( "OMP enabled:", "yes", 1 );
         logger.writeParameter( "OMP threads:", TNL::Devices::Host::getMaxThreadsCount(), 1 );
      }
      else
         logger.writeParameter( "OMP enabled:", "no", 1 );
   }

   logger.writeParameter( "Particle system type:", Model::ParticlesType::writeModelType() );
   logger.writeParameter( "SPH model:", Model::writeModelType() );
   logger.writeParameter( "Verbose:", verbose );
   logger.writeParameter( "Output directory:", outputDirectory );
   logger.writeParameter( "Particles format", particlesFormat );
#ifdef HAVE_MPI
   if( TNL::MPI::isInitialized() ){
      logger.writeParameter( "Load balancing step interval: ", loadBalancingStepInterval );
      logger.writeParameter( "Load balancing measure:", loadBalancingMeasure );
      if( loadBalancingMeasure == "computationalTime" )
         logger.writeParameter( "Comp. time fraction difference to balance [-]:",
             fluidSets[ 0 ]->getDistributedParticles()->getCompTimeResizePercentageTrashold() );
      else if( loadBalancingMeasure == "numberOfParticles" )
         logger.writeParameter( "Particles count fraction difference to balance [-]:",
             fluidSets[ 0 ]->getDistributedParticles()->getParticlesCountResizeTrashold() );
   }
#endif
   writePrologModel( logger, modelParams );

   for( int i = 0; i < numberOfSubsets; i++ ){
      logger.writeHeader( "Fluid " + std::to_string( i ) + " object information." );
      fluidSets[ i ]->writeProlog( logger );
      logger.writeHeader( "Boundary " + std::to_string( i ) + " object information:" );
      boundarySets[ i ]->writeProlog( logger );
   }

   if( openBoundaryPatches.size() > 0 ) {
      for( long unsigned int i = 0; i < openBoundaryPatches.size(); i++ ) {
         logger.writeHeader( "Open boundary buffer" + std::to_string( i + 1 ) + "." );
         openBoundaryPatches[ i ]->writeProlog( logger );
         logger.writeSeparator();
         openBoundaryPatches[ i ]->config.writeProlog( logger );
      }
   }

   if( fluidSets[ 0 ]->periodicPatches.size() > 0 ) {
      const long unsigned int numberOfPeriodicPatches = fluidSets[ 0 ]->periodicPatches.size();
      for( long unsigned int i = 0; i < numberOfPeriodicPatches; i++ ) {
         logger.writeHeader( "Periodic boundary patch " + std::to_string( i + 1 ) + "." );
         fluidSets[ 0 ]->periodicPatches[ i ]->writeProlog( logger );
      }
   }

   logger.writeHeader( "System information." );
   if( writeSystemInformation ) {
      logger.writeSystemInformation( printGPUs );
      logger.writeSeparator();
      logger.writeCurrentTime( "Started at:" );
      logger.writeSeparator();
   }
}

template< typename Model >
template< typename ParameterType >
void
SolverMultiSetBase< Model >::writeLog( const std::string& label, const ParameterType& value, int parameterLevel )
{
   if( verbose == "full" )
      logger.writeParameter( label, value, parameterLevel );
}

template< typename Model >
void
SolverMultiSetBase< Model >::writeInfo() noexcept
{
   logger.writeSeparator();
   logger.writeParameter( "Simulation time: " + std::to_string( timeStepping.getTime() )
                             + " s, simulation step: " + std::to_string( timeStepping.getStep() ),
                          "" );
   logger.writeCurrentTime( "Current time:" );
   for( int i = 0; i < numberOfSubsets; i++ ){
      logger.writeParameter( "Number of fluid particles:", fluidSets[ i ]->getNumberOfParticles() );
      logger.writeParameter( "Number of allocated fluid particles:", fluidSets[ i ]->getNumberOfAllocatedParticles() );
      logger.writeParameter( "Number of boundary particles:", boundarySets[ i ]->getNumberOfParticles() );
      logger.writeParameter( "Number of allocated boundary particles:", boundarySets[ i ]->getNumberOfAllocatedParticles() );
   }
   if( openBoundaryPatches.size() > 0 ) {
      for( long unsigned int i = 0; i < openBoundaryPatches.size(); i++ )
         logger.writeParameter( "Number of buffer" + std::to_string( i + 1 ) + " particles:",
                                openBoundaryPatches[ i ]->getNumberOfParticles() );
   }
   if( fluidSets[ 0 ]->periodicPatches.size() > 0 ) {
      if( verbose == "full" ) {
         for( long unsigned int i = 0; i < fluidSets[ 0 ]->periodicPatches.size(); i++ )
            logger.writeParameter( "Number of fluid particles in periodic patch " + std::to_string( i + 1 ) + ": ",
                                   fluidSets[ 0 ]->periodicPatches[ i ]->particleZone.getNumberOfParticles() );
         for( long unsigned int i = 0; i < boundarySets[ 0 ]->periodicPatches.size(); i++ )
            logger.writeParameter( "Number of boundary particles in periodic patch " + std::to_string( i + 1 ) + " :",
                                   boundarySets[ 0 ]->periodicPatches[ i ]->particleZone.getNumberOfParticles() );
      }
   }
   logger.writeSeparator();
}

template< typename Model >
void
SolverMultiSetBase< Model >::writeEpilog() noexcept
{
   logger.writeHeader( "SPH simulation successfully finished." );
   logger.writeCurrentTime( "Ended at:" );
   timeMeasurement.writeInfo( logger, timeStepping.getStep() );

   std::string saveTimersOutputName = outputDirectory + "/time_measurements";
   timeMeasurement.writeInfoToJson( saveTimersOutputName, timeStepping.getStep() );
}

#ifdef HAVE_MPI

template< typename Model >
void
SolverMultiSetBase< Model >::synchronizeDistributedSimulation()
{
   writeLog( "Starting synchronization.", "" );
   timeMeasurement.start( "synchronize" );
   for( int i = 0; i < numberOfSubsets; i++ ){
      fluidSets[ i ]->synchronizeObject();
      boundarySets[ i ]->synchronizeObject();
   }
   timeMeasurement.stop( "synchronize" );
   writeLog( "Synchronize...", "Done." );
}

template< typename Model >
void
SolverMultiSetBase< Model >::resetOverlaps()
{
   for( int i = 0; i < numberOfSubsets; i++ ){
      fluidSets[ i ]->getParticles()->removeParitclesOutOfDomain();
      boundarySets[ i ]->getParticles()->removeParitclesOutOfDomain();
   }
}

template< typename Model >
void
SolverMultiSetBase< Model >::performLoadBalancing()
{
   fluidSets[ 0 ]->getDistributedParticles()->setNumberOfParticlesForLoadBalancing( fluidSets[ 0 ]->getNumberOfParticles() );
   fluidSets[ 0 ]->getDistributedParticles()->setCompTimeForLoadBalancing( timeMeasurement.getTotalTime() - subdomainCompTimeBackup );
   fluidSets[ 0 ]->synchronizeBalancingMeasures();

   std::pair< IndexVectorType, VectorType > subdomainAdjustment;
   if( loadBalancingMeasure == "computationalTime" )
      subdomainAdjustment = fluidSets[ 0 ]->getDistributedParticles()->loadBalancingDomainAdjustmentCompTime();
   else if( loadBalancingMeasure == "numberOfParticles" )
      subdomainAdjustment = fluidSets[ 0 ]->getDistributedParticles()->loadBalancingDomainAdjustment();
   else
      std::cerr << "Invalid load balancing metrics. Load balancing metrics is: " << loadBalancingMeasure << "." << std::endl;

   const IndexVectorType gridDimensionsAdjustment = subdomainAdjustment.first;
   const VectorType gridOriginAdjustment = subdomainAdjustment.second * fluidSets[ 0 ]->getParticles()->getSearchRadius();

   const IndexVectorType updatedGridDimensions = fluidSets[ 0 ]->getParticles()->getGridDimensions() + gridDimensionsAdjustment;
   const VectorType updatedGridOrigin = fluidSets[ 0 ]->getParticles()->getGridOrigin() + gridOriginAdjustment;

   logger.writeParameter( "Load balancing - subdomain adjustment: ", "" );
   logger.writeParameter( "Grid dimensions adjustment: ", gridDimensionsAdjustment );
   logger.writeParameter( "Grid origin adjustment: ", gridOriginAdjustment );
   logger.writeParameter( "Old grid dimensions: ", fluidSets[ 0 ]->getParticles()->getGridDimensions() );
   logger.writeParameter( "Old grid origin adjustment: ", fluidSets[ 0 ]->getParticles()->getGridOrigin() );
   logger.writeParameter( "Old firstLastCellParticleList size: ", fluidSets[ 0 ]->getParticles()->getCellFirstLastParticleList().getSize() );

   for( int i = 0; i < numberOfSubsets; i++ ){
      fluidSets[ i ]->getParticles()->setGridDimensions( updatedGridDimensions );
      fluidSets[ i ]->getParticles()->setGridOrigin( updatedGridOrigin );
      boundarySets[ i ]->getParticles()->setGridDimensions( updatedGridDimensions );
      boundarySets[ i ]->getParticles()->setGridOrigin( updatedGridOrigin );
   }

   const IndexVectorType updatedGridOriginGlobalCoords = fluidSets[ 0 ]->getParticles()->getGridOriginGlobalCoords() + subdomainAdjustment.second;
   for( int i = 0; i < numberOfSubsets; i++ ){
      fluidSets[ i ]->getParticles()->setGridOriginGlobalCoords( updatedGridOriginGlobalCoords );
      boundarySets[ i ]->getParticles()->setGridOriginGlobalCoords( updatedGridOriginGlobalCoords );
   }

   logger.writeParameter( "New grid dimensions: ", fluidSets[ 0 ]->getParticles()->getGridDimensions() );
   logger.writeParameter( "New grid origin adjustment: ", fluidSets[ 0 ]->getParticles()->getGridOrigin() );
   logger.writeParameter( "New firstLastCellParticleList size: ", fluidSets[ 0 ]->getParticles()->getCellFirstLastParticleList().getSize() );

   for( int i = 0; i < numberOfSubsets; i++ ){
      fluidSets[ i ]->getDistributedParticles()->updateDistriutedGridParameters( updatedGridDimensions,
                                                                                  updatedGridOrigin,
                                                                                  1,
                                                                                  fluidSets[ i ]->getParticles()->getSearchRadius() );
      boundarySets[ i ]->getDistributedParticles()->updateDistriutedGridParameters( updatedGridDimensions,
                                                                                     updatedGridOrigin,
                                                                                     1,
                                                                                     boundarySets[ i ]->getParticles()->getSearchRadius() );
   }
   writeLoadBalancingInfo( gridDimensionsAdjustment[ 0 ] );
   this->subdomainCompTimeBackup = timeMeasurement.getTotalTime();
}

template< typename Model >
void
SolverMultiSetBase< Model >::writeLoadBalancingInfo( const int gridResize )
{
   const std::string outputPath = outputDirectory + "/loadBalancingMetrics_rank" + std::to_string( TNL::MPI::GetRank() ) + ".dat";
   std::ofstream outfile;
   outfile.open( outputPath, std::ios_base::app );
   outfile << timeStepping.getStep() << " "
           << timeStepping.getTime() << " "
           << fluidSets[ 0 ]->getNumberOfParticles() << " "
           << fluidSets[ 0 ]->getNumberOfAllocatedParticles() << " "
           << boundarySets[ 0 ]->getNumberOfParticles() << " "
           << boundarySets[ 0 ]->getNumberOfAllocatedParticles() << " "
           << timeMeasurement.getTotalTime() - subdomainCompTimeBackup << " "
           << fluidSets[ 0 ]->getParticles()->getCellFirstLastParticleList().getSize() << " "
           << gridResize << std::endl;
}

template< typename Model >
void
SolverMultiSetBase< Model >::balanceSubdomains()
{
   if( ( timeStepping.getStep() % loadBalancingStepInterval == 0 ) && ( timeStepping.getStep() > 1 ) ){

      performNeighborSearch( true );

      logger.writeSeparator();
      writeLog( "Starting load balancing.", "" );
      timeMeasurement.start( "rebalance" );
      performLoadBalancing();
      timeMeasurement.stop( "rebalance" );
      writeLog( "Load balancing...", "Done." );
      logger.writeSeparator();
      TNL::MPI::Barrier( communicator );

      resetOverlaps();
      writeLog( "Reset overlaps...", "Done." );
      TNL::MPI::Barrier( communicator );

      performNeighborSearch( true );
      TNL::MPI::Barrier( communicator );

      synchronizeDistributedSimulation();
      TNL::MPI::Barrier( communicator );
   }
}

#endif

} // namespace SPH
} // namespace TNL
