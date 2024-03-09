#include "SPHMultiset_CFD.h"
#include "SPH/OpenBoundaryConfig.h"
#include "SPH/TimeMeasurement.h"
#include <string>

namespace TNL {
namespace SPH {

template< typename Model >
void
SPHMultiset_CFD< Model >::init( TNL::Config::ParameterContainer& parameters, TNL::Logger& logger  )
{
   logger.writeHeader( "SPH simulation initialization." );

   // compute domain properetis
   const VectorType domainOrigin = parameters.getXyz< VectorType >( "domainOrigin" );
   const VectorType domainSize = parameters.getXyz< VectorType >( "domainSize" );
   const RealType searchRadius = parameters.getParameter< RealType >( "searchRadius" );
   const IndexVectorType gridSize = TNL::ceil( ( domainSize - domainOrigin ) / searchRadius );

   // init fluid
   fluid->initialize( parameters.getParameter< int >( "numberOfParticles" ),
                      parameters.getParameter< int >( "numberOfAllocatedParticles" ),
                      searchRadius,
                      gridSize,
                      domainOrigin );

   // init boundary
   boundary->initialize( parameters.getParameter< int >( "numberOfBoundaryParticles" ),
                         parameters.getParameter< int >( "numberOfAllocatedBoundaryParticles" ),
                         searchRadius,
                         gridSize,
                         domainOrigin );

   // init open boundary patches
   const int numberOfBoundaryPatches = parameters.getParameter< int >( "openBoundaryPatches" );
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ) { //TODO: I dont like this.
      openBoundaryPatches.resize( numberOfBoundaryPatches );
      for( int i = 0; i < numberOfBoundaryPatches; i++ ) {
         std::string prefix = "buffer-" + std::to_string( i + 1 ) + "-";
         openBoundaryPatches[ i ]->config.init( parameters, prefix );
         openBoundaryPatches[ i ]->initialize( parameters.getParameter< int >( prefix + "numberOfParticles" ),
                                               parameters.getParameter< int >( prefix + "numberOfAllocatedParticles" ),
                                               searchRadius,
                                               gridSize,
                                               domainOrigin );
      }
   }

   // init periodic boundary conditions
   const int numberOfPeriodicPatches = parameters.getParameter< int >( "periodicBoundaryPatches" );
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfPeriodicBuffers > 0 ) {//TODO: I dont like this.
      periodicBoundaryPatches.resize( numberOfPeriodicPatches );
      for( int i = 0; i < numberOfPeriodicPatches; i++ ) {
         std::string prefix = "buffer-" + std::to_string( i + 1 ) + "-";
         periodicBoundaryPatches[ i ]->config.init( parameters, prefix );
         periodicBoundaryPatches[ i ]->initialize( parameters,
                                                   prefix,
                                                   parameters.getParameter< int >( prefix + "numberOfParticles" ),
                                                   parameters.getParameter< int >( prefix + "numberOfAllocatedParticles" ),
                                                   searchRadius,
                                                   gridSize,
                                                   domainOrigin );
      }
   }

   // init model parameters
   modelParams.init( parameters );

   // init time stepping
   timeStepping.setTimeStep( parameters.getParameter< RealType >( "initial-time-step" ) );
   timeStepping.setEndTime( parameters.getParameter< RealType >( "final-time" ) );
   timeStepping.addOutputTimer( "save_results", parameters.getParameter< RealType >( "snapshot-period" ) );

   // control
   caseName = parameters.getParameter< std::string >( "case-name" );
   verbose = parameters.getParameter< std::string >( "verbose-intensity" );
   outputDirecotry = parameters.getParameter< std::string >( "output-directory" );
   particlesFormat = parameters.getParameter< std::string >( "particles-format" );

   // read particle data
   logger.writeParameter( "Reading fluid particles:", parameters.getParameter< std::string >( "fluid-particles" ) );
   fluid->template readParticlesAndVariables< SimulationReaderType>(
         parameters.getParameter< std::string >( "fluid-particles" ) );
   logger.writeParameter( "Reading boundary particles:", parameters.getParameter< std::string >( "boundary-particles" ) );
   boundary->template readParticlesAndVariables< SimulationReaderType >(
         parameters.getParameter< std::string >( "boundary-particles" ) );

   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ) { //TODO: I dont like this.
      for( int i = 0; i < numberOfBoundaryPatches; i++ ) {
         std::string prefix = "buffer-" + std::to_string( i + 1 ) + "-";
         openBoundaryPatches[ i ]->template readParticlesAndVariables< SimulationReaderType >(
               parameters.getParameter< std::string >( prefix + "particles" ) );
      }
   }

   // initialize the measuretool
   logger.writeSeparator();
   if( parameters.getParameter< std::string >( "measuretool-config" ) != "" ){
      logger.writeParameter( "Simulation monitor initialization.", "" );
      simulationMonitor.init( parameters, timeStepping, logger );
      logger.writeParameter( "Simulation monitor initialization.", "Done." );
   }

   logger.writeHeader( "SPH simulation successfully initialized." );
}

template< typename Model >
void
SPHMultiset_CFD< Model >::performNeighborSearch( TNL::Logger& logger )
{
   //reset cell indices
   timeMeasurement.start( "search_reset" );
   fluid->particles->resetListWithIndices();
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ){
      for( auto& openBoundaryPatch : openBoundaryPatches )
         openBoundaryPatch->particles->resetListWithIndices();
   }

   if( timeStepping.getStep() == 0 )
      boundary->particles->resetListWithIndices();
   timeMeasurement.stop( "search_reset" );
   writeLog( logger, "Search - reset ...", "Done." );

   //compute cell indices
   timeMeasurement.start( "search_cellIndices" );
   fluid->particles->computeParticleCellIndices();
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ){
      for( auto& openBoundaryPatch : openBoundaryPatches )
         openBoundaryPatch->particles->computeParticleCellIndices();
   }

   if( timeStepping.getStep() == 0 )
      boundary->particles->computeParticleCellIndices();
   timeMeasurement.stop( "search_cellIndices" );
   writeLog( logger, "Search - compute cell indices ...", "Done." );

   //sort particles
   timeMeasurement.start( "search_sort" );
   fluid->sortParticles();
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ){
      for( auto& openBoundaryPatch : openBoundaryPatches )
         openBoundaryPatch->sortParticles();
   }

   if( timeStepping.getStep() == 0 )
      boundary->sortParticles();
   timeMeasurement.stop( "search_sort" );
   writeLog( logger, "Search - sort ...", "Done." );

   //update number of particles TODO: Do this in elegant way.
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ){
      for( auto& openBoundaryPatch : openBoundaryPatches ){
         fluid->particles->setNumberOfParticles( fluid->particles->getNumberOfParticles() - openBoundaryPatch->numberOfFluidParticlesToRemove );
         fluid->particles->setLastActiveParticle( fluid->particles->getLastActiveParticle() - openBoundaryPatch->numberOfFluidParticlesToRemove );
         fluid->setLastActiveParticle( fluid->getLastActiveParticle() - openBoundaryPatch->numberOfFluidParticlesToRemove );
         openBoundaryPatch->numberOfFluidParticlesToRemove = 0;
      }
      writeLog( logger, "Search - resize ...", "Done." );
   }

   //assign particles to cells
   timeMeasurement.start( "search_toCells" );
   fluid->particles->particlesToCells();
   writeLog( logger, "Search - particles to cells - fluid ...", "Done." );
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ){
      for( auto& openBoundaryPatch : openBoundaryPatches )
         openBoundaryPatch->particles->particlesToCells();
   }

   if( timeStepping.getStep() == 0 )
      boundary->particles->particlesToCells();
   timeMeasurement.stop( "search_toCells" );
   writeLog( logger, "Search - particles to cells ...", "Done." );
}

template< typename Model >
template< typename ParticleSetPointer >
void
SPHMultiset_CFD< Model >::performNeighborSearchForObject( ParticleSetPointer& objectPointer )
{
   objectPointer->particles->resetListWithIndices();
   objectPointer->particles->computeParticleCellIndices();
   objectPointer->sortParticles();
   objectPointer->particles->particlesToCells();
}

template< typename Model >
void
SPHMultiset_CFD< Model >::extrapolateOpenBC()
{
   for( long unsigned int i = 0; i < std::size( openBoundaryPatches ); i++ ){
      //TODO Check if open boundary buffer is really open boundary buffer
      model.extrapolateOpenBoundaryData( fluid, openBoundaryPatches[ i ], modelParams, openBoundaryPatches[ i ]->config );
   }
}

template< typename Model >
void
SPHMultiset_CFD< Model >::applyOpenBC()
{
   for( long unsigned int i = 0; i < std::size( openBoundaryPatches ); i++ ){
      //TODO Check if open boundary buffer is really open boundary buffer
      openBoundaryModel.applyOpenBoundary( timeStepping.getTimeStep(),
                                           fluid,
                                           openBoundaryPatches[ i ],
                                           openBoundaryPatches[ i ]->config );
   }
}

template< typename Model >
void
SPHMultiset_CFD< Model >::applyPeriodicBCEnforce()
{
   for( long unsigned int i = 0; i < std::size( periodicBoundaryPatches ); i++ ){
      int pairedPeriodicBuffer = periodicBoundaryPatches[ i ]->config.pairedPeriodicBuffer - 1;
      openBoundaryModel.applyPeriodicBoundary( fluid,
                                               periodicBoundaryPatches[ i ]->fluidPeriodicPatch,
                                               periodicBoundaryPatches[ pairedPeriodicBuffer ]->fluidPeriodicPatch,
                                               periodicBoundaryPatches[ i ]->config,
                                               periodicBoundaryPatches[ pairedPeriodicBuffer ]->config );
      openBoundaryModel.applyPeriodicBoundaryOnBoundary( boundary,
                                                         periodicBoundaryPatches[ i ]->boundaryPeriodicPatch,
                                                         periodicBoundaryPatches[ pairedPeriodicBuffer ]->boundaryPeriodicPatch,
                                                         periodicBoundaryPatches[ i ]->config,
                                                         periodicBoundaryPatches[ pairedPeriodicBuffer ]->config );
   }

   //sort the buffer after updating the periodicity zones
   for( long unsigned int i = 0; i < std::size( periodicBoundaryPatches ); i++ ){
      performNeighborSearchForObject( periodicBoundaryPatches[ i ]->fluidPeriodicPatch );
      performNeighborSearchForObject( periodicBoundaryPatches[ i ]->boundaryPeriodicPatch );
   }
}

template< typename Model >
void
SPHMultiset_CFD< Model >::applyPeriodicBCTransfer()
{
   for( long unsigned int i = 0; i < std::size( periodicBoundaryPatches ); i++ ){
      openBoundaryModel.periodicityParticleTransfer( fluid,
                                                     periodicBoundaryPatches[ i ]->fluidPeriodicPatch,
                                                     periodicBoundaryPatches[ i ]->config );
   }
}

template< typename Model >
void
SPHMultiset_CFD< Model >::interact()
{
   //Interacteraction between fluid and fluid and fluid and boundary
   model.interaction( fluid, boundary, modelParams );
   model.updateSolidBoundary( fluid, boundary, modelParams );

   //Interact between fluid and open boundary patches
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ){
      for( long unsigned int i = 0; i < std::size( openBoundaryPatches ); i++ ){
         openBoundaryPatches[ i ]->zone.updateParticlesInZone( fluid->particles );
         model.interactionWithOpenBoundary( fluid, openBoundaryPatches[ i ], modelParams );
         model.updateSolidBoundaryOpenBoundary( boundary, openBoundaryPatches[ i ], modelParams );
      }
   }

   //Interact between fluid and open boundary patches
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfPeriodicBuffers > 0 ){
      for( long unsigned int i = 0; i < std::size( periodicBoundaryPatches ); i++ ){
         periodicBoundaryPatches[ i ]->fluidPeriodicPatch->zone.updateParticlesInZone( fluid->particles );
         model.interactionWithOpenBoundary( fluid, periodicBoundaryPatches[ i ]->fluidPeriodicPatch, modelParams );
         periodicBoundaryPatches[ i ]->boundaryPeriodicPatch->zone.updateParticlesInZone( fluid->particles ); //fuck this out
         model.interactionWithBoundaryPatches( fluid, periodicBoundaryPatches[ i ]->boundaryPeriodicPatch, modelParams );
         model.updateSolidBoundaryOpenBoundary( boundary, periodicBoundaryPatches[ i ]->fluidPeriodicPatch, modelParams );
      }
   }

   //Finalize the interaction
   model.finalizeInteraction( fluid, boundary, modelParams );
}

template< typename Model >
template< typename SPHKernelFunction, typename EOS >
void
SPHMultiset_CFD< Model >::measure( TNL::Logger& logger )
{
   simulationMonitor.template measure< SPHKernelFunction, EOS >( fluid, boundary, modelParams, timeStepping, logger, verbose );
}

template< typename Model >
void
SPHMultiset_CFD< Model >::save( TNL::Logger& logger, bool writeParticleCellIndex )
{
   if( verbose == "with-snapshot" )
      writeInfo( logger );

   const int step = timeStepping.getStep();

   std::string outputFileNameFluid = outputDirecotry + "/particles" + std::to_string( step ) + "_fluid.vtk";
   fluid->template writeParticlesAndVariables< Writer >( outputFileNameFluid, writeParticleCellIndex );
   logger.writeParameter( "Saved:", outputFileNameFluid );

   std::string outputFileNameBound = outputDirecotry + "/particles" + std::to_string( step ) + "_boundary.vtk";
   boundary->template writeParticlesAndVariables< Writer >( outputFileNameBound, writeParticleCellIndex );
   logger.writeParameter( "Saved:", outputFileNameBound );

   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ){
      for( auto& openBoundaryPatch : openBoundaryPatches ){
         std::string outputFileNameOpenBound = outputDirecotry + "/particles" + std::to_string( step ) +\
                                               "_" + openBoundaryPatch->parameters.identifier + ".vtk";
         openBoundaryPatch->template writeParticlesAndVariables< Writer >( outputFileNameOpenBound, writeParticleCellIndex );
         logger.writeParameter( "Saved:", outputFileNameOpenBound );
      }
   }
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfPeriodicBuffers > 0 ){
      for( auto& periodicBoundaryPatch : periodicBoundaryPatches ){
         std::string outputFileNamePeriodicBound = outputDirecotry + "/particles" + std::to_string( step ) +\
                                                   "_" + periodicBoundaryPatch->fluidPeriodicPatch->parameters.identifier + ".vtk";
         periodicBoundaryPatch->fluidPeriodicPatch->template writeParticlesAndVariables< Writer >( outputFileNamePeriodicBound, writeParticleCellIndex );
         logger.writeParameter( "Saved:", outputFileNamePeriodicBound );

         std::string outputFileNamePeriodicBound_bound = outputDirecotry + "/particles" + std::to_string( step ) +\
                                                         "_" + periodicBoundaryPatch->boundaryPeriodicPatch->parameters.identifier + "_boundary.vtk";
         periodicBoundaryPatch->boundaryPeriodicPatch->template writeParticlesAndVariables< Writer >( outputFileNamePeriodicBound_bound, writeParticleCellIndex );
         logger.writeParameter( "Saved:", outputFileNamePeriodicBound_bound );
      }
   }

   // output simulation sensors to files
   simulationMonitor.save( logger );
}

template< typename Model >
void
SPHMultiset_CFD< Model >::writeProlog( TNL::Logger& logger, bool writeSystemInformation ) const noexcept
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
   logger.writeParameter( "Output directory:", outputDirecotry );
   logger.writeParameter( "Particles format", particlesFormat );
   writePrologModel( logger, modelParams );
   logger.writeHeader( "Fluid object information." );
   fluid->writeProlog( logger );
   logger.writeHeader( "Boundary object information:" );
   boundary->writeProlog( logger );
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ){
      for( long unsigned int i = 0; i < openBoundaryPatches.size(); i++ )
      {
         logger.writeHeader( "Open boundary buffer" + std::to_string( i + 1 ) + "." );
         openBoundaryPatches[ i ]->writeProlog( logger );
         logger.writeSeparator();
         openBoundaryPatches[ i ]->config.writeProlog( logger );
      }
   }
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfPeriodicBuffers > 0 ){
      for( long unsigned int i = 0; i < periodicBoundaryPatches.size(); i++ )
      {
         logger.writeHeader( "Open boundary buffer" + std::to_string( i + 1 ) + "." );
         periodicBoundaryPatches[ i ]->fluidPeriodicPatch->writeProlog( logger );
         logger.writeSeparator();
         periodicBoundaryPatches[ i ]->config.writeProlog( logger );
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
SPHMultiset_CFD< Model >::writeLog( TNL::Logger& logger,
                                    const std::string& label,
                                    const ParameterType& value,
                                    int parameterLevel )
{
   if( verbose == "full" )
      logger.writeParameter( label, value, parameterLevel );
}

template< typename Model >
void
SPHMultiset_CFD< Model >::writeInfo( TNL::Logger& logger ) const noexcept
{
   logger.writeSeparator();
   logger.writeParameter( "Simulation time: " + std::to_string( timeStepping.getTime() ) + " s, simulation step: " \
                          + std::to_string( timeStepping.getStep() ), "" );
   logger.writeCurrentTime( "Current time:" );
   logger.writeParameter( "Number of fluid particles:", fluid->getNumberOfParticles() );
   if( verbose == "full" ){
      logger.writeParameter( "Fluid object - first active particle:", fluid->getFirstActiveParticle(), 1 );
      logger.writeParameter( "Fluid particles - first active particle:", fluid->particles->getFirstActiveParticle(), 1 );
      logger.writeParameter( "Fluid object - last active particle:", fluid->getLastActiveParticle(), 1 );
      logger.writeParameter( "Fluid particles - last active particle:", fluid->particles->getLastActiveParticle(), 1 );
   }
   logger.writeParameter( "Number of boundary particles:", boundary->getNumberOfParticles() );
   if( verbose == "full" ){
      logger.writeParameter( "Boundary object - first active particle:", boundary->getFirstActiveParticle(), 1 );
      logger.writeParameter( "Boundary particles - first active particle:", boundary->particles->getFirstActiveParticle(), 1 );
      logger.writeParameter( "Boundary object - last active particle:", boundary->getLastActiveParticle(), 1 );
      logger.writeParameter( "Boundary particles - last active particle:", boundary->particles->getLastActiveParticle(), 1 );
   }
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfBoundaryBuffers > 0 ){
      for( long unsigned int i = 0; i < openBoundaryPatches.size(); i++ )
      {
         logger.writeParameter( "Number of buffer" + std::to_string( i + 1 ) + " particles:",
                                openBoundaryPatches[ i ]->getNumberOfParticles() );
         if( verbose == "full" ){
            logger.writeParameter( "Patch " + openBoundaryPatches[ i ]->config.identifier + " object - first active particle:", openBoundaryPatches[ i ]->getFirstActiveParticle(), 1 );
            logger.writeParameter( "Patch " + openBoundaryPatches[ i ]->config.identifier + " particles - first active particle:", openBoundaryPatches[ i ]->particles->getFirstActiveParticle(), 1 );
            logger.writeParameter( "Patch " + openBoundaryPatches[ i ]->config.identifier + " object - last active particle:", openBoundaryPatches[ i ]->getLastActiveParticle(), 1 );
            logger.writeParameter( "Patch " + openBoundaryPatches[ i ]->config.identifier + " particles - last active particle:", openBoundaryPatches[ i ]->particles->getLastActiveParticle(), 1 );
         }
      }
   }
   if constexpr( Model::ModelConfigType::SPHConfig::numberOfPeriodicBuffers > 0 ){
      for( long unsigned int i = 0; i < periodicBoundaryPatches.size(); i++ )
      {
         logger.writeParameter( "Number of buffer" + std::to_string( i + 1 ) + " fluid particles:",
                                periodicBoundaryPatches[ i ]->fluidPeriodicPatch->getNumberOfParticles() );
         logger.writeParameter( "Number of buffer" + std::to_string( i + 1 ) + " boundary particles:",
                                periodicBoundaryPatches[ i ]->boundaryPeriodicPatch->getNumberOfParticles() );
         if( verbose == "full" ){
         }
      }
   }
   logger.writeSeparator();
}

template< typename Model >
void
SPHMultiset_CFD< Model >::writeEpilog( TNL::Logger& logger ) const noexcept
{
   logger.writeHeader( "SPH simulation successfully finished." );
   logger.writeCurrentTime( "Ended at:" );
   timeMeasurement.writeInfo( logger, timeStepping.getStep() );
}

} // SPH
} // TNL
