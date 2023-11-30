#include "SPH.h"
#include "SPH/OpenBoundaryConfig.h"
#include "SPH/TimeMeasurement.h"
#include "SPHOpen.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Model >
void
SPHOpenSystem< Model >::performNeighborSearch( GlobalIndexType step, TimerMeasurement& timeMeasurement, TNL::Logger& log )
{
   /**
    * Compute gird nad partice cell indices.
    */
   timeMeasurement.start( "search_reset" );
   fluid->particles->resetListWithIndices();
   for( auto& openBoundaryPatch : openBoundaryPatches )
      openBoundaryPatch->particles->resetListWithIndices();

   if( step == 0 )
      boundary->particles->resetListWithIndices();
   timeMeasurement.stop( "search_reset" );
   log.writeParameter( "Search - reset ...", "Done." );

   /**
    * Compute gird and partice cell indices.
    */
   timeMeasurement.start( "search_cellIndices" );
   fluid->particles->computeParticleCellIndices();
   for( auto& openBoundaryPatch : openBoundaryPatches )
      openBoundaryPatch->particles->computeParticleCellIndices();

   if( step == 0 )
      boundary->particles->computeParticleCellIndices();
   timeMeasurement.stop( "search_cellIndices" );
   log.writeParameter( "Search - compute cell indices ...", "Done." );

   /**
    * Sort particles.
    */
   timeMeasurement.start( "search_sort" );
   fluid->sortParticles();
   for( auto& openBoundaryPatch : openBoundaryPatches )
      openBoundaryPatch->sortParticles();

   if( step == 0 )
   {
      boundary->sortParticles();
   }
   timeMeasurement.stop( "search_sort" );
   log.writeParameter( "Search - sort ...", "Done." );

   /**
    * Update number of fluid particles dute to the removed once.
    */
   fluid->particles->setNumberOfParticles( fluid->particles->getNumberOfParticles() - openBoundaryPatches[ 1 ]->numberOfFluidParticlesToRemove );
   fluid->particles->setLastActiveParticle( fluid->particles->getLastActiveParticle() - openBoundaryPatches[ 1 ]->numberOfFluidParticlesToRemove );
   fluid->setLastActiveParticle( fluid->getLastActiveParticle() - openBoundaryPatches[ 1 ]->numberOfFluidParticlesToRemove );
   openBoundaryPatches[ 1 ]->numberOfFluidParticlesToRemove = 0;

   /**
    * Bucketing, particles to cells.
    */
   timeMeasurement.start( "search_toCells" );
   fluid->particles->particlesToCells();
   for( auto& openBoundaryPatch : openBoundaryPatches )
      openBoundaryPatch->particles->particlesToCells();

   if( step == 0 )
      boundary->particles->particlesToCells();
   timeMeasurement.stop( "search_toCells" );
   log.writeParameter( "Search - particles to cells ...", "Done." );
}

template< typename Model >
template< typename SPHKernelFunction, typename EOS, typename SPHState >
void
SPHOpenSystem< Model >::extrapolateOpenBC( SPHState& sphState, std::vector< OpenBoundaryConfigType >& bufferParams )
{
   for( int i = 0; i < std::size( bufferParams ); i++ )
      model->template extrapolateOpenBoundaryData< FluidPointer, OpenBoundaryPointer, SPHKernelFunction, EOS >(
            fluid, openBoundaryPatches[ i ], sphState, bufferParams[ i ] );
}

template< typename Model >
template< typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm, typename EOS, typename SPHState >
void
SPHOpenSystem< Model >::interact( SPHState& sphState )
{

   model->template interaction< FluidPointer, BoundaryPointer, SPHKernelFunction, DiffusiveTerm, ViscousTerm, EOS >(
         fluid, boundary, sphState );
   model->template updateSolidBoundary< FluidPointer, BoundaryPointer, SPHKernelFunction, DiffusiveTerm, ViscousTerm, EOS >(
         fluid, boundary, sphState );

   //Interact buffers
   model->template interactionWithOpenBoundary<
      FluidPointer, BoundaryPointer, OpenBoundaryPointer, SPHKernelFunction, DiffusiveTerm, ViscousTerm, EOS >(
         fluid, boundary, openBoundaryPatches[ 0 ], sphState );
   model->template interactionWithOpenBoundary<
      FluidPointer, BoundaryPointer, OpenBoundaryPointer, SPHKernelFunction, DiffusiveTerm, ViscousTerm, EOS >(
         fluid, boundary, openBoundaryPatches[ 1 ], sphState );

   model->finalizeInteraction( fluid, boundary, sphState );

}

template< typename Model >
template< typename SPHOpenSystemInit >
void
SPHOpenSystem< Model >::addOpenBoundaryPatch( SPHOpenSystemInit sphConfig )
{
   for( int i = 0; i < std::size( sphConfig.numberOfOpenBoundaryParticles ); i++ ){
      openBoundaryPatches.emplace_back( sphConfig.numberOfOpenBoundaryParticles[ i ],
                                        sphConfig.numberOfAllocatedOpenBoundaryParticles[ i ],
                                        sphConfig.searchRadius,
                                        sphConfig.numberOfGridCells );

      openBoundaryPatches[ i ]->particles->setGridSize( sphConfig.gridSize );
      openBoundaryPatches[ i ]->particles->setGridOrigin( sphConfig.gridOrigin );
   }

}

template< typename Model >
template< typename SPHOpenSystemInit >
void
SPHOpenSystem< Model >::addOpenBoundaryPatch( SPHOpenSystemInit sphConfig, std::vector< OpenBoundaryConfigType >& bufferParams )
{

   //using TestBC = OpenBoundaryConfig< TNL::ParticleSystem::SPH::One, SPHConfig >;
   //openBoundaryPatchesConfigs.push_back( std::make_shared< TestBC >() );

   for( int i = 0; i < std::size( sphConfig.numberOfOpenBoundaryParticles ); i++ ){
      openBoundaryPatches.emplace_back( sphConfig.numberOfOpenBoundaryParticles[ i ],
                                        sphConfig.numberOfAllocatedOpenBoundaryParticles[ i ],
                                        sphConfig.searchRadius,
                                        sphConfig.numberOfGridCells );

      openBoundaryPatches[ i ]->particles->setGridSize( sphConfig.gridSize );
      openBoundaryPatches[ i ]->particles->setGridOrigin( sphConfig.gridOrigin );

   }

}

template< typename Model >
template< typename Writer >
void
SPHOpenSystem< Model >::save( const std::string& outputFileName, const int step, bool writeParticleCellIndex )
{
   std::string outputFileNameFluid = outputFileName + std::to_string( step ) + "_fluid.vtk";
   fluid->template writeParticlesAndVariables< Writer >( outputFileNameFluid );

   std::string outputFileNameBound = outputFileName + std::to_string( step ) + "_boundary.vtk";
   boundary->template writeParticlesAndVariables< Writer >( outputFileNameBound );

   for( auto& openBoundaryPatch : openBoundaryPatches ){
      std::string outputFileNameOpenBound = outputFileName + std::to_string( step ) +\
                                            "_" + openBoundaryPatch->parameters.identifier + ".vtk";
      openBoundaryPatch->template writeParticlesAndVariables< Writer >( outputFileNameOpenBound );
   }
}

template< typename Model >
void
SPHOpenSystem< Model >::writeProlog( TNL::Logger& logger ) const noexcept
{

   logger.writeHeader( "SPH Open: Initial simulation configuration." );
   logger.writeParameter( "Number of fluid particles:", this->fluid->particles->getNumberOfParticles() );
   logger.writeParameter( "Number of alloc. fluid particles:", this->fluid->particles->getNumberOfAllocatedParticles() );

   logger.writeParameter( "Number of boundary particles:", this->boundary->particles->getNumberOfParticles() );
   logger.writeParameter( "Number of alloc. boundary particles:", this->boundary->particles->getNumberOfParticles() );

   logger.writeParameter( "Particle grid size: ", this->fluid->particles->getGridSize() );
   logger.writeParameter( "Particle grid origin: ", this->fluid->particles->getGridOrigin() );

   logger.writeParameter( "Particle boundary grid size: ", this->boundary->particles->getGridSize() );
   logger.writeParameter( "Particle boundary grid origin: ", this->boundary->particles->getGridOrigin() );
   logger.writeSeparator();
}

} // SPH
} // ParticleSystem
} // TNL

