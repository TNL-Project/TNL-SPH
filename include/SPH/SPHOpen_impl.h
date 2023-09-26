#include "SPH.h"
#include "SPHOpen.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Model >
void
SPHOpenSystem< Model >::PerformNeighborSearch( GlobalIndexType step, TNL::Timer& timer_reset, TNL::Timer& timer_cellIndices, TNL::Timer& timer_sort, TNL::Timer& timer_toCells )
{
   /**
    * Compute gird nad partice cell indices.
    */
   timer_reset.start();
   fluid->particles->resetListWithIndices();
   //openBoundaryPatches[ 0 ]->neighborSearch->resetListWithIndices();
   for( auto& openBoundaryPatch : openBoundaryPatches )
      openBoundaryPatch->particles->resetListWithIndices();

   if( step == 0 )
      boundary->particles->resetListWithIndices();
   timer_reset.stop();
   std::cout << " - neighborSearch->resetListWithIndices();... done" << std::endl;

   //if( step == 0 ) //TODO: do this better, i. e. move this to constructor
   //fluid->particles->computeGridCellIndices(); //with current settings, I dont need to do this

   timer_cellIndices.start();
   fluid->particles->computeParticleCellIndices();
   //openBoundaryPatches[ 0 ]->particles->computeParticleCellIndices();
   for( auto& openBoundaryPatch : openBoundaryPatches )
      openBoundaryPatch->particles->computeParticleCellIndices();

   if( step == 0 )
      boundary->particles->computeParticleCellIndices();
   timer_cellIndices.stop();
   std::cout << " - particles->computeParticleCellIndices();... done " << std::endl;

   /**
    * Sort particles.
    */
   timer_sort.start();
   fluid->sortParticles();
   //openBoundaryPatches[ 0 ]->sortParticles();
   for( auto& openBoundaryPatch : openBoundaryPatches )
      openBoundaryPatch->sortParticles();

   if( step == 0 )
   {
      boundary->sortParticles();
   }
   timer_sort.stop();
   std::cout << " - model->sortParticlesAndVariables();... done " << std::endl;

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
   timer_toCells.start();
   fluid->particles->particlesToCells();
   //openBoundaryPatches[ 0 ]->neighborSearch->particlesToCells();
   for( auto& openBoundaryPatch : openBoundaryPatches )
      openBoundaryPatch->particles->particlesToCells();
   timer_toCells.stop();

   if( step == 0 )
      boundary->particles->particlesToCells();
   timer_toCells.stop();
   std::cout << " - neighborSearch->particlesToCells();... done " << std::endl;
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

   //model->template interactionWithOpenBoundary< FluidPointer, OpenBoundaryPointer, SPHKernelFunction, DiffusiveTerm, ViscousTerm, EOS >(
   //      fluid, openBoundaryPatches[ 0 ], sphState );

   model->template interactionWithOpenBoundary< FluidPointer, BoundaryPointer, OpenBoundaryPointer, SPHKernelFunction, DiffusiveTerm, ViscousTerm, EOS >(
         fluid, boundary, openBoundaryPatches[ 0 ], sphState );
   model->template interactionWithOpenBoundary< FluidPointer, BoundaryPointer, OpenBoundaryPointer, SPHKernelFunction, DiffusiveTerm, ViscousTerm, EOS >(
         fluid, boundary, openBoundaryPatches[ 1 ], sphState );


   //Interact buffers
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
      //openBoundaryPatches[ 0 ]->particles->setGridSize( sphConfig.gridSize );
      //openBoundaryPatches[ 0 ]->particles->setGridOrigin( sphConfig.gridOrigin );

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
   logger.writeParameter( "Number of fluid particles:", this->fluid->particles->getNumberOfParticles() );
   logger.writeParameter( "Number of alloc. fluid particles:", this->fluid->particles->getNumberOfAllocatedParticles() );

   logger.writeParameter( "Number of boundary particles:", this->boundary->particles->getNumberOfParticles() );
   logger.writeParameter( "Number of alloc. boundary particles:", this->boundary->particles->getNumberOfParticles() );

   logger.writeParameter( "Particle grid size: ", this->fluid->particles->getGridSize() );
   logger.writeParameter( "Particle grid origin: ", this->fluid->particles->getGridOrigin() );

   logger.writeParameter( "Particle boundary grid size: ", this->boundary->particles->getGridSize() );
   logger.writeParameter( "Particle boundary grid origin: ", this->boundary->particles->getGridOrigin() );
}

} // SPH
} // ParticleSystem
} // TNL

