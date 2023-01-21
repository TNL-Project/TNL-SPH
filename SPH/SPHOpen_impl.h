#include "SPH.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Variables, typename ParticleSystem, typename NeighborSearch >
void
SPHOpenSystem< Variables, ParticleSystem, NeighborSearch >::PerformNeighborSearch( GlobalIndexType step, TNL::Timer& timer_reset, TNL::Timer& timer_cellIndices, TNL::Timer& timer_sort, TNL::Timer& timer_toCells )
{
   /**
    * Compute gird nad partice cell indices.
    */
   timer_reset.start();
   fluid->neighborSearch->resetListWithIndices();
   if( step == 0 )
      boundary->neighborSearch->resetListWithIndices();
   openBoundaryPatches[ 0 ]->neighborSearch->resetListWithIndices();
   timer_reset.stop();
   std::cout << " - neighborSearch->resetListWithIndices();... done" << std::endl;

   if( step == 0 ) //TODO: do this better, i. e. move this to constructor
   fluid->particles->computeGridCellIndices(); //with current settings, I dont need to do this

   timer_cellIndices.start();
   fluid->particles->computeParticleCellIndices();
   openBoundaryPatches[ 0 ]->particles->computeParticleCellIndices();
   if( step == 0 )
      boundary->particles->computeParticleCellIndices();
   timer_cellIndices.stop();
   std::cout << " - particles->computeParticleCellIndices();... done " << std::endl;

   /**
    * Sort particles.
    */
   timer_sort.start();
   fluid->sortParticles();
   openBoundaryPatches[ 0 ]->sortParticles();
   if( step == 0 )
   {
      boundary->sortParticles();
   }
   timer_sort.stop();
   std::cout << " - model->sortParticlesAndVariables();... done " << std::endl;

   /**
    * Bucketing, particles to cells.
    */
   timer_toCells.start();
   fluid->neighborSearch->particlesToCells();
   openBoundaryPatches[ 0 ]->neighborSearch->particlesToCells();
   if( step == 0 )
      boundary->neighborSearch->particlesToCells();
   timer_toCells.stop();
   std::cout << " - neighborSearch->particlesToCells();... done " << std::endl;
}

template< typename Variables, typename ParticleSystem, typename NeighborSearch >
template< typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm, typename EOS >
void
SPHOpenSystem< Variables, ParticleSystem, NeighborSearch >::Interact()
{
   model->template Interaction<
      FluidPointer, BoundaryPointer, OpenBoudaryPatchPointer, NeighborSearchPointer,
      SPHKernelFunction, DiffusiveTerm, ViscousTerm, EOS >( fluid, boundary, openBoundaryPatches[ 0 ] );
}

template< typename Variables, typename ParticleSystem, typename NeighborSearch >
void
SPHOpenSystem< Variables, ParticleSystem, NeighborSearch >::addOpenBoundaryPatch( GlobalIndexType size_buffer, GlobalIndexType sizeAllocated_buffer, RealType h, GlobalIndexType numberOfCells )
{
   openBoundaryPatches.emplace_back( size_buffer, sizeAllocated_buffer, h, numberOfCells );
}

} // SPH
} // ParticleSystem
} // TNL

