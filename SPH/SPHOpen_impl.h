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
   openBoundaryPatch->neighborSearch->resetListWithIndices();
   timer_reset.stop();
   std::cout << " - neighborSearch->resetListWithIndices();... done" << std::endl;

   if( step == 0 ) //TODO: do this better, i. e. move this to constructor
   fluid->particles->computeGridCellIndices(); //with current settings, I dont need to do this

   timer_cellIndices.start();
   fluid->particles->computeParticleCellIndices();
   openBoundaryPatch->particles->computeParticleCellIndices();
   if( step == 0 )
      boundary->particles->computeParticleCellIndices();
   timer_cellIndices.stop();
   std::cout << " - particles->computeParticleCellIndices();... done " << std::endl;

   /**
    * Sort particles.
    */
   timer_sort.start();
   model->sortParticlesAndVariablesThrust( fluid->particles, fluid->variables, model->swapFluid );
   //fluid->sortParticles();
   //model->sortVariables( fluid->particles, fluid->variables, fluid->sortPermutations );

   model->sortParticlesAndVariablesThrust( openBoundaryPatch->particles, openBoundaryPatch->variables, model->swapInlet );
   integrator->sortIntegratorArrays( fluid->particles->getNumberOfParticles() ); //TODO: NumberOfPtcs.
   if( step == 0 )
   {
      //model->sortParticlesAndVariablesThrust( model->boundaryParticles, model->BoundaryVariables, model->swapBoundary );
      model->sortBoundaryParticlesAndVariablesThrust( boundary->particles, boundary->variables, model->swapBoundary );
      integrator->sortIntegratorBoundaryArrays( boundary->particles->getNumberOfParticles() ); //TODO: NumberOfPtcs.
   }
   timer_sort.stop();
   std::cout << " - model->sortParticlesAndVariables();... done " << std::endl;

   /**
    * Bucketing, particles to cells.
    */
   timer_toCells.start();
   fluid->neighborSearch->particlesToCells();
   openBoundaryPatch->neighborSearch->particlesToCells();
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
      SPHKernelFunction, DiffusiveTerm, ViscousTerm, EOS >( fluid, boundary, openBoundaryPatch );
}


} // SPH
} // ParticleSystem
} // TNL

