#include "SPH.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Variables, typename ParticleSystem, typename NeighborSearch >
void
SPHSimpleFluid< Variables, ParticleSystem, NeighborSearch >::PerformNeighborSearch( GlobalIndexType step, TNL::Timer& timer_reset, TNL::Timer& timer_cellIndices, TNL::Timer& timer_sort, TNL::Timer& timer_toCells )
{
   /**
    * Compute gird nad partice cell indices.
    */
   timer_reset.start();
   neighborSearch->resetListWithIndices();
   if( step == 0 )
   neighborSearch_bound->resetListWithIndices();
   timer_reset.stop();
   std::cout << " - neighborSearch->resetListWithIndices();... done" << std::endl;

   if( step == 0 ) //TODO: do this better
   particles->computeGridCellIndices();

   timer_cellIndices.start();
   particles->computeParticleCellIndices();
   if( step == 0 )
   particles_bound->computeParticleCellIndices();
   timer_cellIndices.stop();
   std::cout << " - particles->computeParticleCellIndices();... done " << std::endl;

   timer_sort.start();
   model->sortParticlesAndVariablesThrust( model->particles, model->FluidVariables, model->swapFluid );
   integrator->sortIntegratorArrays();
   if( step == 0 )
   {
      model->sortParticlesAndVariablesThrust( model->boundaryParticles, model->BoundaryVariables, model->swapBoundary );
      integrator->sortIntegratorBoundaryArrays();
   }
   timer_sort.stop();
   std::cout << " - model->sortParticlesAndVariables();... done " << std::endl;

   timer_toCells.start();
   neighborSearch->particlesToCells();
   if( step == 0 )
   neighborSearch_bound->particlesToCells();
   timer_toCells.stop();
   std::cout << " - neighborSearch->particlesToCells();... done " << std::endl;
}

template< typename Variables, typename ParticleSystem, typename NeighborSearch >
template< typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm, typename EOS >
void
SPHSimpleFluid< Variables, ParticleSystem, NeighborSearch >::InteractModel()
{
   model->template Interaction< NeighborSearchPointer, SPHKernelFunction, DiffusiveTerm, ViscousTerm, EOS >( neighborSearch, neighborSearch_bound );
}


} // SPH
} // ParticleSystem
} // TNL

