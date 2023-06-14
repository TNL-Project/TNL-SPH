#include "SPH.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Variables, typename ParticleSystem >
void
SPHSimpleFluid< Variables, ParticleSystem >::PerformNeighborSearch( GlobalIndexType step, TNL::Timer& timer_reset, TNL::Timer& timer_cellIndices, TNL::Timer& timer_sort, TNL::Timer& timer_toCells )
{
   /**
    * Compute gird nad partice cell indices.
    */
   timer_reset.start();
   fluid->particles->resetListWithIndices();
   if( step == 0 )
      boundary->particles->resetListWithIndices();
   timer_reset.stop();
   std::cout << " - neighborSearch->resetListWithIndices();... done" << std::endl;

   //deprecated: if( step == 0 ) //TODO: do this better, i. e. move this to constructor
   //deprecated: fluid->particles->computeGridCellIndices(); //with current settings, I dont need to do this

   timer_cellIndices.start();
   fluid->particles->computeParticleCellIndices();
   if( step == 0 )
      boundary->particles->computeParticleCellIndices();
   timer_cellIndices.stop();
   std::cout << " - particles->computeParticleCellIndices();... done " << std::endl;

   /**
    * Sort particles.
    */
   timer_sort.start();
   fluid->sortParticles();
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
   fluid->particles->particlesToCells();
   if( step == 0 )
      boundary->particles->particlesToCells();
   timer_toCells.stop();
   std::cout << " - neighborSearch->particlesToCells();... done " << std::endl;
}

template< typename Variables, typename ParticleSystem >
template< typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm, typename EOS, typename SPHState >
void
SPHSimpleFluid< Variables, ParticleSystem >::Interact( SPHState& sphState )
{
   model->template Interaction< FluidPointer, BoundaryPointer, SPHKernelFunction, DiffusiveTerm, ViscousTerm, EOS >(
         fluid, boundary, sphState );
}

template< typename Variables, typename ParticleSystem >
template< typename SPHKernelFunction, typename RiemannSolver, typename EOS, typename SPHState >
void
SPHSimpleFluid< Variables, ParticleSystem >::Interact( SPHState& sphState )
{
   model->template Interaction< FluidPointer, BoundaryPointer, SPHKernelFunction, RiemannSolver, EOS >(
         fluid, boundary, sphState );
}

template< typename Variables, typename ParticleSystem >
template< typename Writer >
void
SPHSimpleFluid< Variables, ParticleSystem >::save( const std::string& outputFileName, const int step )
{
   std::string outputFileNameFluid = outputFileName + std::to_string( step ) + "_fluid.vtk";
   fluid->template writeParticlesAndVariables< Writer >( outputFileNameFluid );

   std::string outputFileNameBound = outputFileName + std::to_string( step ) + "_boundary.vtk";
   boundary->template writeParticlesAndVariables< Writer >( outputFileNameBound );
}

template< typename Variables, typename ParticleSystem >
void
SPHSimpleFluid< Variables, ParticleSystem >::writeProlog( TNL::Logger& logger ) const noexcept
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

