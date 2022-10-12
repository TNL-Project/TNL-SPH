#include "SPH.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Variables, typename ParticleSystem, typename NeighborSearch >
void
SPHSimulation< Variables, ParticleSystem, NeighborSearch >::PerformNeighborSearch()
{
   /**
    * Compute gird nad partice cell indices.
    */
   particles.computeGridCellIndices(); //I DONT NEED TO REPEAT THIS!
   particles.computeParticleCellIndices();

   /**
    * Assign particles to neighborSearch arrays.
    */
   neighborSearch.particlesToCells();

   /**
    * TEMP
    * Find neigbors.
    */
   neighborSearch.runCycleOverGrid();
}

template< typename Variables, typename ParticleSystem, typename NeighborSearch >
void
SPHSimulation< Variables, ParticleSystem, NeighborSearch >::ProcessOneParticle(GlobalIndexType index_i)
{
  const LocalIndexType numberOfNeigbors = particles.getNeighborsCount( index_i );
  printf(" Particle i: %d has %d of nbs.\n", index_i, numberOfNeigbors);

  /* test */
  auto fetch = [=] __cuda_callable__ ( int i ) -> double { return particles.getNeighbor( index_i, i ); };
  auto reduction = [] __cuda_callable__ ( const double& a, const double& b ) { return a + b; };

  /* real interaction */
  auto fetch_interatcion = [=] __cuda_callable__ ( int i ) -> double
  {
     return particles.getNeighbor( index_i, i );
  };

  auto reduction_interaction = [] __cuda_callable__ ( const double& a, const double& b )
  {
     return a + b;
  };

  model.vars.p[index_i] = Algorithms::reduce< DeviceType >( 0, numberOfNeigbors, fetch, reduction, 0.0 );
}

template< typename Variables, typename ParticleSystem, typename NeighborSearch >
void
SPHSimulation< Variables, ParticleSystem, NeighborSearch >::Interact()
{
    auto init = [=] __cuda_callable__ ( int i ) mutable
    {
       ProcessOneParticle( i );
    };
    Algorithms::ParallelFor< DeviceType >::exec( 0, particles.getNumberOfParticles(), init );
}

} // SPH
} // ParticleSystem
} // TNL
