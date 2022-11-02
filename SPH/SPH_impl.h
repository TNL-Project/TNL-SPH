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
   model.sortParticlesAndVariables(); //I DONT NEED TO DO THIS IN EACH STEP!
   //particles.sortParticles();

   /**
    * Find neigbors.
    */
   neighborSearch.searchForNeighbors();
}

template< typename Variables, typename ParticleSystem, typename NeighborSearch >
void
SPHSimulation< Variables, ParticleSystem, NeighborSearch >::ProcessOneParticle(GlobalIndexType index_i)
{
  const LocalIndexType numberOfNeigbors = particles.getNeighborsCount( index_i );
  printf(" Particle i: %d has %d of nbs.\n", index_i, numberOfNeigbors);

  auto fetch = [=] __cuda_callable__ ( int i ) -> InteractionResultType
  {

    //This has to be moved smehod to model site. (The function processOneParticle should be in model class.)
    GlobalIndexType index_j = particles.getNeighbor( index_i, i );

    if( model.vars.type[ index_j ] == 0 )
      return model.template PerformParticleInteractionFF< WendlandKernel >( index_i , index_j);
    else
      return model.template PerformParticleInteractionFB< WendlandKernel >( index_i , index_j);

  };

  auto reduction = [] __cuda_callable__ ( const InteractionResultType& a, const InteractionResultType& b )
  {
     return a + b;
  };

  model.vars.DrhoDv[ index_i ] = Algorithms::reduce< DeviceType, int, InteractionResultType >( 0, numberOfNeigbors, fetch, reduction, {0., 0., -9.81} );
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
