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

//{ dep }: template< typename Variables, typename ParticleSystem, typename NeighborSearch >
//{ dep }: void
//{ dep }: SPHSimulation< Variables, ParticleSystem, NeighborSearch >::ProcessOneParticle(GlobalIndexType index_i)
//{ dep }: {
//{ dep }:   const LocalIndexType numberOfNeigbors = particles.getNeighborsCount( index_i );
//{ dep }:   printf(" Particle i: %d has %d of nbs.\n", index_i, numberOfNeigbors);
//{ dep }:
//{ dep }:   auto fetch = [=] __cuda_callable__ ( int i ) -> InteractionResultType
//{ dep }:   {
//{ dep }:
//{ dep }:     GlobalIndexType index_j = particles.getNeighbor( index_i, i );
//{ dep }:
//{ dep }:     if( model.vars.type[ index_j ] == 0 )
//{ dep }:       return model.template PerformParticleInteractionFF< WendlandKernel >( index_i , index_j);
//{ dep }:     else
//{ dep }:       return model.template PerformParticleInteractionFB< WendlandKernel >( index_i , index_j);
//{ dep }:
//{ dep }:   };
//{ dep }:
//{ dep }:   auto reduction = [] __cuda_callable__ ( const InteractionResultType& a, const InteractionResultType& b )
//{ dep }:   {
//{ dep }:      return a + b;
//{ dep }:   };
//{ dep }:
//{ dep }:   model.vars.DrhoDv[ index_i ] = Algorithms::reduce< DeviceType, int, InteractionResultType >( 0, numberOfNeigbors, fetch, reduction, {0., 0., -9.81} );
//{ deprecated }: }

template< typename Variables, typename ParticleSystem, typename NeighborSearch >
void
SPHSimulation< Variables, ParticleSystem, NeighborSearch >::Interact()
{
    auto init = [=] __cuda_callable__ ( int i ) mutable
    {
       model.ProcessOneParticle( i );
    };
    Algorithms::ParallelFor< DeviceType >::exec( 0, particles.getNumberOfParticles(), init );
}

} // SPH
} // ParticleSystem
} // TNL
