#include "SPH.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Variables, typename ParticleSystem, typename NeighborSearch >
void
SPHSimulation< Variables, ParticleSystem, NeighborSearch >::PerformNeighborSearch( GlobalIndexType step )
{
   /**
    * Compute gird nad partice cell indices.
    */
   particles.resetNeighborList();
   particles.computeGridCellIndices(); //I DONT NEED TO REPEAT THIS!
   std::cout << "SPHSimulation::PerformNeighborSearch(): ... OK" << std::endl; //debug
   particles.computeParticleCellIndices();
   std::cout << "SPHSimulation::computeParticleCellIndices(): ... OK" << std::endl; //debug

	 if( step % 100 == 0 )
	 {
    	model.sortParticlesAndVariables(); //I DONT NEED TO DO THIS IN EACH STEP!
    	std::cout << "SPHSimulation::sortParticlesAndVariables(): ... OK" << std::endl; //debug
    	//particles.sortParticles();
	 }

   //std::cout << particles.getParticleCellIndices() << std::endl; //debug

   /**
    * Find neigbors.
    */
   neighborSearch.searchForNeighbors();
   std::cout << "neighborSearch.searchForNeighbors() ... OK" << std::endl; //debug
}

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

