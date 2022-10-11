#pragma once

#include <TNL/Containers/Vector.h>
#include <TNL/Algorithms/reduce.h>

#include "../Particles/Particles.h"

#include "SPHFluidTraits.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template < typename Variables, typename ParticleSystem, typename NeighborSearch >
class SPHSimulation
{
public:

  using LocalIndexType = typename ParticleSystem::LocalIndexType;
  using GlobalIndexType = typename ParticleSystem::GlobalIndexType;
  using DeviceType = typename ParticleSystem::Device;

  SPHSimulation() = default;

  SPHSimulation(GlobalIndexType size, float h)
  : particles(size, h), neighborSearch(particles, 100), vars(size) {};


//protected:
  void PerformNeighborSearch()
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

  void ProcessOneParticle(GlobalIndexType index_i)
  {

    const LocalIndexType numberOfNeigbors = particles.getNeighborsCount( index_i );
    printf(" Particle i: %d has %d of nbs.", index_i, numberOfNeigbors);

    auto fetch = [=] __cuda_callable__ ( int i ) -> double { return particles.getNeighbor( index_i, i ); };
    auto reduction = [] __cuda_callable__ ( const double& a, const double& b ) { return a + b; };

    vars.p[index_i] = Algorithms::reduce< DeviceType >( 0, numberOfNeigbors, fetch, reduction, 0.0 );

  }

  //INTERACET
  // - cycle over the neighbors
  void Interact()
  {

    // - for ALL PTCS
    //    - for ALL NBCS of PTC

    auto init = [=] __cuda_callable__ ( int i ) mutable
    {
       ProcessOneParticle( i );
    };
    Algorithms::ParallelFor< DeviceType >::exec( 0, particles.getNumberOfParticles(), init );

  }

  //Particles

  ParticleSystem particles;
  NeighborSearch neighborSearch;

  Variables vars;

private:

};

} // SPH
} // ParticleSystem
} // TNL
