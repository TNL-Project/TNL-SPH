#pragma once

#include <TNL/Containers/Vector.h>
#include <TNL/Algorithms/reduce.h>

#include "../Particles/Particles.h"

#include "SPHFluidTraits.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Variables, typename ParticleSystem, typename NeighborSearch >
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

  /**
   * Perform neighbors search and fill neighborsList in Particle system variable.
   */
  void PerformNeighborSearch();

  /**
   * Proces one particle (i.e. loop over all its neighbors and perform interactions).
   */
  void ProcessOneParticle(GlobalIndexType index_i);

  /**
   * Perform cycle over all particles. For each of them load  all the neighbors and
   * perform the interactions.
   */
  void Interact();


  ParticleSystem particles;
  NeighborSearch neighborSearch;

  Variables vars;

private:

};

} // SPH
} // ParticleSystem
} // TNL

#include "SPH_impl.h"
