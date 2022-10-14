#pragma once

#include <TNL/Containers/Vector.h>
#include <TNL/Algorithms/reduce.h>

#include "../Particles/Particles.h"

#include "SPHFluidTraits.h"
#include "SPHFluidVariables.h"
#include "SPHInteractions.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Model, typename ParticleSystem, typename NeighborSearch >
class SPHSimulation
{
public:

  using LocalIndexType = typename ParticleSystem::LocalIndexType;
  using GlobalIndexType = typename ParticleSystem::GlobalIndexType;
  using DeviceType = typename ParticleSystem::Device;

  using RealType = typename ParticleSystem::RealType;

  SPHSimulation() = default;

  SPHSimulation( GlobalIndexType size, float h )
  : particles( size, h ), neighborSearch( particles, 100 ), model( size, particles.getPoints(), particles ) {};

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
  Model model;

private:

};

} // SPH
} // ParticleSystem
} // TNL

#include "SPH_impl.h"
