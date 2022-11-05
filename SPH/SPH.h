#pragma once

#include <TNL/Containers/Vector.h>
#include <TNL/Algorithms/reduce.h>

#include "../Particles/Particles.h"

#include "Kernels.h"

//Model
#include "SPHTraits.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Model, typename ParticleSystem, typename NeighborSearch >
class SPHSimulation
{
public:

  using DeviceType = typename ParticleSystem::Device;

  using LocalIndexType = typename ParticleSystem::LocalIndexType;
  using GlobalIndexType = typename ParticleSystem::GlobalIndexType;
  using RealType = typename ParticleSystem::RealType;
  using InteractionResultType = typename Model::InteractionResultType;

  SPHSimulation() = default;

  SPHSimulation( GlobalIndexType size, RealType h, GlobalIndexType numberOfCells)
  : particles( size, h ), neighborSearch( particles, numberOfCells ), model( size, particles.getPoints(), particles ) {};

//protected:

  /**
   * Perform neighbors search and fill neighborsList in Particle system variable.
   */
  void PerformNeighborSearch();

  /**
   * Proces one particle (i.e. loop over all its neighbors and perform interactions).
   */
  void ProcessOneParticle( GlobalIndexType index_i );

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
