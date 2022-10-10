#pragma once

#include "../Particles/Particles.h"

#include "SPHFluidTraits.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template < typename Variables, typename ParticleSystem >
class SPHSimulation
{
public:


  SPHSimulation() = default;

  SPHSimulation(int size, float h)
  : particles(size, h), vars(size) {};


//protected:

  //INTERACET

  //Particles

  ParticleSystem particles;
  Variables vars;

private:

};

} // SPH
} // ParticleSystem
} // TNL
