#pragma once

namespace TNL {
namespace ParticleSystem {


class Particle
{
public:
   getIndex();

   template< typename Function, typename... FunctionArgs >
   loopOverNeighbors();

protected:
   const ParticlesPointer* particlesPointer = nullptr;
};

} //namespace Particles
} //namespace TNL

