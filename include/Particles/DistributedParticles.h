#include "GhostZone.h"

namespace TNL {
namespace SPH {

template< typename ParticleSystem >
class DistributedParticleSystem
{
public:

   using RealType = typename ParticleSystem::RealType;
   using IndexArrayType = typename ParticleSystem::IndexArrayType;
   using ParticleZone = ParticleZone< typename ParticleSystem::Config >;
   using GridType = Grid< typename ParticleSystem::Config::spaceDimension, RealType, Device, Index >
   using DistributedGridType = DistributedMesh< GridType >;

protected:


   Containers::Array< GlobalIndexType, Devices::Host, int > innerOverlaps;
   IndexArrayType innerOverlapsLinearized;


};

}  //namespace SPH
}  //namespace TNL

