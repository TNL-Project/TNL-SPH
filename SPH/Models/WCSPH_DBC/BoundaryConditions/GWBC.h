#pragma once

#include "../Variables.h"
#include "../../../SPHTraits.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {
namespace WCSPH {
namespace BoundaryConditions {

template< typename SPHConfig >
class GWBCVariables : public SPHFluidVariables< SPHConfig >
{
   public:
   using BaseType = SPHFluidVariables< SPHConfig >;
   using GlobalIndexType = typename BaseType::GlobalIndexType;

   GWBCVariables( GlobalIndexType size )
   : SPHFluidVariables< SPHConfig >( size ) {};
};

template< typename SPHConfig >
class GWBC
{
   public:

};

} // BoundaryConditions
} // WSPH
} // SPH
} // ParticleSystem
} // TNL

