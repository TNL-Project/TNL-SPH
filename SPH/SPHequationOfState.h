#pragma once

#include "SPHFluidTraits.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

class TaitWeaklyCompressibleEOS
{
  public:
  using RealType = float; //fix this

  __cuda_callable__
  static RealType
  DensityToPressure( RealType rho )
  {
    const RealType gamma = 7.;
    const RealType rho0 = 1000.;
    const RealType relativeDensity = rho / rho0;
    const RealType coefB = 1107129;

    return  coefB * ( pow( relativeDensity, gamma ) - 1 );
  }

};

} // SPH
} // ParticleSystem
} // TNL
