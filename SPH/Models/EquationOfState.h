#pragma once

namespace TNL {
namespace ParticleSystem {
namespace SPH {

//template< typename SPHCaseConfig >
class TaitWeaklyCompressibleEOS
{
  public:
  using RealType = float; //fix this

  __cuda_callable__
  static RealType
  DensityToPressure( RealType rho )
  {
    const RealType gamma = 7.; //fix this
    const RealType rho0 = 1000.; //fix this
    const RealType coefB = 1107129; //fix this

    const RealType relativeDensity = rho / rho0;

    return  coefB * ( pow( relativeDensity, gamma ) - 1 );
  }
};

} // SPH
} // ParticleSystem
} // TNL
