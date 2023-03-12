#pragma once

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHCaseConfig >
class TaitWeaklyCompressibleEOS
{
public:
   using RealType = typename SPHCaseConfig::RealType; //fix this

   static constexpr RealType rho0 = SPHCaseConfig::rho0;
   static constexpr RealType coefB = SPHCaseConfig::coefB;

   __cuda_callable__
   static RealType
   DensityToPressure( RealType rho )
   {
      const RealType gamma = 7.f;
      const RealType relativeDensity = rho / rho0;
      return coefB * ( powf( relativeDensity, gamma ) - 1.f );
   }
};

} // SPH
} // ParticleSystem
} // TNL

