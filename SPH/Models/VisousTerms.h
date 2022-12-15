#pragma once

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHCaseConfig >
class ArtificialViscosity
{
public:
   using RealType = typename SPHCaseConfig::RealType; //fix this

   static constexpr RealType coefAV = ( -2. ) * SPHCaseConfig::alpha * SPHCaseConfig::speedOfSound; //?
   static constexpr RealType h = SPHCaseConfig::h;
   static constexpr RealType epsilon = SPHCaseConfig::eps;

   __cuda_callable__
   static RealType
   Pi( RealType rhoI, RealType rhoJ, RealType drs, RealType drdv )
   {
      const RealType mu = h * drdv / ( drs * drs + epsilon * h * h );
      const RealType pi = ( drdv < 0 ) ? ( coefAV * mu / ( rhoI + rhoJ ) ) : (0.);
      return pi;
   }
};

} // SPH
} // ParticleSystem
} // TNL

