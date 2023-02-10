#pragma once

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHCaseConfig >
class ArtificialViscosity
{
public:
   using RealType = typename SPHCaseConfig::RealType; //fix this

   static constexpr RealType coefAV = ( -2.f ) * SPHCaseConfig::alpha * SPHCaseConfig::speedOfSound; //?
   static constexpr RealType h = SPHCaseConfig::h;
   static constexpr RealType epsilon = SPHCaseConfig::eps;
   static constexpr RealType preventZero = h * h * epsilon;

   __cuda_callable__
   static RealType
   Pi( const RealType& rhoI, const RealType& rhoJ, const RealType& drs, const RealType& drdv )
   {
      const RealType mu = h * drdv / ( drs * drs + preventZero );
      return ( drdv < 0.f ) ? ( coefAV * mu / ( rhoI + rhoJ ) ) : ( 0.f );
   }
};

} // SPH
} // ParticleSystem
} // TNL

