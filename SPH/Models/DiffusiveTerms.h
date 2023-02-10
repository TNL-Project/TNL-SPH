#pragma once

namespace TNL {
namespace ParticleSystem {
namespace SPH {

class NoneDiffusive
{

};

template< typename SPHCaseConfig >
class MolteniDiffusiveTerm
{
   public:
   using RealType = typename SPHCaseConfig::RealType; //fix this

   static constexpr RealType coefDT = ( 2.f ) * SPHCaseConfig::h * SPHCaseConfig::delta * SPHCaseConfig::speedOfSound;

   __cuda_callable__
   static RealType
   Psi( const RealType& rhoI, const RealType& rhoJ, const RealType& drs )
   {
      return coefDT * ( rhoJ - rhoI ) / ( drs * drs );
   }
};

} // SPH
} // ParticleSystem
} // TNL

