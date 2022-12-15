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

   static constexpr RealType coefDT = 2 * SPHCaseConfig::h * SPHCaseConfig::delta * SPHCaseConfig::speedOfSound;

   __cuda_callable__
   static RealType
   Psi( RealType rhoI, RealType rhoJ, RealType drs )
   {
      const RealType psi = coefDT * ( rhoJ - rhoI ) / ( drs * drs );
      return psi;
   }
};

} // SPH
} // ParticleSystem
} // TNL

