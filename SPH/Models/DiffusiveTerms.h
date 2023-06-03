#pragma once

namespace TNL {
namespace ParticleSystem {
namespace SPH {

/**
 * \brief Template for disabled diffusive term.
 */

template< typename SPHCaseConfig >
class NoneDiffusiveTerm
{
   public:
   using RealType = typename SPHCaseConfig::RealType; //fix this

   struct ParamsType
   {
     template< typename SPHState >
     ParamsType( SPHState sphState ) {}
   };

   __cuda_callable__
   static RealType
   Psi( const RealType& rhoI, const RealType& rhoJ, const RealType& drs, const ParamsType& params )
   {
      return 0.f;
   }

};

/**
 * \brief Diffusive term proposed by Molteni & Colagrossi.
 */
template< typename SPHCaseConfig >
class MolteniDiffusiveTerm
{
   public:
   using RealType = typename SPHCaseConfig::RealType; //fix this

   struct ParamsType
   {
     template< typename SPHState >
     ParamsType( SPHState sphState )
     : coefDT( ( 2.f ) * sphState.h * sphState.delta * sphState.speedOfSound ) {}

     const RealType coefDT;
   };

   __cuda_callable__
   static RealType
   Psi( const RealType& rhoI, const RealType& rhoJ, const RealType& drs, const ParamsType& params )
   {
      return params.coefDT * ( rhoJ - rhoI ) / ( drs * drs );
   }
};

} // SPH
} // ParticleSystem
} // TNL

