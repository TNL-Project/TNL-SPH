#pragma once

namespace TNL {
namespace SPH {
namespace DiffusiveTerms {

/**
 * \brief Template for disabled diffusive term.
 */
template< typename SPHCaseConfig >
class NoneDiffusiveTerm
{
   public:
   using RealType = typename SPHCaseConfig::RealType;

   struct ParamsType
   {
     template< typename SPHState >
     __cuda_callable__
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
 * \brief Diffusive term proposed by Molteni & Colagrossi (Molteni & Colagrossi, 2009).
 *
 * \tparam SPHCaseConfig is a default config definig all data types.
 */
template< typename SPHCaseConfig >
class MolteniDiffusiveTerm
{
   public:
   using RealType = typename SPHCaseConfig::RealType;

   struct ParamsType
   {
     template< typename SPHState >
      __cuda_callable__
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

} // DiffusiveTerms
} // SPH
} // TNL

