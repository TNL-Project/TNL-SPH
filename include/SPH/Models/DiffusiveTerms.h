#pragma once

#include "../SPHTraits.h"

namespace TNL {
namespace SPH {
namespace DiffusiveTerms {

/**
 * \brief Template for disabled diffusive term.
 */
template< typename SPHCaseConfig >
class None
{
   public:
   using SPHTraitsType = SPHFluidTraits< SPHCaseConfig >;
   using RealType = typename SPHCaseConfig::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   struct ParamsType
   {
     template< typename SPHState >
     __cuda_callable__
     ParamsType( SPHState sphState ) {}
   };

   __cuda_callable__
   static RealType
   Psi( const RealType& rhoI, const RealType& rhoJ, const VectorType& r_ij, const RealType& drs, const ParamsType& params )
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
   using SPHTraitsType = SPHFluidTraits< SPHCaseConfig >;
   using RealType = typename SPHCaseConfig::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

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
   Psi( const RealType& rhoI, const RealType& rhoJ, const VectorType& r_ij, const RealType& drs, const ParamsType& params )
   {
      return params.coefDT * ( rhoJ - rhoI ) / ( drs * drs );
   }
};

/**
 * \brief Diffusive term proposed by Fourtakas (Fourtakas, 2019).
 *
 * \tparam SPHCaseConfig is a default config definig all data types.
 */
template< typename SPHCaseConfig >
class FourtakasDiffusiveTerm
{
   public:
   using SPHTraitsType = SPHFluidTraits< SPHCaseConfig >;
   using RealType = typename SPHCaseConfig::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   struct ParamsType
   {
     template< typename SPHState >
     __cuda_callable__
     ParamsType( SPHState sphState )
     : gravity( sphState.gravity ),
       rho0( sphState.rho0 ),
       coefB( sphState.coefB ),
       coefDT( ( 2.f ) * sphState.h * sphState.delta * sphState.speedOfSound ) {}

     const VectorType gravity;
     const RealType rho0;
     const RealType coefDT;
     const RealType coefB;
   };

   __cuda_callable__
   static RealType
   Psi( const RealType& rhoI, const RealType& rhoJ, const VectorType& r_ij, const RealType& drs, const ParamsType& params )
   {
      const RealType gamma = 7.f;
      const RealType gravMagnitude = l2Norm( params.gravity );
      const VectorType gravDirection = params.gravity / gravMagnitude;

      const RealType r_z_ij = ( r_ij, gravDirection );
      const RealType p_hydrostatic_ij = 1.f + r_z_ij * params.rho0 * gravMagnitude / params.coefB;
      const RealType rho_hydrostatic_ij = params.rho0 * powf( p_hydrostatic_ij, 1.f / gamma ) - params.rho0;

      return params.coefDT * ( ( rhoJ - rhoI ) - rho_hydrostatic_ij ) / ( drs * drs );
   }
};

} // DiffusiveTerms
} // SPH
} // TNL

