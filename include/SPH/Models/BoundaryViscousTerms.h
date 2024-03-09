#pragma once

#include "../SPHTraits.h"

namespace TNL {
namespace SPH {
namespace BoundaryViscousTerms {


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
     ParamsType( SPHState sphState ) {}
   };

   __cuda_callable__
   static VectorType
   Xi( const VectorType& r_ik,  const VectorType& v_ik, const VectorType& n_k, const ParamsType& params  )
   {
      return 0.f;
   }

};

template< typename SPHCaseConfig >
class NewtonViscousLaw
{
   public:
   using SPHTraitsType = SPHFluidTraits< SPHCaseConfig >;
   using RealType = typename SPHCaseConfig::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   struct ParamsType
   {
     template< typename SPHState >
     ParamsType( SPHState sphState )
     : searchRadius( sphState.searchRadius ),
       dynamicViscosity( sphState.dynamicViscosity ),
       m( sphState.mass ) {}

     const RealType searchRadius;
     const RealType dynamicViscosity;
     const RealType m;
   };

   __cuda_callable__
   static VectorType
   Xi( const VectorType& r_ik,  const VectorType& v_ik, const VectorType& n_k, const ParamsType& params )
   {
      const VectorType t_k = { -1.f, 0.f };
      //if( n_k[ 1 ] > 0.5f )
      //   t_k = { -n_k[ 1 ], n_k[ 0 ] };
      //else if( n_k[ 1 ] < -0.5f )
      //   t_k = { n_k[ 1 ], -n_k[ 0 ] };

      const RealType intersection = 2 * sqrt( pow( params.searchRadius, 2 ) - pow( ( n_k, r_ik ), 2 ) );
      const RealType normalDerivative = ( -1.0f ) * ( v_ik, t_k ) / ( n_k, r_ik );

      return ( params.dynamicViscosity * intersection * normalDerivative / params.m ) * t_k;
   }
};

} // BoundaryViscousTerms
} // SPH
} // TNL

