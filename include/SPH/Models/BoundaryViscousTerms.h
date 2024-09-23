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
       scaleBVTCoef( sphState.scaleBVTCoef ),
       m( sphState.mass ) {}

     const RealType searchRadius;
     const RealType dynamicViscosity;
     const RealType scaleBVTCoef;
     const RealType m;
   };

   __cuda_callable__
   static VectorType
   Xi( const VectorType& r_ik,  const VectorType& v_ik, const VectorType& n_k, const ParamsType& params )
   {
      //FIXME: Normal has to be defined as particle filed passed by initial boundary conditoo
      const VectorType t_k = { -1.f, 0.f };

      const RealType intersection = 2 * sqrt( pow( params.searchRadius, 2 ) - pow( ( n_k, r_ik ), 2 ) );
      const RealType normalDerivative = ( -1.0f ) * ( v_ik, t_k ) / ( n_k, r_ik );

      return ( params.scaleBVTCoef * params.dynamicViscosity * intersection * normalDerivative / params.m ) * t_k;
   }
};

template< typename SPHCaseConfig >
class PhysicalViscosity_MVT
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
     : dynamicViscosity( sphState.dynamicViscosity ),
       preventZero( sphState.h * sphState.h * sphState.eps ) {}

     const RealType dynamicViscosity;
     const RealType preventZero;
   };

   __cuda_callable__
   static VectorType
   Xi( const RealType& drs,
       const VectorType& r_ik,
       const VectorType v_ik,
       const RealType& W,
       const VectorType& n_k,
       const RealType& s_k,
       const ParamsType& params )
   {
      return params.dynamicViscosity ds_k * v_ik / ( r_ik, n_k ) * W * s_k;
   }
};

} // BoundaryViscousTerms
} // SPH
} // TNL

