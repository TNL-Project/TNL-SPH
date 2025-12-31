#pragma once

#include "../SPHTraits.h"

namespace TNL {
namespace SPH {
namespace EquationsOfState {

/**
 * \brief Template for disabled equation of state.
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
   DensityToPressure( const RealType& rho, const ParamsType& params )
   {
      return 0.f;
   }

   __cuda_callable__
   static RealType
   pressureToDensity( const RealType& p, const ParamsType& params )
   {
      return 0.f;
   }
};

template< typename SPHCaseConfig >
class TaitWeaklyCompressibleEOS
{
public:
   using RealType = typename SPHCaseConfig::RealType; //fix this

   struct ParamsType
   {
     template< typename SPHState >
      __cuda_callable__
     ParamsType( SPHState sphState )
     : rho0( sphState.rho0 ),
       coefB( sphState.coefB ),
       p0( sphState.p0 ){}

     const RealType rho0;
     const RealType coefB;
     const RealType p0;
   };

   __cuda_callable__
   static RealType
   DensityToPressure( const RealType& rho, const ParamsType& params )
   {
      const RealType gamma = 7.f;
      const RealType relativeDensity = rho / params.rho0;
      return params.coefB * ( powf( relativeDensity, gamma ) - 1.f ) + params.p0;
   }

   __cuda_callable__
   static RealType
   pressureToDensity( const RealType& p, const ParamsType& params )
   {
      const RealType gamma = 7.f;
      const RealType relativeDensity = powf( ( ( p - params.p0 ) / params.coefB ) + 1, 1.f / gamma );
      return params.rho0 * relativeDensity;
   }
};

template< typename SPHCaseConfig >
class TaitLinearizedWeaklyCompressibleEOS
{
public:
   using RealType = typename SPHCaseConfig::RealType; //fix this

   struct ParamsType
   {
     template< typename SPHState >
     __cuda_callable__
     ParamsType( SPHState sphState )
     : rho0( sphState.rho0 ),
       c0( sphState.speedOfSound ),
       p0 ( sphState.p0 ){}

     const RealType rho0;
     const RealType c0;
     const RealType p0;
   };

   __cuda_callable__
   static RealType
   DensityToPressure( const RealType& rho, const ParamsType& params )
   {
      const RealType c0 = params.c0;
      const RealType rho0 = params.rho0;
      const RealType p0 = params.p0;
      return c0 * c0 * ( rho - rho0 ) + p0;
   }

   __cuda_callable__
   static RealType
   pressureToDensity( const RealType& p, const ParamsType& params )
   {
      const RealType c0 = params.c0;
      const RealType rho0 = params.rho0;
      const RealType p0 = params.p0;
      return ( p - p0 ) / ( c0 * c0 ) + rho0;
   }
};


} // EquationsOfState
} // SPH
} // TNL

