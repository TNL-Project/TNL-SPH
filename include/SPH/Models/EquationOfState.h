#pragma once

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHCaseConfig >
class TaitWeaklyCompressibleEOS
{
public:
   using RealType = typename SPHCaseConfig::RealType; //fix this

   struct ParamsType
   {
     template< typename SPHState >
     ParamsType( SPHState sphState )
     : rho0( sphState.rho0 ),
       coefB( sphState.coefB ) {}

     const RealType rho0;
     const RealType coefB;
   };

   __cuda_callable__
   static RealType
   DensityToPressure( const RealType& rho, const ParamsType& params )
   {
      const RealType gamma = 7.f;
      const RealType relativeDensity = rho / params.rho0;
      return params.coefB * ( powf( relativeDensity, gamma ) - 1.f );
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
     ParamsType( SPHState sphState )
     : rho0( sphState.rho0 ),
       c0( sphState.speedOfSound ) {}

     const RealType rho0;
     const RealType c0;
   };

   __cuda_callable__
   static RealType
   DensityToPressure( const RealType& rho, const ParamsType& params )
   {
      const RealType c0 = params.c0;
      const RealType rho0 = params.rho0;
      return c0 * c0 * ( rho - rho0 );
   }
};

} // SPH
} // ParticleSystem
} // TNL

