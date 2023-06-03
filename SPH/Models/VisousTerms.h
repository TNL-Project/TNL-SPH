#pragma once

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHCaseConfig >
class ArtificialViscosity
{
   public:
   using RealType = typename SPHCaseConfig::RealType; //fix this

   struct ParamsType
   {
     template< typename SPHState >
     ParamsType( SPHState sphState )
     : h( sphState.h ),
       coefAV( ( -2.f ) * sphState.alpha * sphState.speedOfSound ),
       preventZero( sphState.h * sphState.h * sphState.eps ) {}

     const RealType h;
     const RealType coefAV;
     const RealType preventZero;
   };

   __cuda_callable__
   static RealType
   Pi( const RealType& rhoI, const RealType& rhoJ, const RealType& drs, const RealType& drdv, const ParamsType& params )
   {
      const RealType mu = params.h * drdv / ( drs * drs + params.preventZero );
      return ( drdv < 0.f ) ? ( params.coefAV * mu / ( rhoI + rhoJ ) ) : ( 0.f );
   }
};

} // SPH
} // ParticleSystem
} // TNL

