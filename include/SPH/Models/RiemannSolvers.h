#pragma once

#include <algorithm>
namespace TNL {
namespace SPH {
namespace RiemannSolvers {

template< typename SPHCaseConfig >
class RoeLinearized
{
   public:
   using RealType = typename SPHCaseConfig::RealType;

   struct ParamsType
   {
     template< typename SPHState >
     ParamsType( SPHState sphState )
     : eta( sphState.limiterEta ), speedOfSound( sphState.speedOfSound ) {}

     const RealType eta;
     const RealType speedOfSound;
   };

   __cuda_callable__
   static RealType
   statePressure( const RealType& vL,
                  const RealType& vR,
                  const RealType& pL,
                  const RealType& pR,
                  const RealType& rhoAverage,
                  const ParamsType& params )
   {
      const RealType beta = std::min( params.eta * max( vL - vR, 0.f ), params.speedOfSound );
      return 0.5f * ( pL + pR ) + 0.5f * ( vL - vR ) * rhoAverage * beta;
   }

   __cuda_callable__
   static RealType
   stateVelocity( const RealType& vL,
                  const RealType& vR,
                  const RealType& pL,
                  const RealType& pR,
                  const RealType& rhoAverage,
                  const ParamsType& params )
   {
      return 0.5f * ( vR + vL ) + 0.5f * ( pL - pR ) / ( rhoAverage * params.speedOfSound );
   }
};

} // RiemannSolvers
} // SPH
} // TNL

