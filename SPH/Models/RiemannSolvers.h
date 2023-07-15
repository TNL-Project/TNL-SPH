#pragma once

namespace TNL {
namespace ParticleSystem {
namespace SPH {

#define MIN( a, b ) ( ( ( a ) < ( b ) ) ? ( a ) : ( b ) ) //TNL::min with speed of sounds returns error.

template< typename SPHCaseConfig >
class RiemanSolverLinearized
{
   public:
   using RealType = typename SPHCaseConfig::RealType;

   static constexpr RealType eta = SPHCaseConfig::etaLimiter;
   static constexpr RealType speedOfSound = SPHCaseConfig::speedOfSound;

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
      //const RealType beta = TNL::min( eta * max( vL - vR, 0 ), speedOfSound );
      const RealType beta = MIN( params.eta * max( vL - vR, 0.f ), params.speedOfSound );
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

} // SPH
} // ParticleSystem
} // TNL

