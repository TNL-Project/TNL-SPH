#pragma once

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHCaseConfig >
class RiemanSolverLinearized
{
public:
   using RealType = typename SPHCaseConfig::RealType; //fix this

   static constexpr eta = SPHCaseConfig::etaLimiter;
   static constexpr speedOfSound = SPHCaseConfig::speedOfSound;

   __cuda_callable__
   static RealType
   statePressure( const RealType& vL, const RealType& vR, const RealType& pL, const RealType& pR, const RealType& rhoAverage )
   {
      const RealType beta = min( eta * max( vL - vR, 0 ), c ); //limiter
      return 0.5 * ( pL + pR ) + 0.5 * ( vL - vR ) * rhoAvg * beta;
   }

   __cuda_callable__
   static RealType
   stateVelocity( const RealType& vL, const RealType& vR, const RealType& pL, const RealType& pR, const RealType& rhoAverage )
   {
      return 0.5 * ( vR + vL ) + 0.5 * ( pL - pR ) / ( rhoAvg * c );
   }
};

} // SPH
} // ParticleSystem
} // TNL

