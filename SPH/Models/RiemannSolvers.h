#pragma once

namespace TNL {
namespace ParticleSystem {
namespace SPH {

#define MIN( a, b ) ( ( ( a ) < ( b ) ) ? ( a ) : ( b ) ) //TNL::min with speed of sounds returns error.

template< typename SPHCaseConfig >
class RiemanSolverLinearized
{
public:
   using RealType = typename SPHCaseConfig::RealType; //fix this

   static constexpr RealType eta = SPHCaseConfig::etaLimiter;
   static constexpr RealType speedOfSound = SPHCaseConfig::speedOfSound;

   __cuda_callable__
   static RealType
   statePressure( const RealType& vL, const RealType& vR, const RealType& pL, const RealType& pR, const RealType& rhoAverage )
   {
      //const RealType beta = TNL::min( eta * max( vL - vR, 0 ), speedOfSound ); // limiter
      const RealType beta = MIN( eta * max( vL - vR, 0.f ), speedOfSound ); //limiter
      return 0.5f * ( pL + pR ) + 0.5f * ( vL - vR ) * rhoAverage * beta;
   }

   __cuda_callable__
   static RealType
   stateVelocity( const RealType& vL, const RealType& vR, const RealType& pL, const RealType& pR, const RealType& rhoAverage )
   {
      return 0.5f * ( vR + vL ) + 0.5f * ( pL - pR ) / ( rhoAverage * speedOfSound );
   }
};

} // SPH
} // ParticleSystem
} // TNL

