#pragma once

#include "../../SPHTraits.h"
#include "details.h"

namespace TNL {
namespace SPH {
namespace Stress {

template< typename SPHCaseConfig >
class FluidStress
{
public:
   using TraitsType = SPHFluidTraits< SPHCaseConfig >;
   using RealType = typename SPHCaseConfig::RealType;
   using VectorType = typename TraitsType::VectorType;
   using MatrixType = typename TraitsType::MatrixType;

   struct ParamsType
   {
     template< typename SPHState >
     __cuda_callable__
     ParamsType( SPHState sphState )
     : rho0( sphState.rho0 ),
       cl( sphState.speedOfSound_bulk ),
       cs( sphState.speedOfSound_shear ),
       antiClumpFact( sphState.antiClumpFactor ){}

     const RealType rho0;
     const RealType cs;
     const RealType cl;
     const RealType antiClumpFact;
   };

   __cuda_callable__
   static MatrixType
   distortionToStress( const MatrixType& A, const RealType& rho, const ParamsType& params )
   {
      const RealType cl = params.cl;
      const RealType cs = params.cs;
      const RealType rho0 = params.rho0;
      const RealType antiClumpFact = params.antiClumpFact;

      const MatrixType E = details::template unitMatrix< MatrixType >(); //FIXME: Do unit matrix constructor!
      const MatrixType tAA = A * transpose( A );

      return powf( cl, 2 ) * ( rho - rho0 / ( 1.f + antiClumpFact ) ) * E + powf( cs, 2 ) * rho * details::deviator( tAA ) * tAA;

   }
};

} // Stress
} // SPH
} // TNL

