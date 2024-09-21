#pragma once

namespace TNL {
namespace SPH {
namespace ViscousTerms {

template< typename SPHCaseConfig >
class ArtificialViscosity
{
   public:
   using RealType = typename SPHCaseConfig::RealType;

   struct ParamsType
   {
     template< typename SPHState >
     __cuda_callable__
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

template< typename SPHCaseConfig >
class PhysicalViscosity
{
   public:
   using RealType = typename SPHCaseConfig::RealType;

   struct ParamsType
   {
     template< typename SPHState >
     __cuda_callable__
     ParamsType( SPHState sphState )
     : h( sphState.h ),
       viscoValue( sphState.dynamicViscosity ),
       preventZero( sphState.h * sphState.h * sphState.eps ) {}

     const RealType h;
     const RealType viscoValue;
     const RealType dimensionCoef = ( 2.f + SPHCaseConfig::spaceDimension ) * 2.f;
     const RealType preventZero;
   };

   __cuda_callable__
   static RealType
   Pi( const RealType& rhoI, const RealType& rhoJ, const RealType& drs, const RealType& drdv, const ParamsType& params )
   {
      const RealType viscoCoef = params.dimensionCoef * params.viscoValue / ( rhoI * rhoJ );
      return viscoCoef * drdv / ( drs * drs + params.preventZero );
   }
};

template< typename SPHCaseConfig >
class CombinedViscosity
{
   public:
   using RealType = typename SPHCaseConfig::RealType;

   struct ParamsType
   {
     template< typename SPHState >
     __cuda_callable__
     ParamsType( SPHState sphState )
     : h( sphState.h ),
       coefAV( ( -2.f ) * sphState.alpha * sphState.speedOfSound ),
       preventZero( sphState.h * sphState.h * sphState.eps ),
       kinematicViscosity( sphState.dynamicViscosity / sphState.rho0 ) {}

     const RealType h;
     const RealType coefAV;
     const RealType kinematicViscosity;
     const RealType preventZero;
   };

   __cuda_callable__
   static RealType
   Pi( const RealType& rhoI, const RealType& rhoJ, const RealType& drs, const RealType& drdv, const ParamsType& params )
   {
      const RealType mu = params.h * drdv / ( drs * drs + params.preventZero );
      return ( drdv < 0.f ) ? ( params.coefAV * mu / ( rhoI + rhoJ ) + mu * params.kinematicViscosity) : ( params.kinematicViscosity * mu );
   }
};

} // ViscousTerms
} // SPH
} // TNL

