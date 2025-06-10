#pragma once

#include "../SPHTraits.h"

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

namespace BIViscousTerms {

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
   static VectorType
   Pi( const RealType& drs,
       const VectorType& r_ij,
       const VectorType& v_ij,
       const RealType& rho_i,
       const RealType& rho_j,
       const VectorType& gradWV_j,
       const ParamsType& params )
   {
      const VectorType zeroVector = 0.f;
      return zeroVector;
   }

   __cuda_callable__
   static VectorType
   BI_Pi( const RealType& drs,
          const VectorType& r_ik,
          const VectorType v_ik,
          const RealType& rho_i,
          const RealType& rho_j,
          const VectorType& n_k,
          const RealType& WS_k,
          const ParamsType& params )
   {
      const VectorType zeroVector = 0.f;
      return zeroVector;
   }
};

template< typename SPHCaseConfig >
class ArtificialViscosity
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
     : h( sphState.h ),
       coefAV( ( 2.f ) * sphState.alpha * sphState.speedOfSound ),
       preventZero( sphState.h * sphState.h * sphState.eps ) {}

     const RealType h;
     const RealType coefAV;
     const RealType preventZero;
   };

   __cuda_callable__
   static VectorType
   Pi( const RealType& drs,
       const VectorType& r_ij,
       const VectorType& v_ij,
       const RealType& rho_i,
       const RealType& rho_j,
       const VectorType& gradWV_j,
       const ParamsType& params )
   {
      const RealType drdv = ( r_ij, v_ij );
      const RealType mu = params.h * drdv / ( drs * drs + params.preventZero );
      const VectorType gradWm_j = rho_j * gradWV_j;
      const VectorType zeroVector = 0.f;
      return ( drdv < 0.f ) ? ( params.coefAV * mu / ( rho_i + rho_j ) * gradWm_j ) : ( zeroVector );
   }

   __cuda_callable__
   static VectorType
   BI_Pi( const RealType& drs,
          const VectorType& r_ik,
          const VectorType v_ik,
          const RealType& rho_i,
          const RealType& rho_j,
          const VectorType& n_k,
          const RealType& WS_ik,
          const ParamsType& params )
   {
      const RealType drdv = ( r_ik, v_ik );
      const RealType mu = params.h * drdv / ( drs * drs + params.preventZero );
      const VectorType zeroVector = 0.f;
      return ( drdv < 0.f ) ? ( params.coefAV * mu / ( rho_i + rho_j ) * rho_j * n_k * WS_ik ) : ( zeroVector );

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
   Pi( const RealType& drs,
       const VectorType& r_ij,
       const VectorType& v_ij,
       const RealType& rho_i,
       const RealType& rho_j,
       const VectorType& gradWV_j,
       const ParamsType& params )
   {
      const RealType viscoCoef = params.dynamicViscosity / rho_i;
      return 2.f * viscoCoef * v_ij * ( r_ij, gradWV_j ) / ( drs * drs + params.preventZero );
   }

   __cuda_callable__
   static VectorType
   BI_Pi( const RealType& drs,
          const VectorType& r_ik,
          const VectorType v_ik,
          const RealType& rho_i,
          const RealType& rho_j,
          const VectorType& n_k,
          const RealType& WS_ik,
          const ParamsType& params )
   {
      const RealType viscoCoef = params.dynamicViscosity / rho_i;
      return 2.f * viscoCoef * v_ik / ( r_ik, n_k ) * WS_ik;
   }
};

template< typename SPHCaseConfig >
class PhysicalViscosity_MGVT
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
     : dp( sphState.dp ),
       dynamicViscosity( sphState.dynamicViscosity ),
       preventZero( sphState.h * sphState.h * sphState.eps ) {}

     const RealType dp;
     const RealType dynamicViscosity;
     const RealType dimensionCoef = ( 2.f + SPHCaseConfig::spaceDimension ) * 2.f;
     const RealType preventZero;
   };

   __cuda_callable__
   static VectorType
   Pi( const RealType& drs,
       const VectorType& r_ij,
       const VectorType& v_ij,
       const RealType& rho_i,
       const RealType& rho_j,
       const VectorType& gradWV_j,
       const ParamsType& params )
   {
      const RealType viscoCoef = params.dimensionCoef * params.dynamicViscosity / rho_i;
      const RealType pi = ( r_ij, v_ij ) / ( drs * drs + params.preventZero );

      return viscoCoef * pi * gradWV_j;
   }

   __cuda_callable__
   static VectorType
   BI_Pi( const RealType& drs,
          const VectorType& r_ik,
          const VectorType v_ik,
          const RealType& rho_i,
          const RealType& rho_j,
          const VectorType& n_k,
          const RealType& WS_ik,
          const ParamsType& params )
   {
      const RealType viscoCoef = params.dimensionCoef * params.dynamicViscosity / rho_i;

      const RealType pi = ( r_ik, v_ik ) / ( drs * drs + params.preventZero );
      const RealType r_ik_n = std::max( std::abs( ( r_ik, n_k ) ), params.dp );
      const VectorType v_ik_t = v_ik - ( v_ik, n_k ) * n_k;

      const VectorType lap_v_n = pi * rho_j * n_k * WS_ik;
      const VectorType lap_v_t = 2.f * v_ik_t / ( r_ik_n ) * WS_ik;

      return viscoCoef * ( lap_v_n + lap_v_t );
   }
};

} // BIViscousTerms

} // SPH
} // TNL

