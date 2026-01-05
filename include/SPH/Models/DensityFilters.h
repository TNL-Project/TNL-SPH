#pragma once

#include "../SPHTraits.h"

namespace TNL {
namespace SPH {
namespace DensityFilters {

template< typename ParticlesType, typename ModelConfig >
class None
{
public:

   template< typename FluidPointer, typename ModelParams >
   static void
   filterDensity( FluidPointer& fluid, ModelParams& modelParams )
   {}
};

template< typename ParticlesType, typename ModelConfig >
class ShepardFilter
{
public:

   using SPHTraitsType = SPHFluidTraits< typename ModelConfig::SPHConfig >;
   using IndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;
   using KernelFunction = typename ModelConfig::KernelFunction;

   template< typename FluidPointer, typename ModelParams >
   static void
   filterDensity( FluidPointer& fluid, ModelParams& modelParams )
   {
      auto searchInFluid = fluid->getParticles()->getSearchToken( fluid->getParticles() );

      const auto view_points = fluid->getParticles()->getPoints().getConstView();
      auto view_rho = fluid->getVariables()->rho.getView();
      const RealType searchRadius = fluid->getParticles()->getSearchRadius();
      const RealType m = modelParams.mass;
      const RealType h = modelParams.h;
      const RealType rho0 = modelParams.rho0;

      auto DensityFilter = [=] __cuda_callable__ ( IndexType i,
                                                   IndexType j,
                                                   VectorType& r_i,
                                                   RealType* rho_i,
                                                   RealType* gamma_i ) mutable
      {
         const VectorType r_j = view_points[ j ];
         const VectorType r_ij = r_i - r_j;
         const RealType drs = l2Norm( r_ij );
         if (drs <= searchRadius )
         {
            const RealType rho_j = view_rho[ j ];
            const RealType W = KernelFunction::W( drs, h );

            *rho_i += W * m;
            *gamma_i += W * m / rho_j;
         }
      };

      auto particleLoop = [=] __cuda_callable__ ( IndexType i ) mutable
      {
         const VectorType r_i = view_points[ i ];
         RealType rho_i = 0.f;
         RealType gamma_i = 0.f;
         ParticlesType::NeighborsLoop::exec( i, r_i, searchInFluid, DensityFilter, &rho_i,  &gamma_i );

         if( gamma_i > 0.01f ){
            view_rho[ i ] = rho_i / gamma_i;
         }
         else
            view_rho[ i ] = rho0;
      };
      fluid->getParticles()->forAll( particleLoop );
   }
};


} // DensityFilter
} // SPH
} // TNL

