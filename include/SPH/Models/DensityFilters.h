#pragma once

#include "../SPHTraits.h"

namespace TNL {
namespace SPH {
namespace DensityFilters {

class None
{
   public:

   template< typename FluidPointer, typename ModelParams >
   static void
   filterDensity( FluidPointer& fluid, ModelParams& modelParams ) {}
};

template< typename SPHCaseConfig, typename KernelFunction >
class ShepardFilter
{
   public:
   using DeviceType = typename SPHCaseConfig::DeviceType;
   using SPHTraitsType = SPHFluidTraits< SPHCaseConfig >;
   using LocalIndexType = typename SPHCaseConfig::LocalIndexType;
   using GlobalIndexType = typename SPHCaseConfig::GlobalIndexType;
   using RealType = typename SPHCaseConfig::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   template< typename Particles, typename FluidPointer, typename ModelParams >
   static void
   filterDensity( FluidPointer& fluid, ModelParams& modelParams )
   {
      typename Particles::NeighborsLoopParams searchInFluid( fluid->particles );

      const auto view_points = fluid->getParticles()->getPoints().getConstView();
      auto view_rho = fluid->getVariables()->rho.getView();
      const RealType searchRadius = fluid->getParticles()->getSearchRadius();
      const RealType m = modelParams.mass;
      const RealType h = modelParams.h;
      const RealType rho0 = modelParams.rho0;

      auto DensityFilter = [=] __cuda_callable__ ( LocalIndexType i,
                                                   LocalIndexType j,
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

      auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
      {
         const VectorType r_i = view_points[ i ];
         RealType rho_i = 0.f;
         RealType gamma_i = 0.f;
         Particles::NeighborsLoop::exec( i, r_i, searchInFluid, DensityFilter, &rho_i,  &gamma_i );

         if( gamma_i > 0.01f ){
            view_rho[ i ] = rho_i / gamma_i;
         }
         else
            view_rho[ i ] = rho0;
      };
      Algorithms::parallelFor< DeviceType >(
            fluid->getFirstActiveParticle(), fluid->getLastActiveParticle() + 1, particleLoop );
   }

   template< typename Particles, typename FluidPointer, typename ModelParams >
   static void
   filterDensityOverlaps( FluidPointer& fluid, ModelParams& modelParams )
   {
      typename Particles::NeighborsLoopParams searchInFluid( fluid->particles );

      const auto view_points = fluid->getParticles()->getPoints().getConstView();
      auto view_rho = fluid->getVariables()->rho.getView();
      auto view_gamma = fluid->getVariables()->gamma.getView();
      const RealType searchRadius = fluid->getParticles()->getSearchRadius();
      const RealType m = modelParams.mass;
      const RealType h = modelParams.h;
      const RealType rho0 = modelParams.rho0;

      auto DensityFilter = [=] __cuda_callable__ ( LocalIndexType i,
                                                   LocalIndexType j,
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

      auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
      {
         const VectorType r_i = view_points[ i ];
         RealType rho_i = 0.f;
         RealType gamma_i = 0.f;
         Particles::NeighborsLoop::exec( i, r_i, searchInFluid, DensityFilter, &rho_i,  &gamma_i );

         view_rho[ i ] = rho_i ;
         view_gamma[ i ] = gamma_i;
      };
      Algorithms::parallelFor< DeviceType >(
            fluid->getFirstActiveParticle(), fluid->getLastActiveParticle() + 1, particleLoop );

      //TODO: Remove this ridiculous periodic update
      if constexpr( SPHCaseConfig::numberOfPeriodicBuffers > 0 ) {
         for( long unsigned int i = 0; i < std::size( fluid->periodicPatches ); i++ ) {
            const auto zoneParticleIndices_view = fluid->periodicPatches[ i ]->particleZone.getParticlesInZone().getConstView();
            const GlobalIndexType numberOfZoneParticles = fluid->periodicPatches[ i ]->particleZone.getNumberOfParticles();
            const VectorType shift = fluid->periodicPatches[ i ]->config.shift;

            auto periodicParticleLoop = [=] __cuda_callable__ ( LocalIndexType j ) mutable
            {
               const GlobalIndexType p = zoneParticleIndices_view[ j ];
               const VectorType r_i = view_points[ p ] + shift;
               RealType rho_i = 0.f;
               RealType gamma_i = 0.f;
               Particles::NeighborsLoop::exec( p, r_i, searchInFluid, DensityFilter, &rho_i,  &gamma_i );

               view_rho[ p ] += rho_i;
               view_gamma[ p ] += gamma_i;
            };
            Algorithms::parallelFor< DeviceType >( 0, numberOfZoneParticles, periodicParticleLoop );
         }
      }

      // finalize interaction
      auto finalizeParticleLoop = [ = ] __cuda_callable__( LocalIndexType i ) mutable
      {
         const RealType gamma_i = view_gamma[ i ];
         const RealType rho_i = view_rho[ i ];
         if( gamma_i > 0.01f )
            view_rho[ i ] = view_rho[ i ] / gamma_i;
         else
            view_rho[ i ] = rho0;
      };
      TNL::Algorithms::parallelFor< DeviceType >(
            fluid->getFirstActiveParticle(), fluid->getLastActiveParticle() + 1, finalizeParticleLoop );
   }

};


} // ViscousTerms
} // SPH
} // TNL

