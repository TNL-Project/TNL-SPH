#pragma once

#include "../SPHTraits.h"

namespace TNL {
namespace SPH {
namespace VelocityFilters {

class None
{
   public:

   template< typename FluidPointer, typename ModelParams >
   static void
   filterVelocity( FluidPointer& fluid, ModelParams& modelParams ) {}
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
   filterVelocity( FluidPointer& fluid, ModelParams& modelParams )
   {
      typename Particles::NeighborsLoopParams searchInFluid( fluid->particles );

      const auto view_points = fluid->getParticles()->getPoints().getConstView();
      auto view_v = fluid->getVariables()->v.getView();
      const auto view_rho = fluid->getVariables()->v.getRhoView();
      const RealType searchRadius = fluid->getParticles()->getSearchRadius();
      const RealType m = modelParams.mass;
      const RealType h = modelParams.h;
      const RealType rho0 = modelParams.rho0;

      auto VelocityFilter = [=] __cuda_callable__ ( LocalIndexType i,
                                                    LocalIndexType j,
                                                    VectorType& r_i,
                                                    VectorType* v_i,
                                                    RealType* gamma_i ) mutable
      {
         const VectorType r_j = view_points[ j ];
         const VectorType r_ij = r_i - r_j;
         const RealType drs = l2Norm( r_ij );
         if (drs <= searchRadius )
         {
            const RealType rho_j = view_rho[ j ];
            const VectorType v_j = view_v[ j ];
            const RealType W = KernelFunction::W( drs, h );
            const RealType V_j = m / rho_j;

            *v_i += v_j* W * V_j;
            *gamma_i += W * V_j;
         }
      };

      auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
      {
         const VectorType r_i = view_points[ i ];
         VectorType v_i = 0.f;
         RealType gamma_i = 0.f;
         Particles::NeighborsLoop::exec( i, r_i, searchInFluid, VelocityFilter, &v_i,  &gamma_i );

         if( gamma_i > 0.01f )
            view_v[ i ] = v_i / gamma_i;
      };
      fluid->getParticles()->forAll( particleLoop );
   }

   template< typename Particles, typename FluidPointer, typename ModelParams >
   static void
   filterVelocityOverlaps( FluidPointer& fluid, ModelParams& modelParams )
   {
      typename Particles::NeighborsLoopParams searchInFluid( fluid->particles );

      const auto view_points = fluid->getParticles()->getPoints().getConstView();
      const auto view_rho = fluid->getVariables()->rho.getView();
      auto view_v = fluid->getVariables()->v.getView();
      //TEMP: Store velocity into acceleration due this overlap incident
      auto view_a = fluid->getVariables()->a.getView();
      auto view_gamma = fluid->getVariables()->gamma.getView();
      const RealType searchRadius = fluid->getParticles()->getSearchRadius();
      const RealType m = modelParams.mass;
      const RealType h = modelParams.h;
      const RealType rho0 = modelParams.rho0;

      auto VelocityFilter = [=] __cuda_callable__ ( LocalIndexType i,
                                                    LocalIndexType j,
                                                    VectorType& r_i,
                                                    VectorType* v_i,
                                                    RealType* gamma_i ) mutable
      {
         const VectorType r_j = view_points[ j ];
         const VectorType r_ij = r_i - r_j;
         const RealType drs = l2Norm( r_ij );
         if (drs <= searchRadius )
         {
            const RealType rho_j = view_rho[ j ];
            const VectorType v_j = view_v[ j ];
            const RealType W = KernelFunction::W( drs, h );
            const RealType V_j = m / rho_j;

            *v_i += v_j* W * V_j;
            *gamma_i += W * V_j;
         }
      };

      auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
      {
         const VectorType r_i = view_points[ i ];
         VectorType v_i = 0.f;
         RealType gamma_i = 0.f;
         Particles::NeighborsLoop::exec( i, r_i, searchInFluid, VelocityFilter, &v_i,  &gamma_i );

         view_a[ i ] = v_i;
         view_gamma[ i ] = gamma_i;
      };
      fluid->getParticles()->forAll( particleLoop );

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
               VectorType v_i = 0.f;
               RealType gamma_i = 0.f;
               Particles::NeighborsLoop::exec( p, r_i, searchInFluid, VelocityFilter, &v_i,  &gamma_i );

               view_a[ p ] += v_i;
               view_gamma[ p ] += gamma_i;
            };
            Algorithms::parallelFor< DeviceType >( 0, numberOfZoneParticles, periodicParticleLoop );
         }
      }

      // finalize interaction
      auto finalize = [ = ] __cuda_callable__( LocalIndexType i ) mutable
      {
         const RealType gamma_i = view_gamma[ i ];
         if( gamma_i > 0.01f )
            view_v[ i ] = view_a[ i ] / gamma_i;
      };
      TNL::Algorithms::parallelFor< DeviceType >(
            fluid->getFirstActiveParticle(), fluid->getLastActiveParticle() + 1, finalize );
   }

};


} // ViscousTerms
} // SPH
} // TNL

