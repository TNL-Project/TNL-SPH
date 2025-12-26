#pragma once

#include "../SPHTraits.h"

namespace TNL {
namespace SPH {

template< typename ParticlesType, typename ModelConfig >
class PST
{
public:

   using SPHTraitsType = SPHFluidTraits< typename ModelConfig::SPHConfig >;
   using IndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;
   using KernelFunction = typename ModelConfig::KernelFunction;

   using ScalarArrayViewType = typename Containers::VectorView< RealType, TNL::Devices::Cuda >;
   using VectorTypeArrayViewType = typename Containers::VectorView< VectorType, TNL::Devices::Cuda >;

   template< typename FluidPointer, typename BoundaryPointer, typename ModelParams >
   static void
   shift( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams, const RealType dt )
   {
      auto searchInFluid = fluid->getParticles()->getSearchToken( fluid->getParticles() );
      auto searchInBound = boundary->getParticles()->getSearchToken( boundary->getParticles() );

      const RealType searchRadius = fluid->getParticles()->getSearchRadius();
      const RealType h = modelParams.h;
      const RealType m = modelParams.mass;

      // PST "magical" constants
      const RealType const_A = 2.f;
      const RealType const_A_fsm = ParticlesType::spaceDimension;
      //const RealType const_A_tfs = ( ParticlesType::spaceDimension == 2 ) ? 1.5 : 2.75;
      const RealType const_A_tfs = 0.f;

      auto view_points = fluid->getParticles()->getPoints().getView();
      const auto view_v = fluid->getVariables()->v.getConstView();
      auto view_rho = fluid->getVariables()->rho.getView();

      auto view_points_boundary = boundary->getParticles()->getPoints().getView();
      //auto view_rho_boundary = boundary->getVariables()->rho.getView();

      auto getGradC = [=] __cuda_callable__ (
            IndexType i,
            IndexType j,
            VectorType& r_i,
            RealType* gamma_i,
            RealType* div_r_i,
            VectorType* gradC_i,
            VectorTypeArrayViewType& points_view,
            ScalarArrayViewType& rho_view ) mutable
      {
         const VectorType r_j = view_points[ j ];
         const VectorType r_ij = r_i - r_j;
         const RealType drs = l2Norm( r_ij );
         if( drs <= searchRadius )
         {
            const RealType rho_j = view_rho[ j ];
            const RealType WV_j = KernelFunction::W( drs, h ) * m / rho_j;
            const VectorType gradWV_j = r_ij * KernelFunction::F( drs, h ) * m / rho_j;

            *gamma_i += WV_j;
            *gradC_i += gradWV_j;
            *div_r_i += ( -1.f ) * ( r_ij, gradWV_j );
         }
      };

      auto cancelNearBoundary = [=] __cuda_callable__ (
            IndexType i,
            IndexType j,
            VectorType& r_i,
            bool* near_boundary ) mutable
      {
         const VectorType r_j = view_points_boundary[ j ];
         const VectorType r_ij = r_i - r_j;
         const RealType drs = l2Norm( r_ij );
         if( drs <= searchRadius )
         {
            *near_boundary = true;
         }
      };

      auto particleLoop = [=] __cuda_callable__ ( IndexType i ) mutable
      {
         const VectorType r_i = view_points[ i ];
         const VectorType v_i = view_v[ i ];
         RealType gamma_i = 0.f;
         RealType div_r_i = 0.f;
         VectorType gradC_i = 0.f;
         bool near_boundary = false;

         ParticlesType::NeighborsLoop::exec(
               i, r_i, searchInFluid, getGradC, &gamma_i, &div_r_i, &gradC_i, view_points, view_rho );
         //ParticlesType::NeighborsLoop::exec(
         //      i, r_i, searchInBound, getGradC, &gamma_i, &div_r_i, &gradC_i, view_points_boundary, view_rho_boundary );
         ParticlesType::NeighborsLoop::exec( i, r_i, searchInBound, cancelNearBoundary, &near_boundary );

         const RealType const_A_fsc = ( div_r_i - const_A_tfs )  / ( const_A_fsm - const_A_tfs );
         const VectorType delta_r_i_core = ( -1.f ) * const_A * h * l2Norm( v_i ) * gradC_i * dt;

         VectorType delta_r_i = 0;
         if ( div_r_i > const_A_tfs )
            delta_r_i = const_A_fsc * delta_r_i_core;
         else
            delta_r_i = 0;

         if( near_boundary == true )
            delta_r_i = 0;

         view_points[ i ] += delta_r_i;

         //printf(" [COEF: %f]", const_A * h * dt );
         //printf("[shifttfs: %f, coefumagn: %f, umagnpre: %f, rs.w: %f, coeftfs: %f, umagn: %f ]", shifttfs, coefumagn, umagnpre, rs.w, coeftfs, umagn );

         //const RealType query = div_r_i - const_A_tfs;
         //if( query > 0 )
         //   printf( " [ %.12f, %.12f ] ",  delta_r_i_core[ 0 ], delta_r_i_core[ 1 ] );
      };
      fluid->getParticles()->forAll( particleLoop );
   }

};

} // SPH
} // TNL

