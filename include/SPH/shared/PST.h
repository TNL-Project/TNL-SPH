#pragma once

#include "../SPHTraits.h"

namespace TNL {
namespace SPH {

template< typename ParticlesType, typename SPHConfig >
class PST
{
public:

   template< typename FluidPointer, typename ModelParams >
   void
   shift( FluidPointer& fluid, ModelParams& modelParams )
   {
      auto searchInFluid = fluid->getParticles()->getSearchToken( fluid->getParticles() );

      const RealType searchRadius = fluid->getParticles()->getSearchRadius();
      const RealType h = modelParams.h;
      const RealType m = modelParams.mass;

      // PST "magical" constants
      const RealType const_A = 2.f;
      const RealType const_A_fsm = ParticlesType::spaceDimension;
      const RealType const_A_fst = ( ParticlesType::spaceDimension == 2 ) ? 1.5 : 2.75;

      auto view_points = fluid->getParticles()->getPoints().getConstView();
      const auto view_v = fluid->getVariables()->v.getConstView();
      const auto view_rho = fluid->getVariables()->rho.getConstView();

      auto getGradC = [=] __cuda_callable__ (
            IndexType i,
            IndexType j,
            VectorType& r_i
            RealType* gamma_i,
            VectorType* gradC_i,
            VectorType* div_r_i ) mutable
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
            *div_r_x += ( -1.f ) * ( r_xj, gradW ) * V_j;
         }
      };

      auto particleLoop = [=] __cuda_callable__ ( IndexType i ) mutable
      {
         const VectorType r_i = view_points[ i ] + shift;
         const VectorType v_i = view_v[ i ];
         RealType gamma_i = 0.f;
         VectorType gradC_i = 0.f;
         VectorType div_r_i = 0.f;

         ParticlesType::NeighborsLoop::exec( i, r_i, searchInFluid, getGradC, &gamma_i, &gradC_i, &div_r_i );

         const RealType const_A_fsc = ( div_r_i - const_A_fts )  / ( const_A_fsm - const_A_fst);
         const RealType delta_r_i_core= * const_A * h * l2Norm( v_i ) * gradC_i * dt;

         VectorType delta_r_i = 0;
         if ( ( div_r_i - const_A_fts  ) < 0 )
            delta_r_i =  const_A_fsc * delta_r_i;
         else if ( ( div_r_i - const_A_fts  ) = 0 )
            delta_r_i = delta_r_i_core;

         view_points[ i ] += delta_r_i;
      };
      fluid->getParticles()->forAll( particleLoop );
   }
};

} // PST
} // SPH
} // TNL

