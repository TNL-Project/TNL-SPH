#pragma once

#include <TNL/Containers/Array.h>
#include <TNL/Containers/ArrayView.h>
#include <TNL/Pointers/SharedPointer.h>

#include <TNL/Particles/details/thrustExecPolicySelector.h>
#include <thrust/sort.h>
#include <thrust/gather.h>

#include "../../../SPHTraits.h"

namespace TNL {
namespace SPH {
namespace IntegrationSchemes {

template< typename SPHConfig >
class MidpointWithEnergySecantIntegrationSchemeVariables
{
   public:
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using ScalarArrayType = typename SPHTraitsType::ScalarArrayType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;
   using IndexArrayType = typename SPHTraitsType::IndexArrayType;
   using IndexArrayTypePointer = typename Pointers::SharedPointer< IndexArrayType, typename SPHConfig::DeviceType >;

   MidpointWithEnergySecantIntegrationSchemeVariables() = default;

   void
   setSize( const GlobalIndexType& size )
   {
      dvdt_in.setSize( size );
      drhodt_in.setSize( size );
      r_in.setSize( size );
      v_in.setSize( size );
      rho_in.setSize( size );
      residua.setSize( size );
      eps_k.setSize( size );
      eps_c.setSize( size );
      alpha_k.setSize( size );
      alpha_c.setSize( size );

      dvdt_in_swap.setSize( size );
      drhodt_in_swap.setSize( size );
      r_in_swap.setSize( size );
      v_in_swap.setSize( size );
      rho_in_swap.setSize( size );
      residua_swap.setSize( size );
      eps_k_swap.setSize( size );
      eps_c_swap.setSize( size );

      //FIXME:
      dvdt_in = 0.f;
      dvdt_in_swap = 0.f;
      drhodt_in = 0.f;
      drhodt_in_swap = 0.f;

      r_in = 0.f;
      r_in_swap = 0.f;

      v_in = 0.f;
      v_in_swap = 0.f;

      rho_in = 1000.f;
      rho_in_swap = 1000.f;

      residua = 0.f;
      residua_swap = 0.f;

      //Energy swap
      eps_k = 0.f;
      eps_k_swap = 0.f;
      eps_c = 0.f;
      eps_c_swap = 0.f;

      //FIXME: disgusting out-place swap
      swapScalar.setSize( size );
      swapVector.setSize( size );

      //FIXME: disgusting out-place swap
      rho_old.setSize( size );
      v_old.setSize( size );
   }

   template< typename ParticlesPointer >
   void
   sortVariables( ParticlesPointer& particles )
   {
      particles->reorderArray( dvdt_in, dvdt_in_swap );
      particles->reorderArray( drhodt_in, drhodt_in_swap );
      particles->reorderArray( r_in, r_in_swap );
      particles->reorderArray( v_in, v_in_swap );
      particles->reorderArray( rho_in, rho_in_swap );
      particles->reorderArray( residua, residua_swap );
      particles->reorderArray( eps_k, eps_k_swap );
      particles->reorderArray( eps_c, eps_c_swap );
   }

   // Midpoint integration scheme fields
   VectorArrayType dvdt_in;
   ScalarArrayType drhodt_in;

   VectorArrayType r_in;
   VectorArrayType v_in;
   ScalarArrayType rho_in;

   // backup data from previous iteration step
   VectorArrayType r_prev;
   VectorArrayType v_prev;
   ScalarArrayType rho_prev;

   ScalarArrayType residua;
   ScalarArrayType eps_k;
   ScalarArrayType eps_c;

   // relaxation factors
   ScalarArrayType alpha_k;
   ScalarArrayType alpha_c;

   // swap variables
   //FIXME
   VectorArrayType dvdt_in_swap;
   ScalarArrayType drhodt_in_swap;

   VectorArrayType r_in_swap;
   VectorArrayType v_in_swap;
   ScalarArrayType rho_in_swap;

   ScalarArrayType residua_swap;
   ScalarArrayType eps_k_swap;
   ScalarArrayType eps_c_swap;

   //FIXME: disgusting out-place swap
   ScalarArrayType swapScalar;
   VectorArrayType swapVector;

   //FIXME: ugly temp work around to deal with inlet buffer feading previous time steps
   ScalarArrayType rho_old;
   VectorArrayType v_old;
};

template< typename SPHConfig >
class MidpointSchemeWithEnergySecant
{
public:

   using DeviceType = typename SPHConfig::DeviceType;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using PairIndexType = Containers::StaticVector< 2, GlobalIndexType >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;
   using IntegrationSchemeVariablesType = MidpointWithEnergySecantIntegrationSchemeVariables< SPHConfig >;

   template< typename FluidPointer >
   void
   predictor( RealType dt, FluidPointer& fluid )
   {
      fluid->getIntegratorVariables()->dvdt_in = fluid->getVariables()->a;
      fluid->getIntegratorVariables()->drhodt_in = fluid->getVariables()->drho;

      fluid->getIntegratorVariables()->r_in = fluid->getParticles()->getPoints();
      fluid->getIntegratorVariables()->v_in = fluid->getVariables()->v;
      fluid->getIntegratorVariables()->rho_in = fluid->getVariables()->rho;
   }

   template< typename FluidPointer >
   void
   midpointUpdateVariables( RealType dt, FluidPointer& fluid )
   {
      // backup derivatives
      fluid->getIntegratorVariables()->drhodt_in = fluid->getVariables()->drho;
      fluid->getIntegratorVariables()->dvdt_in = fluid->getVariables()->a;

      auto v_view = fluid->getVariables()->v.getView();
      const auto v_in_view = fluid->getIntegratorVariables()->v_in.getConstView();
      const auto dvdt_view = fluid->getVariables()->a.getConstView();
      auto rho_view = fluid->getVariables()->rho.getView();
      const auto rho_in_view = fluid->getIntegratorVariables()->rho_in.getConstView();
      const auto drhodt_view = fluid->getVariables()->drho.getConstView();

      const RealType dt05 = 0.5f * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         v_view[ i ] = v_in_view[ i ] + dt05 * dvdt_view[ i ];
         rho_view[ i ] = rho_in_view[ i ] + dt05 * drhodt_view[ i ];
      };
      fluid->getParticles()->forAll( init );
   }

   template< typename FluidPointer >
   void
   midpointUpdatePositions( RealType dt, FluidPointer& fluid )
   {
      auto r_view = fluid->getParticles()->getPoints().getView();
      const auto r_in_view = fluid->getIntegratorVariables()->r_in.getConstView();
      const auto v_view = fluid->getVariables()->v.getConstView();

      const RealType dt05 = 0.5f * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         r_view[ i ] = r_in_view[ i ] + dt05 * v_view[ i ];
      };
      fluid->getParticles()->forAll( init );
   }

   template< typename FluidPointer, typename ModelParams >
   void
   relax( FluidPointer& fluid, ModelParams& modelParams )
   {
      RealType relaxMidpoint;
      if( midpointIteration == 0 )
         relaxMidpoint = modelParams.midpointRelaxCoef_0;
      else if( midpointIteration == modelParams.midpointMaxInterations )
         relaxMidpoint = 0.f;
      else
         relaxMidpoint = midpointRelaxCoef;

      fluid->getVariables()->a = relaxMidpoint * fluid->getIntegratorVariables()->dvdt_in + ( 1.f - relaxMidpoint ) * fluid->getVariables()->a;
      fluid->getVariables()->drho = relaxMidpoint * fluid->getIntegratorVariables()->drhodt_in + ( 1.f - relaxMidpoint ) * fluid->getVariables()->drho;
   }

   template< typename FluidPointer, typename ModelParams >
   RealType
   midpointResiduals( FluidPointer& fluid, ModelParams& modelParams )
   {
      using EOS = typename ModelParams::EOS;

      const RealType m = modelParams.mass;
      typename EOS::ParamsType eosParams( modelParams );

      auto residua_view = fluid->getIntegratorVariables()->residua.getView();
      auto eps_k_view = fluid->getIntegratorVariables()->eps_k.getView();
      auto eps_c_view = fluid->getIntegratorVariables()->eps_c.getView();
      auto alpha_k_view = fluid->getIntegratorVariables()->alpha_k.getView();
      auto alpha_c_view = fluid->getIntegratorVariables()->alpha_c.getView();

      const auto dvdt_view = fluid->getVariables()->a.getConstView();
      const auto dvdt_in_view = fluid->getIntegratorVariables()->dvdt_in.getConstView();
      const auto drhodt_view = fluid->getVariables()->drho.getConstView();
      const auto drhodt_in_view = fluid->getIntegratorVariables()->drhodt_in.getConstView();

      const auto v_view = fluid->getVariables()->v.getConstView();
      const auto rho_view = fluid->getVariables()->rho.getConstView();

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         const VectorType v_i = v_view[ i ];
         const RealType rho_i = rho_view[ i ];
         const RealType rho2_i = rho_i * rho_i;
         const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );

         const VectorType v_prev_i = v_prev_view[ i ];
         const RealType rho_prev_i = rho_prev_view[ i ];
         const RealType rho2_prev_i = rho_prev_i * rho_prev_i;
         const RealType p_prev_i = EOS::DensityToPressure( rho_prev_i, eosParams );

         const VectorType dvdt_i = dvdt_view[ i ];
         const RealType drhodt_i = drhodt_view[ i ];

         //const RealType res_dv_dt = m * std::abs( ( v_view[ i ], dvdt_view[ i ] - dvdt_in_view[ i ] ) );
         //const RealType res_drho_dt = m * std::abs( ( p_i / rho2_i ) * ( drhodt_view[ i ] - drhodt_in_view[ i ] ) );

         const RealType p_k_prev = m * std::abs( ( v_prev_i, dvdt_i ) );
         const RealType p_k_new = m * std::abs( ( v_i, dvdt_i ) );

         const RealType p_c_prev = m * std::abs( ( p_prev_i / rho2_prev_i ) * drhodt_i );
         const RealType p_c_new = m * std::abs( ( p_i / rho2_i ) * drhodt_i );


         // compute eps

         // compute alpha

         // update eps

         //residua_view[ i ] = res_dv_dt + res_drho_dt;
      };
      fluid->getParticles()->forAll( init );

      return TNL::sum( fluid->getIntegratorVariables()->residua );
   }

   template< typename FluidPointer >
   void
   corrector( RealType dt, FluidPointer& fluid )
   {
      auto r_view = fluid->getParticles()->getPoints().getView();
      const auto r_in_view = fluid->getIntegratorVariables()->r_in.getConstView();
      auto v_view = fluid->getVariables()->v.getView();
      const auto v_in_view = fluid->getIntegratorVariables()->v_in.getConstView();
      const auto dvdt_view = fluid->getVariables()->a.getConstView();
      auto rho_view = fluid->getVariables()->rho.getView();
      const auto rho_in_view = fluid->getIntegratorVariables()->rho_in.getConstView();
      const auto drhodt_view = fluid->getVariables()->drho.getConstView();

      const RealType dtdt05 = 0.5 * dt * dt;
      const RealType dt2 = 2 * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         r_view[ i ] = r_in_view[ i ] + dt * v_in_view[ i ] + dtdt05 * dvdt_view[ i ];
         v_view[ i ] = v_in_view[ i ] + dt * dvdt_view[ i ];
         rho_view[ i ] = rho_in_view[ i ] + dt * drhodt_view[ i ];
      };
      fluid->getParticles()->forAll( init );
   }

   template< typename ModelParams >
   bool
   runMidpointSubiteration( ModelParams& modelParams )
   {
      if( midpointIteration == 0 ){
         midpointIteration = 0;
         residual = 0.f;
         residualPrevious = 0.f;
         midpointRelaxCoef = modelParams.midpointRelaxCoef;
      }
      // backup iteration due to midpoint output log;
      midpointIterationBackup = midpointIteration + 1;


      if( midpointIteration < modelParams.midpointMaxInterations ) {
         midpointIteration++;
         return true;
      }
      else {
         midpointIteration = 0;
         return false;
      }


   }

   template< typename FluidPointer, typename ModelParams >
   void
   computeResiduals( FluidPointer& fluid, ModelParams& modelParams )
   {
      // compute residuals and control
      residual = midpointResiduals( fluid, modelParams );
      if( residual < modelParams.midpointResidualTolerance )
         this->midpointIteration = modelParams.midpointMaxInterations;
   }

   template< typename ModelParams >
   void
   updateRelaxationFactor( ModelParams& modelParams )
   {
      // control residuals decay (NOTE: -1 since midpointIteration are icreased at the start of the loop)
      if( ( midpointIteration - 1 )  > 0 )
         if( residual / residualPrevious > modelParams.midpointResidualMinimalDecay )
            midpointRelaxCoef = modelParams.midpointRelaxCoefIncrement
                              + ( 1.0 - modelParams.midpointRelaxCoefIncrement ) * midpointRelaxCoef;

      // backup the residuals
      residualPrevious = residual;
   }


   // Midpoint integration scheme  variables
   int midpointIteration = 0;
   int midpointIterationBackup = 0; // temp iteraction backup just due to outpu
   RealType residual = 0.f;
   RealType residualPrevious = 0.f;
   RealType midpointRelaxCoef = 0.f;
};

} // IntegrationSchemes
} // SPH
} // TNL

