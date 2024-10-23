#pragma once

#include <TNL/Containers/Array.h>
#include <TNL/Containers/ArrayView.h>
#include <TNL/Pointers/SharedPointer.h>

#include "../../../shared/thrustExecPolicySelector.h"
#include <thrust/sort.h>
#include <thrust/gather.h>

#include "../../../SPHTraits.h"

namespace TNL {
namespace SPH {
namespace IntegrationSchemes {

template< typename SPHConfig >
class MidpointIntegrationSchemeVariables
{
   public:
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using ScalarArrayType = typename SPHTraitsType::ScalarArrayType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;
   using IndexArrayType = typename SPHTraitsType::IndexArrayType;
   using IndexArrayTypePointer = typename Pointers::SharedPointer< IndexArrayType, typename SPHConfig::DeviceType >;

   MidpointIntegrationSchemeVariables() = default;

   void
   setSize( const GlobalIndexType& size )
   {
      dvdt_in.setSize( size );
      drhodt_in.setSize( size );
      r_in.setSize( size );
      v_in.setSize( size );
      rho_in.setSize( size );
      residua.setSize( size );

      dvdt_in_swap.setSize( size );
      drhodt_in_swap.setSize( size );
      r_in_swap.setSize( size );
      v_in_swap.setSize( size );
      rho_in_swap.setSize( size );
      residua_swap.setSize( size );

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

      //FIXME: disgusting out-place swap
      swapScalar.setSize( size );
      swapVector.setSize( size );
   }

   /*
   void
   sortVariables( IndexArrayTypePointer& map, GlobalIndexType numberOfParticles )
   {
      auto view_map = map->getView();

      using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< typename SPHConfig::DeviceType >;
      ThrustDeviceType thrustDevice;
      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            dvdt_in.getArrayData(), swapVector.getArrayData() );
      dvdt_in.swap( swapVector );
      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            drhodt_in.getArrayData(), swapScalar.getArrayData() );
      drhodt_in.swap( swapScalar );
      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            r_in.getArrayData(), swapVector.getArrayData() );
      r_in.swap( swapVector );
      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            v_in.getArrayData(), swapVector.getArrayData() );
      v_in.swap( swapVector );
      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            rho_in.getArrayData(), swapScalar.getArrayData() );
      rho_in.swap( swapScalar );
      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            residua.getArrayData(), swapScalar.getArrayData() );
      residua.swap( swapScalar );
   }
   */

   void
   sortVariables( IndexArrayTypePointer& map, GlobalIndexType numberOfParticles )
   {
      auto view_map = map->getView();

      using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< typename SPHConfig::DeviceType >;
      ThrustDeviceType thrustDevice;
      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            dvdt_in.getArrayData(), dvdt_in_swap.getArrayData() );
      dvdt_in.swap( dvdt_in_swap );

      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            drhodt_in.getArrayData(), drhodt_in_swap.getArrayData() );
      drhodt_in.swap( drhodt_in_swap );

      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            r_in.getArrayData(), r_in_swap.getArrayData() );
      r_in.swap( r_in_swap );

      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            v_in.getArrayData(), v_in_swap.getArrayData() );
      v_in.swap( v_in_swap );

      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            rho_in.getArrayData(), rho_in_swap.getArrayData() );
      rho_in.swap( rho_in_swap );

      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            residua.getArrayData(), residua_swap.getArrayData() );
      residua.swap( residua_swap );
   }

   VectorArrayType dvdt_in;
   ScalarArrayType drhodt_in;

   VectorArrayType r_in;
   VectorArrayType v_in;
   ScalarArrayType rho_in;

   ScalarArrayType residua;

   //FIXME
   VectorArrayType dvdt_in_swap;
   ScalarArrayType drhodt_in_swap;

   VectorArrayType r_in_swap;
   VectorArrayType v_in_swap;
   ScalarArrayType rho_in_swap;

   ScalarArrayType residua_swap;

   //FIXME: disgusting out-place swap
   ScalarArrayType swapScalar;
   VectorArrayType swapVector;
};

template< typename SPHConfig >
class MidpointScheme
{
public:

   using DeviceType = typename SPHConfig::DeviceType;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using PairIndexType = Containers::StaticVector< 2, GlobalIndexType >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;
   using IntegrationSchemeVariablesType = MidpointIntegrationSchemeVariables< SPHConfig >;

   template< typename FluidPointer >
   void
   predictor( RealType dt, FluidPointer& fluid )
   {
      fluid->integratorVariables->dvdt_in = fluid->variables->a;
      fluid->integratorVariables->drhodt_in = fluid->variables->drho;

      fluid->integratorVariables->r_in = fluid->getParticles()->getPoints();
      fluid->integratorVariables->v_in = fluid->variables->v;
      fluid->integratorVariables->rho_in = fluid->variables->rho;
   }

   template< typename FluidPointer >
   void
   midpointUpdateVariables( RealType dt, FluidPointer& fluid )
   {
      auto v_view = fluid->variables->v.getView();
      const auto v_in_view = fluid->integratorVariables->v_in.getConstView();
      const auto dvdt_view = fluid->variables->a.getConstView();
      auto rho_view = fluid->variables->rho.getView();
      const auto rho_in_view = fluid->integratorVariables->rho_in.getConstView();
      const auto drhodt_view = fluid->variables->drho.getConstView();

      const RealType dt05 = 0.5f * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         v_view[ i ] = v_in_view[ i ] + dt05 * dvdt_view[ i ];
         rho_view[ i ] = rho_in_view[ i ] + dt05 * drhodt_view[ i ];
      };
      fluid->particles->forAll( init );
   }

   template< typename FluidPointer >
   void
   midpointUpdatePositions( RealType dt, FluidPointer& fluid )
   {
      auto r_view = fluid->getParticles()->getPoints().getView();
      const auto r_in_view = fluid->integratorVariables->r_in.getConstView();
      const auto v_view = fluid->variables->v.getConstView();

      const RealType dt05 = 0.5f * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         r_view[ i ] = r_in_view[ i ] + 0.5f * v_view[ i ];
      };
      fluid->particles->forAll( init );
   }

   template< typename FluidPointer >
   void
   relax( FluidPointer& fluid, const RealType relaxMidpoint )
   {
      fluid->variables->a = relaxMidpoint * fluid->integratorVariables->dvdt_in + ( 1.f - relaxMidpoint ) * fluid->variables->a;
      fluid->variables->drho = relaxMidpoint * fluid->integratorVariables->drhodt_in + ( 1.f - relaxMidpoint ) * fluid->variables->drho;

      std::cout << "relaxed with relax midpoint: " << relaxMidpoint << std::endl;
   }

   template< typename FluidPointer, typename ModelParams >
   void
   midpointResiduals( FluidPointer& fluid, ModelParams& modelParams )
   {
      using EOS = typename ModelParams::EOS;

      const RealType m = modelParams.mass;
      typename EOS::ParamsType eosParams( modelParams );

      auto residua_view = fluid->integratorVariables->residua.getView();
      const auto v_view = fluid->variables->v.getConstView();
      const auto dvdt_view = fluid->variables->a.getConstView();
      const auto dvdt_in_view = fluid->integratorVariables->dvdt_in.getConstView();
      const auto rho_view = fluid->variables->rho.getConstView();
      const auto drhodt_view = fluid->variables->drho.getConstView();
      const auto drhodt_in_view = fluid->integratorVariables->drhodt_in.getConstView();

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         const RealType rho_i = rho_view[ i ];
         const RealType rho2_i = rho_i * rho_i;
         const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );

         const RealType res_dv_dt = m * std::abs( ( v_view[ i ], dvdt_view[ i ] - dvdt_in_view[ i ] ) );
         const RealType res_drho_dt = m * std::abs( ( p_i / rho2_i ) * ( drhodt_view[ i ] - drhodt_in_view[ i ] ) );
         residua_view[ i ] = res_dv_dt + res_drho_dt;
      };
      fluid->particles->forAll( init );
   }

   //TODO: This can be merged together with midpointResiduals function
   template< typename FluidPointer >
   const RealType
   getMaxResidua( FluidPointer& fluid )
   {
      return TNL::sum( fluid->integratorVariables->residua );
   }

   template< typename FluidPointer >
   void
   corrector( RealType dt, FluidPointer& fluid )
   {
      auto r_view = fluid->getParticles()->getPoints().getView();
      const auto r_in_view = fluid->integratorVariables->r_in.getConstView();
      auto v_view = fluid->variables->v.getView();
      const auto v_in_view = fluid->integratorVariables->v_in.getConstView();
      const auto dvdt_view = fluid->variables->a.getConstView();
      auto rho_view = fluid->variables->rho.getView();
      const auto rho_in_view = fluid->integratorVariables->rho_in.getConstView();
      const auto drhodt_view = fluid->variables->drho.getConstView();

      const RealType dtdt05 = 0.5 * dt * dt;
      const RealType dt2 = 2 * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         r_view[ i ] = r_in_view[ i ] + dt * v_in_view[ i ] + dtdt05 * dvdt_view[ i ];
         v_view[ i ] = v_in_view[ i ] + dt * dvdt_view[ i ];
         rho_view[ i ] = rho_in_view[ i ] + dt * drhodt_view[ i ];
      };
      fluid->particles->forAll( init );
   }
};

} // IntegrationSchemes
} // SPH
} // TNL

