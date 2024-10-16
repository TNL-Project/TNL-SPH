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

class IntegrationSchemeVariables
{
   public:
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using ScalarArrayType = typename SPHTraitsType::ScalarArrayType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;
   using IndexArrayType = typename SPHTraitsType::IndexArrayType;
   using IndexArrayTypePointer = typename Pointers::SharedPointer< IndexArrayType, typename SPHConfig::DeviceType >;

   IntegrationSchemeVariables() = default;

   void
   setSize( const GlobalIndexType& size )
   {
      dvdt_in.setSize( size );
      drhodt_in.setSize( size );
      r_in.setSize( size );
      v_in.setSize( size );
      rho_in.setSize( size );
      residua.setSize( size );

      //FIXME: disgusting out-place swap
      swapScalar.setSize( size );
      swapVector.setSize( size );
   }

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

   VectorArrayType dvdt_in;
   ScalarArrayType drhodt_in;

   VectorArrayType r_in;
   VectorArrayType v_in;
   ScalarArrayType rho_in;

   ScalarArrayType residua;

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
   using IntegrationSchemeVariablesType = IntegrationSchemeVariables< SPHConfig >;


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
      const RealType dt05 = 0.5f * dt;
      fluid->variables->v = fluid->integratorVariables->v_in + dt05 * fluid->variables->a;
      fluid->variables->rho = fluid->integratorVariables->rho_in + dt05 * fluid->variables->drho;
   }

   template< typename FluidPointer >
   void
   midpointUpdatePositions( RealType dt, FluidPointer& fluid )
   {
      const RealType dt05 = 0.5f * dt;
      fluid->getParticles()->getPoints() = fluid->integratorVariables.r_in + dt05 * fluid->variables->v;
   }

   template< typename FluidPointer >
   void
   relax( RealType dt, FluidPointer& fluid, const RealType relaxMidpoint )
   {
      fluid->variables->a = relaxMidpoint * fluid->integratorVariables->dvdt_in + ( 1 - relaxMidpoint ) * fluid->variables->a;
      fluid->variables->drho = relaxMidpoint * fluid->integratorVariables->drhodt_in + ( 1 - relaxMidpoint ) * fluid->variables->drho;
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

         const RealType res_dv_dt = m * std::abs( ( v_view[ i ], dvdt_view - dvdt_in_view[ i ] ) )
         const RealType res_drho_dt = m * std::abs( p_i / rho2_i ) * ( drhodt_view[ i ] - drhodt_in_view[ i ] );
         residua_view[ i ] = res_dv_dt + res_drho_dt;
      };
      fluid->particles->forAll( init );
   }

   const RealType
   getMaxResidua()
   {
      return TNL::Max( residua );
   }

   template< typename FluidPointer >
   void
   corrector( RealType dt, FluidPointer& fluid )
   {
      const RealType dtdt05 = 0.5 * dt * dt;
      const RealType dt2 = 2 * dt;

      fluid->getParticles()->getPoints() = fluid->integratorVariables->r_in + dt * fluid->integratorVariables->v_in + dtdt05 * fluid->variables->a;
      fluid->variables->v = fluid->integratorVariables->v_in + dt * fluid->variables->a;
      fluid->variables->rho = fluid->integratorVariables->rho_in + dt * fluid->variables->drho;
   }
};

} // IntegrationSchemes
} // SPH
} // TNL

