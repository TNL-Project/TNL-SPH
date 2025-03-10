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
class RK45IntegrationSchemeVariables
{
   public:
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using ScalarArrayType = typename SPHTraitsType::ScalarArrayType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;
   using IndexArrayType = typename SPHTraitsType::IndexArrayType;
   using IndexArrayTypePointer = typename Pointers::SharedPointer< IndexArrayType, typename SPHConfig::DeviceType >;

   RK45IntegrationSchemeVariables() = default;

   void
   setSize( const GlobalIndexType& size )
   {
      dvdt_in.setSize( size );
      drhodt_in.setSize( size );
      v_staged.setSize( size );
      r_in.setSize( size );
      v_in.setSize( size );
      rho_in.setSize( size );

      dvdt_in_swap.setSize( size );
      drhodt_in_swap.setSize( size );
      v_staged_swap.setSize( size );
      r_in_swap.setSize( size );
      v_in_swap.setSize( size );
      rho_in_swap.setSize( size );

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
            v_staged.getArrayData(), v_staged_swap.getArrayData() );
      v_staged.swap( v_staged_swap );

      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            r_in.getArrayData(), r_in_swap.getArrayData() );
      r_in.swap( r_in_swap );

      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            v_in.getArrayData(), v_in_swap.getArrayData() );
      v_in.swap( v_in_swap );

      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            rho_in.getArrayData(), rho_in_swap.getArrayData() );
      rho_in.swap( rho_in_swap );

   }

   VectorArrayType dvdt_in;
   ScalarArrayType drhodt_in;
   VectorArrayType v_staged;

   VectorArrayType r_in;
   VectorArrayType v_in;
   ScalarArrayType rho_in;

   ScalarArrayType residua;

   //FIXME
   VectorArrayType dvdt_in_swap;
   ScalarArrayType drhodt_in_swap;
   VectorArrayType v_staged_swap;

   VectorArrayType r_in_swap;
   VectorArrayType v_in_swap;
   ScalarArrayType rho_in_swap;

   //FIXME: disgusting out-place swap
   ScalarArrayType swapScalar;
   VectorArrayType swapVector;
};

template< typename SPHConfig >
class RK45Scheme
{
public:

   using DeviceType = typename SPHConfig::DeviceType;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using PairIndexType = Containers::StaticVector< 2, GlobalIndexType >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;
   using IntegrationSchemeVariablesType = RK45IntegrationSchemeVariables< SPHConfig >;

   template< typename FluidPointer >
   void
   predictor( RealType dt, FluidPointer& fluid, const int predictorStep )
   {
      auto r_view = fluid->getParticles()->getPoints().getView();
      const auto r_in_view = fluid->getIntegratorVariables()->r_in.getConstView();
      auto v_view = fluid->getVariables()->v.getView();
      const auto v_in_view = fluid->getIntegratorVariables()->v_in.getConstView();
      const auto dvdt_view = fluid->getVariables()->a.getConstView();
      auto dvdt_view_in = fluid->getIntegratorVariables()->dvdt_in.getView();
      auto rho_view = fluid->getVariables()->rho.getView();
      const auto rho_in_view = fluid->getIntegratorVariables()->rho_in.getConstView();
      const auto drhodt_view = fluid->getVariables()->drho.getConstView();
      auto drhodt_view_in = fluid->getIntegratorVariables()->drhodt_in.getView();
      auto v_staged_view = fluid->getIntegratorVariables()->v_staged.getView();

      auto step1 = [=] __cuda_callable__ ( int i ) mutable
      {
         //w^(1) = w^(0) + F(w^(0)) * dt/2
         r_view[ i ] = r_in_view[ i ] + 0.5f * dt * v_view[ i ];
         v_view[ i ] = v_in_view[ i ] + 0.5f * dt * dvdt_view[ i ];
         rho_view[ i ] = rho_in_view[ i ] + 0.5f * dt * drhodt_view[ i ];
      };

      auto step2 = [=] __cuda_callable__ ( int i ) mutable
      {
         //w^(2) = w^(0) + F(w^(1)) * dt/2
         r_view[ i ] = r_in_view[ i ] + 0.5f * dt * v_view[ i ];
         v_view[ i ] = v_in_view[ i ] + 0.5f * dt * dvdt_view[ i ];
         rho_view[ i ] = rho_in_view[ i ] + 0.5f * dt * drhodt_view[ i ];

         //F(w^(0)) + 2 * F(w^(1))
         dvdt_view_in[ i ] += 2 * dvdt_view[ i ];
         drhodt_view_in[ i ] += 2 * drhodt_view[ i ];
         v_staged_view[ i ] += 2 * v_view[ i ];
      };

      auto step3 = [=] __cuda_callable__ ( int i ) mutable
      {
         //w^(3) = w^(0) + F(w^(1)) * dt/2 + dt * F(w^(2))
         r_view[ i ] = r_in_view[ i ] + dt * v_view[ i ];
         v_view[ i ] = v_in_view[ i ] + dt * dvdt_view[ i ];
         rho_view[ i ] = rho_in_view[ i ] + dt * drhodt_view[ i ];

         //F(w^(0)) + 2 * F(w^(1)) + 2 * F(w^(2))
         dvdt_view_in[ i ] += 2 * dvdt_view[ i ];
         drhodt_view_in[ i ] += 2 * drhodt_view[ i ];
         v_staged_view[ i ] += 2 * v_view[ i ];
      };

      if( predictorStep == 0 ){
         //w^(0) = w^n
         fluid->getIntegratorVariables()->dvdt_in = fluid->getVariables()->a;
         fluid->getIntegratorVariables()->drhodt_in = fluid->getVariables()->drho;
         fluid->getIntegratorVariables()->r_in = fluid->getParticles()->getPoints();
         fluid->getIntegratorVariables()->v_in = fluid->getVariables()->v;
         fluid->getIntegratorVariables()->rho_in = fluid->getVariables()->rho;

         fluid->getIntegratorVariables()->v_staged = fluid->getVariables()->v;
      }
      else if( predictorStep == 1 ){
         fluid->getParticles()->forAll( step1 );
      }
      else if( predictorStep == 2 ){
         fluid->getParticles()->forAll( step2 );
      }
      else if( predictorStep == 3 ){
         fluid->getParticles()->forAll( step3 );
      }

   }

   template< typename FluidPointer >
   void
   corrector( const RealType dt, FluidPointer& fluid )
   {
      auto r_view = fluid->getParticles()->getPoints().getView();
      const auto r_in_view = fluid->getIntegratorVariables()->r_in.getConstView();
      auto v_view = fluid->getVariables()->v.getView();
      const auto v_in_view = fluid->getIntegratorVariables()->v_in.getConstView();
      const auto dvdt_view = fluid->getVariables()->a.getConstView();
      auto dvdt_view_in = fluid->getIntegratorVariables()->dvdt_in.getView();
      auto rho_view = fluid->getVariables()->rho.getView();
      const auto rho_in_view = fluid->getIntegratorVariables()->rho_in.getConstView();
      const auto drhodt_view = fluid->getVariables()->drho.getConstView();
      auto drhodt_view_in = fluid->getIntegratorVariables()->drhodt_in.getView();
      const auto v_staged_view = fluid->getIntegratorVariables()->v_staged.getConstView();

      const RealType dtdt05 = 0.5f * dt * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         //w^(4) = w^(0) + (dt / 6) * [F(w^(0)) + 2 * F(w^(1)) + 2 * F(w^(2)) + F(w^(3))]
         r_view[ i ] = r_in_view[ i ] + ( dt / 6.f ) * ( v_staged_view[ i ] + v_view[ i ] );
         v_view[ i ] = v_in_view[ i ] + ( dt / 6.f ) * ( dvdt_view_in[ i ] + dvdt_view[ i ] );
         rho_view[ i ] = rho_in_view[ i ] + ( dt / 6.f ) * ( drhodt_view_in[ i ] + drhodt_view[ i ] );
      };
      fluid->getParticles()->forAll( init );
   }

};

} // IntegrationSchemes
} // SPH
} // TNL

