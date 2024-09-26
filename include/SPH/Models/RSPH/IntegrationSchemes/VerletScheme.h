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

   IntegrationSchemeVariables( GlobalIndexType size )
   : rho_old( size ), v_old( size ), rho_old_swap( size ), v_old_swap( size ) {}

   void
   setSize( const GlobalIndexType& size )
   {
      rho_old.setSize( size );
      v_old.setSize( size );
      rho_old_swap.setSize( size );
      v_old_swap.setSize( size );
   }

   void
   sortVariables( IndexArrayTypePointer& map, GlobalIndexType numberOfParticles )
   {
      auto view_map = map->getView();

      auto view_rho_old = rho_old.getView();
      auto view_v_old = v_old.getView();
      auto view_rho_old_swap = rho_old_swap.getView();
      auto view_v_old_swap = v_old_swap.getView();

      using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< typename SPHConfig::DeviceType >;
      ThrustDeviceType thrustDevice;
      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_rho_old.getArrayData(), view_rho_old_swap.getArrayData() );
      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_v_old.getArrayData(), view_v_old_swap.getArrayData() );

      rho_old.swap( rho_old_swap );
      v_old.swap( v_old_swap );
   }

   void
   sortVariables( IndexArrayTypePointer& map, GlobalIndexType numberOfParticles, GlobalIndexType firstActiveParticle )
   {
      auto view_map = map->getView();

      auto view_rho_old = rho_old.getView();
      auto view_v_old = v_old.getView();
      auto view_rho_old_swap = rho_old_swap.getView();
      auto view_v_old_swap = v_old_swap.getView();

      using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< typename SPHConfig::DeviceType >;
      ThrustDeviceType thrustDevice;
      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_rho_old.getArrayData() + firstActiveParticle, view_rho_old_swap.getArrayData() + firstActiveParticle );
      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_v_old.getArrayData() + firstActiveParticle, view_v_old_swap.getArrayData() + firstActiveParticle );

      rho_old.swap( rho_old_swap );
      v_old.swap( v_old_swap );
   }

   ScalarArrayType rho_old;
   VectorArrayType v_old;

   ScalarArrayType rho_old_swap;
   VectorArrayType v_old_swap;
};

template< typename SPHConfig >
class VerletScheme
{
public:

   using DeviceType = typename SPHConfig::DeviceType;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using PairIndexType = Containers::StaticVector< 2, GlobalIndexType >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;
   using IntegrationSchemeVariablesType = IntegrationSchemeVariables< SPHConfig >;

   VerletScheme() = default;

   template< typename FluidPointer >
   void
   integrateVerlet( RealType dt, FluidPointer& fluid )
   {
      auto v_view = fluid->variables->v.getView();
      auto r_view = fluid->particles->getPoints().getView();
      auto rho_old_view = fluid->integratorVariables->rho_old.getView();
      auto v_old_view = fluid->integratorVariables->v_old.getView();
      const auto drho_view = fluid->variables->drho.getConstView();
      const auto a_view = fluid->variables->a.getConstView();

      const RealType dtdt05 = 0.5 * dt * dt;
      const RealType dt2 = 2 * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         r_view[ i ] += v_view[ i ] * dt + a_view[ i ] * dtdt05;
         v_old_view[ i ] += a_view[ i ] * dt2;
         rho_old_view[ i ] += drho_view[ i ] * dt2;
      };
      fluid->particles->forAll( init );

      fluid->variables->v.swap( fluid->integratorVariables->v_old );
      fluid->variables->rho.swap( fluid->integratorVariables->rho_old );
   }

   template< typename BoundaryPointer >
   void
   integrateVerletBoundary( RealType dt, BoundaryPointer& boundary )
   {
      auto rho_old_view = boundary->integratorVariables->rho_old.getView();
      const auto drho_view = boundary->variables->drho.getConstView();

      const RealType dtdt05 = 0.5 * dt * dt;
      const RealType dt2 = 2 * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         rho_old_view[ i ] += drho_view[ i ] * dt2;
      };
      boundary->particles->forAll( init );

      boundary->variables->rho.swap( boundary->integratorVariables->rho_old );
   }

   template< typename FluidPointer >
   void
   integrateEuler( RealType dt, FluidPointer& fluid )
   {
      auto rho_view = fluid->variables->rho.getView();
      auto v_view = fluid->variables->v.getView();
      auto r_view = fluid->particles->getPoints().getView();
      auto rho_old_view = fluid->integratorVariables->rho_old.getView();
      auto v_old_view = fluid->integratorVariables->v_old.getView();
      const auto drho_view = fluid->variables->drho.getConstView();
      const auto a_view = fluid->variables->a.getConstView();

      const RealType dtdt05 = 0.5 * dt * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         r_view[ i ] += v_view[ i ] * dt + a_view[ i ] * dtdt05;
         v_old_view[ i ] = v_view[ i ];
         v_view[ i ] += a_view[ i ] * dt;
         rho_old_view[ i ] = rho_view[ i ];
         rho_view[ i ] += drho_view[ i ] * dt;
      };
      fluid->particles->forAll( init );
   }

   template< typename BoundaryPointer >
   void
   integrateEulerBoundary( RealType dt, BoundaryPointer& boundary )
   {
      auto rho_view = boundary->variables->rho.getView();
      auto rho_old_view = boundary->integratorVariables->rho_old.getView();

      const auto drho_view = boundary->variables->drho.getView();

      RealType dtdt05 = 0.5 * dt * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         rho_old_view[ i ] = rho_view[ i ];
         rho_view[ i ] += drho_view[ i ] * dt;
      };
      boundary->particles->forAll( init );
   }

   template< typename BoundaryPointer >
   void
   correctBoundaryDensity( BoundaryPointer& boundary )
   {
      auto rho_view = boundary->variables->rho.getView();

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         if( rho_view[ i ] < 1000.f )
            rho_view[ i ] = 1000.f;
      };
      boundary->particles->forAll( init );
   }

   template< typename FluidPointer, typename BoundaryPointer, typename TimeStepping >
   void
   integratStepVerlet( FluidPointer& fluid, BoundaryPointer& boundary, TimeStepping& timeStepping )
   {
      if( timeStepping.getStep() % 20 == 0 ) {
         integrateEuler( timeStepping.getTimeStep(), fluid ); //TODO: Timer!
      }
      else {
         integrateVerlet( timeStepping.getTimeStep(), fluid );
      }
      //correctBoundaryDensity( boundary );
   }

};

} // IntegrationSchemes
} // SPH
} // TNL

