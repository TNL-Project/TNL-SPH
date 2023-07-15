#pragma once

#include <TNL/Containers/Array.h>
#include <TNL/Containers/ArrayView.h>
#include <thrust/execution_policy.h>
#include <thrust/gather.h>

#include "Variables.h"
#include "../../SPHTraits.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHFluidConfig >
class IntegratorVariables
{
   public:
   using SPHFluidTraitsType = SPHFluidTraits< SPHFluidConfig >;

   using GlobalIndexType = typename SPHFluidTraitsType::GlobalIndexType;
   using RealType = typename SPHFluidTraitsType::RealType;

   using ScalarArrayType = typename SPHFluidTraitsType::ScalarArrayType;
   using VectorArrayType = typename SPHFluidTraitsType::VectorArrayType;

   using IndexArrayType = typename SPHFluidTraitsType::IndexArrayType;
   using IndexArrayTypePointer = typename Pointers::SharedPointer< IndexArrayType, typename SPHFluidConfig::DeviceType >;

   IntegratorVariables( GlobalIndexType size )
   : rho_old( size ), v_old( size ), rho_old_swap( size ), v_old_swap( size )
   {
      rho_old = 1000.; //TODO: Fix this.
      v_old = 0.;
   }

   void
   sortVariables( IndexArrayTypePointer& map, GlobalIndexType numberOfParticles, GlobalIndexType firstActiveParticle )
   {
      auto view_map = map->getView();

      auto view_rho_old = rho_old.getView();
      auto view_v_old = v_old.getView();

      auto view_rho_old_swap = rho_old_swap.getView();
      auto view_v_old_swap = v_old_swap.getView();

      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_rho_old.getArrayData() + firstActiveParticle, view_rho_old_swap.getArrayData() + firstActiveParticle );
      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_v_old.getArrayData() + firstActiveParticle, view_v_old_swap.getArrayData() + firstActiveParticle );

      rho_old.swap( rho_old_swap );
      v_old.swap( v_old_swap );
   }

   ScalarArrayType rho_old;
   VectorArrayType v_old;

   ScalarArrayType rho_old_swap;
   VectorArrayType v_old_swap;
};

template< typename ModelPointer, typename SPHFluidConfig, typename Variables = SPHFluidVariables< SPHFluidConfig > >
class VerletIntegrator
{
public:

   using SPHFluidTraitsType = SPHFluidTraits< SPHFluidConfig >;
   using DeviceType = typename SPHFluidConfig::DeviceType;

   using GlobalIndexType = typename SPHFluidTraitsType::GlobalIndexType;
   using RealType = typename SPHFluidTraitsType::RealType;
   using VectorType = typename SPHFluidTraitsType::VectorType;

   using IntegratorVariablesType = IntegratorVariables< SPHFluidConfig >;
   using IntegratorVariablesPointer = typename Pointers::SharedPointer< IntegratorVariablesType, DeviceType >;

   VerletIntegrator() = default;

   template< typename FluidPointer >
   void
   IntegrateVerlet( RealType dt, FluidPointer& fluid )
   {
      auto rho_view = fluid->variables->rho.getView();
      auto v_view = fluid->variables->v.getView();
      auto r_view = fluid->particles->getPoints().getView();

      auto rho_old_view = fluid->integratorVariables->rho_old.getView();
      auto v_old_view = fluid->integratorVariables->v_old.getView();

      const auto drho_view = fluid->variables->drho.getView();
      const auto a_view = fluid->variables->a.getView();

      RealType dtdt05 = 0.5 * dt * dt;
      RealType dt2 = 2 * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         r_view[ i ] += v_view[ i ] * dt + a_view[ i ] * dtdt05;
         v_old_view[ i ] += a_view[ i ] * dt2;
         rho_old_view[ i ] += drho_view[ i ] * dt2;
      };
      Algorithms::parallelFor< DeviceType >( fluid->getFirstActiveParticle(), fluid->getLastActiveParticle() + 1, init );

      fluid->variables->v.swap( fluid->integratorVariables->v_old );
      fluid->variables->rho.swap( fluid->integratorVariables->rho_old );
   }

   template< typename BoundaryPointer >
   void
   IntegrateVerletBoundary( RealType dt, BoundaryPointer& boundary )
   {
      auto rho_view = boundary->variables->rho.getView();
      auto rho_old_view = boundary->integratorVariables->rho_old.getView();

      const auto drho_view = boundary->variables->drho.getView();

      RealType dtdt05 = 0.5 * dt * dt;
      RealType dt2 = 2 * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         rho_old_view[ i ] += drho_view[ i ] * dt2;
      };
      Algorithms::parallelFor< DeviceType >( boundary->getFirstActiveParticle(), boundary->getLastActiveParticle() + 1, init );

      boundary->variables->rho.swap( boundary->integratorVariables->rho_old );
   }

   template< typename FluidPointer >
   void
   IntegrateEuler( RealType dt, FluidPointer& fluid )
   {
      auto rho_view = fluid->variables->rho.getView();
      auto v_view = fluid->variables->v.getView();
      auto r_view = fluid->particles->getPoints().getView();

      auto rho_old_view = fluid->integratorVariables->rho_old.getView();
      auto v_old_view = fluid->integratorVariables->v_old.getView();

      const auto drho_view = fluid->variables->drho.getView();
      const auto a_view = fluid->variables->a.getView();

      RealType dtdt05 = 0.5 * dt * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         r_view[ i ] += v_view[ i ] * dt + a_view[ i ] * dtdt05;
         v_old_view[ i ] = v_view[ i ];
         v_view[ i ] += a_view[ i ] * dt;
         rho_old_view[ i ] = rho_view[ i ];
         rho_view[ i ] += drho_view[ i ] * dt;
      };
      Algorithms::parallelFor< DeviceType >( fluid->getFirstActiveParticle(), fluid->getLastActiveParticle() + 1, init );
   }

   template< typename BoundaryPointer >
   void
   IntegrateEulerBoundary( RealType dt, BoundaryPointer& boundary )
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
      Algorithms::parallelFor< DeviceType >( boundary->getFirstActiveParticle(), boundary->getLastActiveParticle() + 1, init );
   }

   template< typename FluidPointer, typename BoundaryPointer, typename TimeStepping >
   void
   integratStepVerlet( FluidPointer& fluid, BoundaryPointer& boundary, TimeStepping& timeStepping )
   {
      if( timeStepping.getStep() % 20 == 0 ) {
         IntegrateEuler( timeStepping.getTimeStep(), fluid ); //TODO: Timer!
         //IntegrateEulerBoundary( timeStepping.getTimeStep(), boundary );
      }
      else {
         IntegrateVerlet( timeStepping.getTimeStep(), fluid );
         //IntegrateVerletBoundary( timeStepping.getTimeStep(), boundary );
      }
   }


};

} // SPH
} // ParticleSystem
} // TNL

