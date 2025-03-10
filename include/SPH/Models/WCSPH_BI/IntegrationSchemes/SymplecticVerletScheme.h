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
class SymplecticVerletSchemeVariables
{
   public:
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using ScalarArrayType = typename SPHTraitsType::ScalarArrayType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;
   using IndexArrayType = typename SPHTraitsType::IndexArrayType;
   using IndexArrayTypePointer = typename Pointers::SharedPointer< IndexArrayType, typename SPHConfig::DeviceType >;

   SymplecticVerletSchemeVariables() = default;

   SymplecticVerletSchemeVariables( GlobalIndexType size )
   : rho_old( size ), r_old( size ), v_old( size ), rho_old_swap( size ), r_old_swap( size ), v_old_swap( size ) {}

   void
   setSize( const GlobalIndexType& size )
   {
      rho_old.setSize( size );
      rho_old = 0.f;
      r_old.setSize( size );
      r_old = 0.f;
      v_old.setSize( size );
      v_old = 0.f;
      rho_old_swap.setSize( size );
      rho_old_swap = 0.f;
      r_old_swap.setSize( size );
      r_old_swap = 0.f;
      v_old_swap.setSize( size );
      v_old_swap = 0.f;
   }

   void
   sortVariables( IndexArrayTypePointer& map, GlobalIndexType numberOfParticles )
   {
      auto view_map = map->getView();

      auto view_rho_old = rho_old.getView();
      auto view_r_old = r_old.getView();
      auto view_v_old = v_old.getView();
      auto view_rho_old_swap = rho_old_swap.getView();
      auto view_r_old_swap = r_old_swap.getView();
      auto view_v_old_swap = v_old_swap.getView();

      using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< typename SPHConfig::DeviceType >;
      ThrustDeviceType thrustDevice;
      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_rho_old.getArrayData(), view_rho_old_swap.getArrayData() );
      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_r_old.getArrayData(), view_r_old_swap.getArrayData() );
      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_v_old.getArrayData(), view_v_old_swap.getArrayData() );

      r_old.swap( r_old_swap );
      rho_old.swap( rho_old_swap );
      v_old.swap( v_old_swap );
   }

   void
   sortVariables( IndexArrayTypePointer& map, GlobalIndexType numberOfParticles, GlobalIndexType firstActiveParticle )
   {
      auto view_map = map->getView();

      auto view_rho_old = rho_old.getView();
      auto view_r_old = r_old.getView();
      auto view_v_old = v_old.getView();
      auto view_rho_old_swap = rho_old_swap.getView();
      auto view_r_old_swap = r_old_swap.getView();
      auto view_v_old_swap = v_old_swap.getView();

      using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< typename SPHConfig::DeviceType >;
      ThrustDeviceType thrustDevice;
      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_rho_old.getArrayData() + firstActiveParticle, view_rho_old_swap.getArrayData() + firstActiveParticle );
      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_r_old.getArrayData() + firstActiveParticle, view_r_old_swap.getArrayData() + firstActiveParticle );
      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_v_old.getArrayData() + firstActiveParticle, view_v_old_swap.getArrayData() + firstActiveParticle );

      r_old.swap( r_old_swap );
      rho_old.swap( rho_old_swap );
      v_old.swap( v_old_swap );
   }

   ScalarArrayType rho_old;
   VectorArrayType r_old;
   VectorArrayType v_old;

   ScalarArrayType rho_old_swap;
   VectorArrayType r_old_swap;
   VectorArrayType v_old_swap;
};

template< typename SPHConfig >
class SymplecticVerletScheme
{
public:

   using DeviceType = typename SPHConfig::DeviceType;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using PairIndexType = Containers::StaticVector< 2, GlobalIndexType >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;
   using IntegrationSchemeVariablesType = SymplecticVerletSchemeVariables< SPHConfig >;

   SymplecticVerletScheme() = default;

   /**
    * Print model identifier.
    */
   static std::string
   getSchemeName()
   {
      return "TNL::SPH::IntegrationSchemes::SymplecticVerletScheme";
   }


   template< typename FluidPointer >
   void
   predictorStep( RealType dt, FluidPointer& fluid )
   {
      auto v_view = fluid->getVariables()->v.getView();
      auto r_view = fluid->getParticles()->getPoints().getView();
      auto r_old_view = fluid->getIntegratorVariables()->r_old.getView();
      auto rho_view = fluid->getVariables()->rho.getView();
      auto rho_old_view = fluid->getIntegratorVariables()->rho_old.getView();
      auto v_old_view = fluid->getIntegratorVariables()->v_old.getView();
      const auto drho_view = fluid->getVariables()->drho.getConstView();
      const auto a_view = fluid->getVariables()->a.getConstView();

      const RealType dt05 = 0.5f * dt;

      //TODO: To resolve notes about fucked swap require initialize integrator variables with
      //      the initial condition fields.
      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         //TODO: Workaround due to slow swap, otherwise r_old_viw[ i ] += v_view[ i ] dt05 and r.swap( r_old )
         //r_old_viw[ i ] += v_view[ i ] * dt05

         r_old_view[ i ] = r_view[ i ];
         r_view[ i ] += v_view[ i ] * dt05;

         v_old_view[ i ] = v_view[ i ];
         v_view[ i ] += a_view[ i ] * dt05;

         rho_old_view[ i ] = rho_view[ i ];
         rho_view[ i ] += drho_view[ i ] * dt05;
      };
      fluid->getParticles()->forAll( init );

      //TODO: For some reason, fluid swap is insane slow, due to this, the integrator is build with a workaround
      //fluid->getVariables()->v.swap( fluid->getIntegratorVariables()->v_old );
      //fluid->getVariables()->rho.swap( fluid->getIntegratorVariables()->rho_old );
      //fluid->getParticles()->getPoints().swap( fluid->getIntegratorVariables()->r_old );
   }

   template< typename BoundaryPointer >
   void
   predictorStepBoundary( RealType dt, BoundaryPointer& boundary )
   {
      auto rho_old_view = boundary->getIntegratorVariables()->rho_old.getView();
      const auto drho_view = boundary->getVariables()->drho.getConstView();

      const RealType dt05 = 0.5f * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         rho_old_view[ i ] += drho_view[ i ] * dt05;
      };
      boundary->getParticles()->forAll( init );

      boundary->getVariables()->rho.swap( boundary->getIntegratorVariables()->rho_old );
   }

   template< typename FluidPointer >
   void
   correctorStep( RealType dt, FluidPointer& fluid )
   {
      auto r_view = fluid->getParticles()->getPoints().getView();
      auto r_old_view = fluid->getIntegratorVariables()->r_old.getView();
      auto v_view = fluid->getVariables()->v.getView();
      const auto v_old_view = fluid->getIntegratorVariables()->v_old.getConstView();
      auto rho_view = fluid->getVariables()->rho.getView();
      const auto rho_old_view = fluid->getIntegratorVariables()->rho_old.getConstView();
      const auto drho_view = fluid->getVariables()->drho.getConstView();
      const auto a_view = fluid->getVariables()->a.getConstView();

      const RealType dt05 = 0.5f * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         v_view[ i ] = v_old_view[ i ] + a_view[ i ] * dt;
         r_view[ i ] = r_old_view[ i ] + ( v_view[ i ] + v_old_view[ i ] ) * dt05;
         const RealType epsilon =  ( -1.f ) * ( drho_view[ i ] / rho_view[ i ] ) * dt;
         rho_view[ i ] = rho_old_view[ i ] * ( ( 2.f - epsilon ) / ( 2.f + epsilon ) );
      };
      fluid->getParticles()->forAll( init );
   }

   template< typename BoundaryPointer >
   void
   correctorStepBoundary( RealType dt, BoundaryPointer& boundary )
   {
      auto rho_view = boundary->getVariables()->rho.getView();
      const auto rho_old_view = boundary->getIntegratorVariables()->rho_old.getConstView();
      const auto drho_view = boundary->getVariables()->drho.getConstView();

      const RealType dt05 = 0.5 * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         const RealType epsilon =  - ( drho_view[ i ] / rho_view[ i ] ) * dt;
         rho_view[ i ] = rho_old_view[ i ] * ( ( 2.f - epsilon ) / ( 2.f + epsilon ) );
      };
      boundary->getParticles()->forAll( init );
   }


   template< typename BoundaryPointer >
   void
   correctBoundaryDensity( BoundaryPointer& boundary )
   {
      auto rho_view = boundary->getVariables()->rho.getView();

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         if( rho_view[ i ] < 1000.f )
            rho_view[ i ] = 1000.f; //TODO: Use referential density or move this to different part
      };
      boundary->getParticles()->forAll( init );
   }

   template< typename FluidPointer, typename BoundaryPointer, typename TimeStepping >
   void
   integratePredictorStep( FluidPointer& fluid, BoundaryPointer& boundary, TimeStepping& timeStepping )
   {
      predictorStep( timeStepping.getTimeStep(), fluid ); //TODO: Timer!
      //correctBoundaryDensity( boundary );
   }

   template< typename FluidPointer, typename BoundaryPointer, typename TimeStepping >
   void
   integrateCorrectorStep( FluidPointer& fluid, BoundaryPointer& boundary, TimeStepping& timeStepping )
   {
      correctorStep( timeStepping.getTimeStep(), fluid ); //TODO: Timer!
      //correctBoundaryDensity( boundary );
   }
};

} // IntegrationSchemes
} // SPH
} // TNL

