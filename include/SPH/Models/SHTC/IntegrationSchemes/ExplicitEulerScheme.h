#pragma once

#include "../../../SPHTraits.h"

namespace TNL {
namespace SPH {
namespace IntegrationSchemes {

template< typename SPHConfig >
class ExplicitEulerSchemeVariables
{
   //TODO: Default types?
   public:
   ExplicitEulerSchemeVariables() = default;

   template< typename IndexType >
   void
   setSize( const IndexType& size )
   {}

   template< typename ParticlesPointer >
   void
   sortVariables( ParticlesPointer& particles )
   {}
};

template< typename SPHConfig >
class ExplicitEulerScheme
{
public:

   using DeviceType = typename SPHConfig::DeviceType;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using PairIndexType = Containers::StaticVector< 2, GlobalIndexType >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;
   using IntegrationSchemeVariablesType = ExplicitEulerSchemeVariables< SPHConfig >;

   // distinguish between different integration stages with unified interface
   enum class Stages
   {
      updateVariables,
      updateDensity,
      updateVelocity,
      updateDistortion,
      moveParticles
   };

   ExplicitEulerScheme() = default;

   template< typename FluidPointer >
   void
   updateVariables( RealType dt, FluidPointer& fluid )
   {
      auto view_rho = fluid->getVariables()->rho.getView();
      const auto view_drhodt = fluid->getVariables()->drhodt.getConstView();
      auto view_v = fluid->getVariables()->v.getView();
      const auto view_dvdt = fluid->getVariables()->dvdt.getConstView();
      auto view_A = fluid->getVariables()->A.getView();
      const auto view_dAdt = fluid->getVariables()->dAdt.getConstView();

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         view_rho[ i ] += view_drhodt[ i ] * dt;
         view_v[ i ] += view_dvdt[ i ] * dt;
         view_A[ i ] += view_dAdt[ i ] * dt;
      };
      fluid->getParticles()->forAll( init );
   }

   template< typename FluidPointer >
   void
   moveParticles( RealType dt, FluidPointer& fluid )
   {
      fluid->getPoints() += dt * fluid->getVariables()->v;
   }

   template< typename FluidPointer, typename BoundaryPointer, typename TimeStepping >
   void
   integrate( FluidPointer& fluid, BoundaryPointer& boundary, TimeStepping& timeStepping, Stages stage )
   {
      //not symplectic //TODO: Asserts
      //not symplectic if( stage == Stages::updateVariables )
      //not symplectic    updateVariables( timeStepping.getTimeStep(), fluid );
      //not symplectic else if( stage == Stages::moveParticles )
      //not symplectic    moveParticles( timeStepping.getTimeStep(), fluid );

      //TODO: Asserts
      const RealType dt = timeStepping.getTimeStep();

      if( stage == Stages::updateDensity ){
         fluid->getVariables()->rho += dt * fluid->getVariables()->drhodt;
      }
      else if( stage == Stages::updateVelocity ){
         fluid->getVariables()->v += dt * fluid->getVariables()->dvdt;
      }
      else if( stage == Stages::updateDistortion ){
         auto view_A = fluid->getVariables()->A.getView();
         const auto view_dAdt = fluid->getVariables()->dAdt.getConstView();

         auto init = [=] __cuda_callable__ ( int i ) mutable
         {
            view_A[ i ] += view_dAdt[ i ] * dt;
         };
         fluid->getParticles()->forAll( init );
      }
      else if( stage == Stages::moveParticles ){
         moveParticles( timeStepping.getTimeStep(), fluid );
      }
   }
};

} // IntegrationSchemes
} // SPH
} // TNL

