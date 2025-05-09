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

   template< typename IndexArrayTypePointer, typename IndexType >
   void
   sortVariables( IndexArrayTypePointer& map, const IndexType numberOfParticles )
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
      moveParticles
   };

   ExplicitEulerScheme() = default;

   template< typename FluidPointer >
   void
   updateVariables( RealType dt, FluidPointer& fluid )
   {
      fluid->getVariables()->rho += dt * fluid->getVariables()->drhodt;
      fluid->getVariables()->v += dt * fluid->getVariables()->dvdt;
      fluid->getVariables()->A += dt * fluid->getVariables()->dAdt;
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
      //TODO: Asserts
      if( stage == Stages::updateVariables )
         updateVariables( timeStepping.getTimeStep(), fluid );
      else if( stage == Stages::moveParticles )
         moveParticles( timeStepping.getTimeStep(), fluid );
   }

};

} // IntegrationSchemes
} // SPH
} // TNL

