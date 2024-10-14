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

   VectorArrayType dvdt_in;
   ScalarArrayType drhodt_in;

   VectorArrayType r_in;
   VectorArrayType v_in;
   ScalarArrayType rho_in;


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

   template< typename FluidPointer >
   void
   midpointResiduals( FluidPointer& fluid )
   {

   }

   template< typename FluidPointer >
   void
   coreector( RealType dt, FluidPointer& fluid )
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

