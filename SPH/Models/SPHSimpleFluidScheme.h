#pragma once

#include "../SPHTraits.h"

/**
 * Modules used as default.
 **/
#include "EquationOfState.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHConfig >
class SPHSimpleFluidScheme
{
public:

   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using DeviceType = typename SPHConfig::DeviceType;

   /* Default modules. */
   using EOS = TaitWeaklyCompressibleEOS< SPHConfig >;

   /**
    * Constructor.
    */
   SPHSimpleFluidScheme() = default;

   /**
    * General interaction function.
    */
   template< typename FluidPointer, typename BoudaryPointer, typename NeighborSearchPointer, typename SPHKernelFunction >
   virtual void
   Interaction( FluidPointer& fluid, BoudaryPointer& boundary );

   /**
    * Compute pressure from density.
    * This holds for all the schemes for SimpleFluid model.
    */
   template< typename VariablesPointer& variablesPointer, typename EquationOfState = TaitWeaklyCompressibleEOS< SPHConfig > >
   virtual void
   ComputePressureFromDensity( VariablesPointer& variables, GlobalIndexType numberOfParticles )
   {
      auto view_rho = variables->rho.getView();
      auto view_p = variables->p.getView();

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         view_p[ i ] = EquationOfState::DensityToPressure( view_rho[ i ] );
      };
      Algorithms::parallelFor< DeviceType >( 0, numberOfParticles, init );
   }

};

} // SPH
} // ParticleSystem
} // TNL

