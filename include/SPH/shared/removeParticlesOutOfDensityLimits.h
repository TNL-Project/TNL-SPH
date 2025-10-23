#pragma once
#include <TNL/Algorithms/reduce.h>

namespace TNL {
namespace SPH {
namespace customFunctions {

template< typename FluidPointer, typename ModelParams >
int
removeParticlesOutOfDensityLimits( FluidPointer& fluid, ModelParams& modelParams )
{
   using DeviceType = typename ModelParams::SPHConfig::DeviceType;

   const float rho0 = modelParams.rho0;
   const float rhoUpperTrashold = 1.2f * rho0;
   const float rhoLowerTrashold = 0.8f * rho0;

   auto r_view = fluid->getPoints().getView();
   auto rho_view = fluid->getVariables()->rho.getView();

   auto checkParticleDensity  = [ = ] __cuda_callable__( int i ) mutable
   {
      // if the particle is already removed, skip
      if( r_view[ 0 ] == FLT_MAX )
         return 0;

      if( ( rho_view[ i ] > rhoLowerTrashold ) && ( rho_view[ i ] < rhoUpperTrashold ) ) {
         return 0;
      }
      else {
         r_view[ i ] = FLT_MAX;
         return 1;
      }
   };
   const int numberOfParticlesToRemove =
      TNL::Algorithms::reduce< DeviceType >( 0, fluid->getNumberOfParticles(), checkParticleDensity , TNL::Plus() );
   fluid->getParticles()->setNumberOfParticlesToRemove(
         fluid->getParticles()->getNumberOfParticlesToRemove() + numberOfParticlesToRemove );

   return numberOfParticlesToRemove;
}

}  //namespace customFunctions
}  //namespace SPH
}  //namespace TNL

