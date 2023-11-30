#include "Integrator.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {


template< typename SPHConfig >
template< typename FluidPointer, typename OpenBoundaryPointer >
void
VerletIntegrator< SPHConfig >::applyPeriodicBoundary( FluidPointer& fluid,
                                                      OpenBoundaryPointer& periodicBoundary1,
                                                      OpenBoundaryPointer& periodicBoundary2,
                                                      VectorType shift )
{
   copyGhostParticles( fluid, periodicBoundary1, periodicBoundary2, shift );
   VectorType shift2 =  -1.f * shift;
   copyGhostParticles( fluid, periodicBoundary2, periodicBoundary1, shift2 );
}

template< typename SPHConfig >
template< typename FluidPointer, typename OpenBoundaryPointer >
void
VerletIntegrator< SPHConfig >::copyGhostParticles( FluidPointer& fluid,
                                                   OpenBoundaryPointer& sendingBuffer,
                                                   OpenBoundaryPointer& receivingBuffer,
                                                   VectorType shift )
{
   sendingBuffer->zone.updateParticlesInZone( fluid->particles );

   const auto zoneParticleIndices_view = sendingBuffer->zone.getParticlesInZone().getConstView();
   const GlobalIndexType numberOfZoneParticles = sendingBuffer->zone.getNumberOfParticles();

   auto view_r = fluid->particles->getPoints().getView();
   auto view_v = fluid->variables->v.getView();
   auto view_rho = fluid->variables->rho.getView();

   auto view_r_recBuffer = receivingBuffer->particles->getPoints().getView();
   auto view_v_recBuffer = receivingBuffer->variables->v.getView();
   auto view_rho_recBuffer = receivingBuffer->variables->rho.getView();

   auto copyParticles = [=] __cuda_callable__ ( int i ) mutable
   {
      const GlobalIndexType p = zoneParticleIndices_view[ i ];

      view_r_recBuffer[ i ] = view_r[ p ] + shift;
      view_v_recBuffer[ i ] = view_v[ p ];
      view_rho_recBuffer[ i ] = view_rho[ p ];
   };
   Algorithms::parallelFor< DeviceType >( 0, numberOfZoneParticles, copyParticles );

   receivingBuffer->setLastActiveParticle( numberOfZoneParticles - 1 );
   receivingBuffer->particles->setLastActiveParticle( numberOfZoneParticles - 1 );
   receivingBuffer->particles->setNumberOfParticles( numberOfZoneParticles );
}

template< typename SPHConfig >
template< typename FluidPointer, typename OpenBoundaryPointer >
void
VerletIntegrator< SPHConfig >::periodicityParticleTransfer( FluidPointer& fluid,
                                                            OpenBoundaryPointer& periodicBuffer,
                                                            VectorType& posShift )

{
   const auto zoneParticleIndices_view = periodicBuffer->zone.getParticlesInZone().getConstView();
   const GlobalIndexType numberOfZoneParticles = periodicBuffer->zone.getNumberOfParticles();

   auto view_r_fluid = fluid->particles->getPoints().getView();

   const VectorType bufferPosition = periodicBuffer->parameters.position;
   const VectorType bufferOrientation = periodicBuffer->parameters.orientation;

   auto moveParticles = [=] __cuda_callable__ ( int i ) mutable
   {
      const GlobalIndexType p = zoneParticleIndices_view[ i ];
      const VectorType r = view_r_fluid[ p ];

      const VectorType r_relative = bufferPosition - r;

      if( ( r_relative, bufferOrientation ) > 0.f )
         view_r_fluid[ p ] += posShift;

   };
   Algorithms::parallelFor< DeviceType >( 0, numberOfZoneParticles, moveParticles );

}

} // SPH
} // ParticleSystem
} // TNL
