#include "OpenBoundaryConditions.h"
#include "SPH/Models/WCSPH_BI/BoundaryConditionsTypes.h"

namespace TNL {
namespace SPH {

template< typename SPHConfig >
template< typename OpenBoundaryPointer >
typename OpenBoundaryConditionsBuffers< SPHConfig >::GlobalIndexType
OpenBoundaryConditionsBuffers< SPHConfig >::moveInletBufferParticles( RealType dt, OpenBoundaryPointer& openBoundary )
{
   const GlobalIndexType numberOfBufferParticles = openBoundary->particles->getNumberOfParticles();

   auto view_r_buffer = openBoundary->particles->getPoints().getView();
   auto view_v_buffer = openBoundary->variables->v.getView();
   auto view_inletMark = openBoundary->variables->particleMark.getView();

   const VectorType inletOrientation = openBoundary->parameters.orientation;
   const VectorType bufferPosition = openBoundary->parameters.position;

   view_inletMark = 1;

   auto moveBufferParticles = [ = ] __cuda_callable__( int i ) mutable
   {
      view_r_buffer[ i ] += view_v_buffer[ i ] * dt;
      const VectorType r = view_r_buffer[ i ];
      const VectorType r_relative = bufferPosition - r;

      if( ( r_relative, inletOrientation ) <= 0.f )
         view_inletMark[ i ] = 0;

      return view_inletMark[ i ];
   };
   const GlobalIndexType numberOfNotRetyped =
      Algorithms::reduce< DeviceType >( 0, numberOfBufferParticles, moveBufferParticles, TNL::Plus() );
   const GlobalIndexType numberOfRetyped = numberOfBufferParticles - numberOfNotRetyped;

   return numberOfRetyped;
}

template< typename SPHConfig >
template< typename OpenBoundaryPointer >
typename OpenBoundaryConditionsBuffers< SPHConfig >::GlobalIndexType
OpenBoundaryConditionsBuffers< SPHConfig >::moveOutletBufferParticles( RealType dt, OpenBoundaryPointer& openBoundary )
{
   GlobalIndexType numberOfBufferParticles = openBoundary->particles->getNumberOfParticles();

   auto view_r_buffer = openBoundary->particles->getPoints().getView();
   auto view_v_buffer = openBoundary->variables->v.getView();
   auto view_inletMark = openBoundary->variables->particleMark.getView();

   const VectorType inletOrientation = openBoundary->parameters.orientation;
   const VectorType bufferWidth = openBoundary->parameters.bufferWidth;
   const VectorType bufferPosition = openBoundary->parameters.position;

   view_inletMark = 0;

   auto moveBufferParticles = [ = ] __cuda_callable__( int i ) mutable
   {
      view_r_buffer[ i ] += view_v_buffer[ i ] * dt;
      //view_r_buffer[ i ] += ( view_v_buffer[ i ], orientation ) * orientation * dt;
      const VectorType r = view_r_buffer[ i ];
      const VectorType r_relative = bufferPosition - r;

      if( ( r_relative, inletOrientation ) > bufferWidth[ 0 ] )
         view_inletMark[ i ] = 1;

      return view_inletMark[ i ];
   };
   const GlobalIndexType removeFromBufferCount =
      Algorithms::reduce< DeviceType >( 0, numberOfBufferParticles, moveBufferParticles, TNL::Plus() );

   return removeFromBufferCount;
}

template< typename SPHConfig >
template< typename OpenBoundaryPointer >
void
OpenBoundaryConditionsBuffers< SPHConfig >::sortBufferParticlesByMark( OpenBoundaryPointer& openBoundary )
{
   const GlobalIndexType numberOfBufferParticles = openBoundary->particles->getNumberOfParticles();

   auto view_r_buffer = openBoundary->particles->getPoints().getView();
   auto view_v_buffer = openBoundary->variables->v.getView();
   auto view_rho_buffer = openBoundary->variables->rho.getView();
   auto view_inletMark = openBoundary->variables->particleMark.getView();

   //sort particles by mark
   using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< typename SPHConfig::DeviceType >;
   ThrustDeviceType thrustDevice;
   thrust::sort_by_key( thrustDevice,
                        view_inletMark.getArrayData(),
                        view_inletMark.getArrayData() + numberOfBufferParticles,
                        thrust::make_zip_iterator( thrust::make_tuple(
                           view_r_buffer.getArrayData(), view_v_buffer.getArrayData(), view_rho_buffer.getArrayData() ) ) );
}

template< typename SPHConfig >
template< typename FluidPointer, typename OpenBoundaryPointer >
void
OpenBoundaryConditionsBuffers< SPHConfig >::convertBufferToFluid( FluidPointer& fluid,
                                                                  OpenBoundaryPointer& openBoundary,
                                                                  OpenBoundaryConfig& openBoundaryParams,
                                                                  const GlobalIndexType numberOfRetyped )
{
   const GlobalIndexType numberOfParticle = fluid->particles->getNumberOfParticles();

   auto view_r_buffer = openBoundary->particles->getPoints().getView();
   auto view_v_buffer = openBoundary->variables->v.getView();
   auto view_rho_buffer = openBoundary->variables->rho.getView();

   const VectorType inletOrientation = openBoundary->parameters.orientation;
   const VectorType bufferWidth = openBoundary->parameters.bufferWidth;
   const VectorType inletConstVelocity = openBoundaryParams.velocity;
   const RealType inletConstDensity = openBoundaryParams.density;

   auto view_r_fluid = fluid->particles->getPoints().getView();
   auto view_v_fluid = fluid->variables->v.getView();
   auto view_rho_fluid = fluid->variables->rho.getView();
   auto view_rho_old = fluid->integratorVariables->rho_old.getView();
   auto view_v_old = fluid->integratorVariables->v_old.getView();

   auto createNewFluidParticles = [ = ] __cuda_callable__( int i ) mutable
   {
      view_r_fluid[ numberOfParticle + i ] = view_r_buffer[ i ];
      view_rho_fluid[ numberOfParticle + i ] = view_rho_buffer[ i ];
      view_v_fluid[ numberOfParticle + i ] = view_v_buffer[ i ];

      view_rho_old[ numberOfParticle + i ] = view_rho_buffer[ i ];
      view_v_old[ numberOfParticle + i ] = view_v_buffer[ i ];

      //const VectorType r_relative = bufferPosition - view_r_buffer[ i ];
      //const VectorType newBufferParticle = view_r_buffer[ i ] - ( r_relative, inletOrientation ) * inletOrientation -
      //bufferWidth[ 0 ] * inletOrientation;
      const VectorType newBufferParticle = view_r_buffer[ i ] - bufferWidth[ 0 ] * inletOrientation;

      view_r_buffer[ i ] = newBufferParticle;
      view_v_buffer[ i ] = inletConstVelocity;
      //view_rho_buffer[ i ] = inletConstDensity; //TODO: Keep the density same.
   };
   Algorithms::parallelFor< DeviceType >( 0, numberOfRetyped, createNewFluidParticles );

   //Update number of particles
   fluid->particles->setNumberOfParticles( numberOfParticle + numberOfRetyped );
   fluid->particles->setLastActiveParticle( fluid->particles->getLastActiveParticle() + numberOfRetyped );
   fluid->setLastActiveParticle( fluid->getLastActiveParticle() + numberOfRetyped );
}

template< typename SPHConfig >
template< typename FluidPointer, typename OpenBoundaryPointer >
typename OpenBoundaryConditionsBuffers< SPHConfig >::GlobalIndexType
OpenBoundaryConditionsBuffers< SPHConfig >::getFluidParticlesEnteringOutlet( FluidPointer& fluid,
                                                                             OpenBoundaryPointer& openBoundary )
{
   const GlobalIndexType numberOfBufferParticles = openBoundary->particles->getNumberOfParticles();

   const VectorType inletOrientation = openBoundary->parameters.orientation;
   const VectorType bufferPosition = openBoundary->parameters.position;

   auto view_r_fluid = fluid->particles->getPoints().getView();

   const auto zoneParticleIndices_view = openBoundary->zone.getParticlesInZone().getConstView();
   const GlobalIndexType numberOfZoneParticles = openBoundary->zone.getNumberOfParticles();

   auto receivingParticleMark_view = openBoundary->variables->receivingParticleMark.getView();
   receivingParticleMark_view = INT_MAX;

   auto checkFluidParticles = [ = ] __cuda_callable__( int i ) mutable
   {
      const GlobalIndexType p = zoneParticleIndices_view[ i ];
      const VectorType r = view_r_fluid[ p ];
      const VectorType r_relative = bufferPosition - r;

      if( ( r_relative, inletOrientation ) > 0 ) {
         receivingParticleMark_view[ i ] = p;
         return 1;
      }
      return 0;
   };
   const GlobalIndexType fluidToBufferCount =
      Algorithms::reduce< DeviceType >( 0, numberOfZoneParticles, checkFluidParticles, TNL::Plus() );

   const GlobalIndexType rangeToSort =
      ( numberOfZoneParticles > numberOfBufferParticles ) ? numberOfZoneParticles : numberOfBufferParticles;
   using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< typename SPHConfig::DeviceType >;
   ThrustDeviceType thrustDevice;
   thrust::sort(
      thrustDevice, receivingParticleMark_view.getArrayData(), receivingParticleMark_view.getArrayData() + rangeToSort );

   return fluidToBufferCount;
}

template< typename SPHConfig >
template< typename FluidPointer, typename OpenBoundaryPointer >
void
OpenBoundaryConditionsBuffers< SPHConfig >::convertFluidToBuffer( FluidPointer& fluid,
                                                                  OpenBoundaryPointer& openBoundary,
                                                                  const GlobalIndexType fluidToBufferCount )
{
   const GlobalIndexType numberOfBufferParticles = openBoundary->particles->getNumberOfParticles();
   auto receivingParticleMark_view = openBoundary->variables->receivingParticleMark.getView();

   auto view_r_buffer = openBoundary->particles->getPoints().getView();
   auto view_v_buffer = openBoundary->variables->v.getView();
   auto view_rho_buffer = openBoundary->variables->rho.getView();

   auto view_r_fluid = fluid->particles->getPoints().getView();
   auto view_v_fluid = fluid->variables->v.getView();
   auto view_rho_fluid = fluid->variables->rho.getView();

   auto retypeFluidToOutlet = [ = ] __cuda_callable__( int i ) mutable
   {
      const GlobalIndexType p = receivingParticleMark_view[ i ];

      view_r_buffer[ numberOfBufferParticles + i ] = view_r_fluid[ p ];
      view_rho_buffer[ numberOfBufferParticles + i ] = view_rho_fluid[ p ];
      view_v_buffer[ numberOfBufferParticles + i ] = view_v_fluid[ p ];

      view_r_fluid[ p ] = FLT_MAX;
   };
   Algorithms::parallelFor< DeviceType >( 0, fluidToBufferCount, retypeFluidToOutlet );

   openBoundary->particles->setNumberOfParticles( numberOfBufferParticles + fluidToBufferCount );
   openBoundary->particles->setLastActiveParticle( openBoundary->particles->getLastActiveParticle() + fluidToBufferCount );
   openBoundary->setLastActiveParticle( openBoundary->getLastActiveParticle() + fluidToBufferCount );

   openBoundary->numberOfFluidParticlesToRemove = fluidToBufferCount;
}

template< typename SPHConfig >
template< typename FluidPointer, typename OpenBoundaryPointer >
void
OpenBoundaryConditionsBuffers< SPHConfig >::applyOpenBoundary( RealType dt,
                                                               FluidPointer& fluid,
                                                               OpenBoundaryPointer& openBoundary,
                                                               OpenBoundaryConfig& openBoundaryParams )
{
   if( openBoundaryParams.type == WCSPH_BCTypes::OpenBoundaryConditionsType::Inlet )
      applyInletBoundaryCondition( dt, fluid, openBoundary, openBoundaryParams );
   else if( openBoundaryParams.type == WCSPH_BCTypes::OpenBoundaryConditionsType::Outlet )
      applyOuletBoundaryCondition( dt, fluid, openBoundary, openBoundaryParams );
   else {
      std::cerr << "Invalid open boundary type." << std::endl;
   }
}

template< typename SPHConfig >
template< typename FluidPointer, typename OpenBoundaryPointer >
void
OpenBoundaryConditionsBuffers< SPHConfig >::applyInletBoundaryCondition( RealType dt,
                                                                         FluidPointer& fluid,
                                                                         OpenBoundaryPointer& openBoundary,
                                                                         OpenBoundaryConfig& openBoundaryParams )
{
   const GlobalIndexType bufferToFluidCount = moveInletBufferParticles( dt, openBoundary );

   if( bufferToFluidCount == 0 )
      return;

   sortBufferParticlesByMark( openBoundary );
   convertBufferToFluid( fluid, openBoundary, openBoundaryParams, bufferToFluidCount );
}

template< typename SPHConfig >
template< typename FluidPointer, typename OpenBoundaryPointer >
void
OpenBoundaryConditionsBuffers< SPHConfig >::applyOuletBoundaryCondition( RealType dt,
                                                                         FluidPointer& fluid,
                                                                         OpenBoundaryPointer& openBoundary,
                                                                         OpenBoundaryConfig& openBoundaryParams )
{
   //Remove leaving buffer particles:
   const GlobalIndexType bufferToVoidCount = moveOutletBufferParticles( dt, openBoundary );

   sortBufferParticlesByMark( openBoundary );
   openBoundary->particles->setNumberOfParticles( openBoundary->particles->getNumberOfParticles() - bufferToVoidCount );
   openBoundary->particles->setLastActiveParticle( openBoundary->getLastActiveParticle() - bufferToVoidCount );
   openBoundary->setLastActiveParticle( openBoundary->getLastActiveParticle() - bufferToVoidCount );

   //Convert fluid to buffer;
   const GlobalIndexType fluidToBufferCount = getFluidParticlesEnteringOutlet( fluid, openBoundary );
   if( fluidToBufferCount == 0 )
      return;
   convertFluidToBuffer( fluid, openBoundary, fluidToBufferCount );
}

template< typename SPHConfig >
template< typename FluidPointer, typename PeriodicBoundaryPatch >
void
OpenBoundaryConditionsBuffers< SPHConfig >::periodicityParticleTransfer( FluidPointer& fluid,
                                                                         PeriodicBoundaryPatch& periodicPatch )
{
   const auto zoneParticleIndices_view = periodicPatch->particleZone.getParticlesInZone().getConstView();
   const GlobalIndexType numberOfZoneParticles = periodicPatch->particleZone.getNumberOfParticles();

   auto view_r_fluid = fluid->particles->getPoints().getView();

   const VectorType posShift = periodicPatch->config.shift;
   const VectorType bufferPosition = periodicPatch->config.position;
   const VectorType bufferOrientation = periodicPatch->config.orientation;

   auto moveParticles = [ = ] __cuda_callable__( int i ) mutable
   {
      const GlobalIndexType p = zoneParticleIndices_view[ i ];
      const VectorType r = view_r_fluid[ p ];
      const VectorType r_relative = bufferPosition - r;

      if( ( r_relative, bufferOrientation ) > 0.f )
         view_r_fluid[ p ] += posShift;
   };
   Algorithms::parallelFor< DeviceType >( 0, numberOfZoneParticles, moveParticles );
}

}  //namespace SPH
}  //namespace TNL
