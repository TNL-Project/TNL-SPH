#include "OpenBoundaryConditions.h"
#include "SPH/Models/WCSPH_BI/BoundaryConditionsTypes.h"

namespace TNL {
namespace SPH {

template< typename SPHConfig, typename ModelConfig >
template< typename OpenBoundaryPointer >
typename OpenBoundaryConditionsBuffers< SPHConfig, ModelConfig >::GlobalIndexType
OpenBoundaryConditionsBuffers< SPHConfig, ModelConfig >::moveInletBufferParticles( RealType dt,
                                                                                   OpenBoundaryPointer& openBoundary )
{
   const GlobalIndexType numberOfBufferParticles = openBoundary->getParticles()->getNumberOfParticles();

   auto view_r_buffer = openBoundary->getParticles()->getPoints().getView();
   auto view_v_buffer = openBoundary->getVariables()->v.getView();
   auto view_inletMark = openBoundary->getVariables()->particleMark.getView();

   const VectorType inletOrientation = openBoundary->parameters.orientation;
   const VectorType bufferPosition = openBoundary->parameters.position;

   view_inletMark = 1;

   auto moveBufferParticles = [=] __cuda_callable__ ( int i ) mutable
   {
      view_r_buffer[ i ] += view_v_buffer[ i ] * dt;
      const VectorType r = view_r_buffer[ i ];
      const VectorType r_relative = bufferPosition - r;

      if( ( r_relative, inletOrientation ) <= 0.f )
         view_inletMark[ i ] = 0;

      return view_inletMark[ i ];
   };
   const GlobalIndexType numberOfNotRetyped = Algorithms::reduce< DeviceType >(
         0, numberOfBufferParticles, moveBufferParticles, TNL::Plus() );
   const GlobalIndexType numberOfRetyped = numberOfBufferParticles - numberOfNotRetyped;

   return numberOfRetyped;
}

template< typename SPHConfig, typename ModelConfig >
template< typename OpenBoundaryPointer >
typename OpenBoundaryConditionsBuffers< SPHConfig, ModelConfig >::GlobalIndexType
OpenBoundaryConditionsBuffers< SPHConfig, ModelConfig >::moveOutletBufferParticles( RealType dt,
                                                                                    OpenBoundaryPointer& openBoundary )
{
   GlobalIndexType numberOfBufferParticles = openBoundary->getParticles()->getNumberOfParticles();

   auto view_r_buffer = openBoundary->getParticles()->getPoints().getView();
   auto view_v_buffer = openBoundary->getVariables()->v.getView();
   auto view_inletMark = openBoundary->getVariables()->particleMark.getView();

   const VectorType inletOrientation = openBoundary->parameters.orientation;
   const VectorType bufferWidth = openBoundary->parameters.bufferWidth;
   const VectorType bufferPosition = openBoundary->parameters.position;

   view_inletMark = 0;

   auto moveBufferParticles = [=] __cuda_callable__ ( int i ) mutable
   {
      view_r_buffer[ i ] += view_v_buffer[ i ] * dt; //view_r_buffer[ i ] += ( view_v_buffer[ i ], orientation ) * orientation * dt;
      const VectorType r = view_r_buffer[ i ];
      const VectorType r_relative = bufferPosition - r;

      if( ( r_relative, inletOrientation ) > bufferWidth[ 0 ] )
         view_inletMark[ i ] = 1;

      return view_inletMark[ i ];
   };
   const GlobalIndexType removeFromBufferCount = Algorithms::reduce< DeviceType >(
         0, numberOfBufferParticles, moveBufferParticles, TNL::Plus() );

   return removeFromBufferCount;
}

template< typename SPHConfig, typename ModelConfig >
template< typename OpenBoundaryPointer >
void
OpenBoundaryConditionsBuffers< SPHConfig, ModelConfig >::sortBufferParticlesByMark( OpenBoundaryPointer& openBoundary )
{
   const GlobalIndexType numberOfBufferParticles = openBoundary->getParticles()->getNumberOfParticles();

   auto view_r_buffer = openBoundary->getParticles()->getPoints().getView();
   auto view_v_buffer = openBoundary->getVariables()->v.getView();
   auto view_rho_buffer = openBoundary->getVariables()->rho.getView();
   auto view_inletMark = openBoundary->getVariables()->particleMark.getView();

   //sort particles by mark
   using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< typename SPHConfig::DeviceType >;
   ThrustDeviceType thrustDevice;
   thrust::sort_by_key( thrustDevice,
                        view_inletMark.getArrayData(),
                        view_inletMark.getArrayData() + numberOfBufferParticles,
                        thrust::make_zip_iterator( thrust::make_tuple( view_r_buffer.getArrayData(),
                                                                       view_v_buffer.getArrayData(),
                                                                       view_rho_buffer.getArrayData() ) ) );
}

template< typename SPHConfig, typename ModelConfig >
template< typename FluidPointer, typename OpenBoundaryPointer >
void
OpenBoundaryConditionsBuffers< SPHConfig, ModelConfig >::convertBufferToFluid( FluidPointer& fluid,
                                                                               OpenBoundaryPointer& openBoundary,
                                                                               OpenBoundaryConfig& openBoundaryParams,
                                                                               const GlobalIndexType numberOfRetyped )
{
   const GlobalIndexType numberOfParticle = fluid->getParticles()->getNumberOfParticles();

   auto view_r_buffer = openBoundary->getParticles()->getPoints().getView();
   auto view_v_buffer = openBoundary->getVariables()->v.getView();
   auto view_rho_buffer = openBoundary->getVariables()->rho.getView();

   const VectorType inletOrientation = openBoundary->parameters.orientation;
   const VectorType bufferWidth = openBoundary->parameters.bufferWidth;
   // get last referential index of existing particles
   const GlobalIndexType highestReferentialIdx = fluid->getVariables()->highestReferentialIdx;

   auto view_r_fluid = fluid->getParticles()->getPoints().getView();
   auto view_v_fluid = fluid->getVariables()->v.getView();
   auto view_rho_fluid = fluid->getVariables()->rho.getView();
   auto view_rho_old = fluid->getIntegratorVariables()->rho_old.getView();
   auto view_v_old = fluid->getIntegratorVariables()->v_old.getView();
   auto view_referentialIdx = fluid->getVariables()->referentialIdx.getView();

   auto createNewFluidParticles = [=] __cuda_callable__ ( int i ) mutable
   {
      view_r_fluid[ numberOfParticle + i ] = view_r_buffer[ i ];
      view_rho_fluid[ numberOfParticle + i ] = view_rho_buffer[ i ];
      view_v_fluid[ numberOfParticle + i ] = view_v_buffer[ i ];

      view_rho_old[ numberOfParticle + i ] = view_rho_buffer[ i ];
      view_v_old[ numberOfParticle + i ] = view_v_buffer[ i ];

      view_referentialIdx[ numberOfParticle + i ] = highestReferentialIdx + i;

      //const VectorType r_relative = bufferPosition - view_r_buffer[ i ];
      //const VectorType newBufferParticle = view_r_buffer[ i ] - ( r_relative, inletOrientation ) * inletOrientation - bufferWidth[ 0 ] * inletOrientation;
      const VectorType newBufferParticle = view_r_buffer[ i ] - bufferWidth[ 0 ] * inletOrientation;

      view_r_buffer[ i ] = newBufferParticle;
   };
   Algorithms::parallelFor< DeviceType >( 0, numberOfRetyped, createNewFluidParticles );

   //Update number of particles
   fluid->getParticles()->setNumberOfParticles( numberOfParticle + numberOfRetyped );
   fluid->getVariables()->highestReferentialIdx += numberOfRetyped;
}

template< typename SPHConfig, typename ModelConfig >
template< typename FluidPointer, typename OpenBoundaryPointer >
typename OpenBoundaryConditionsBuffers< SPHConfig, ModelConfig >::GlobalIndexType
OpenBoundaryConditionsBuffers< SPHConfig, ModelConfig >::getFluidParticlesEnteringOutlet( FluidPointer& fluid,
                                                                                          OpenBoundaryPointer& openBoundary )
{
   const GlobalIndexType numberOfBufferParticles = openBoundary->getParticles()->getNumberOfParticles();

   const VectorType inletOrientation = openBoundary->parameters.orientation;
   const VectorType bufferPosition = openBoundary->parameters.position;

   auto view_r_fluid = fluid->getParticles()->getPoints().getView();

   const auto zoneParticleIndices_view = openBoundary->zone.getParticlesInZone().getConstView();
   const GlobalIndexType numberOfZoneParticles = openBoundary->zone.getNumberOfParticles();

   auto receivingParticleMark_view = openBoundary->getVariables()->receivingParticleMark.getView();
   receivingParticleMark_view = INT_MAX;

   auto checkFluidParticles = [=] __cuda_callable__ ( int i ) mutable
   {
      const GlobalIndexType p = zoneParticleIndices_view[ i ];
      const VectorType r = view_r_fluid[ p ];
      const VectorType r_relative = bufferPosition - r;

      if( ( r_relative, inletOrientation ) > 0 ){
         receivingParticleMark_view[ i ] = p;
         return 1;
      }
      return 0;
   };
   const GlobalIndexType fluidToBufferCount = Algorithms::reduce< DeviceType >(
         0, numberOfZoneParticles, checkFluidParticles, TNL::Plus() );

   const GlobalIndexType rangeToSort = ( numberOfZoneParticles > numberOfBufferParticles ) ? numberOfZoneParticles : numberOfBufferParticles;
   using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< typename SPHConfig::DeviceType >;
   ThrustDeviceType thrustDevice;
   thrust::sort( thrustDevice,
                 receivingParticleMark_view.getArrayData(),
                 receivingParticleMark_view.getArrayData() + rangeToSort );

   return fluidToBufferCount;
}

template< typename SPHConfig, typename ModelConfig >
template< typename FluidPointer, typename OpenBoundaryPointer >
void
OpenBoundaryConditionsBuffers< SPHConfig, ModelConfig >::convertFluidToBuffer( FluidPointer& fluid,
                                                                               OpenBoundaryPointer& openBoundary,
                                                                               const GlobalIndexType fluidToBufferCount )
{
   const GlobalIndexType numberOfBufferParticles = openBoundary->getParticles()->getNumberOfParticles();
   auto receivingParticleMark_view = openBoundary->getVariables()->receivingParticleMark.getView();

   auto view_r_buffer = openBoundary->getParticles()->getPoints().getView();
   auto view_v_buffer = openBoundary->getVariables()->v.getView();
   auto view_rho_buffer = openBoundary->getVariables()->rho.getView();

   auto view_r_fluid = fluid->getParticles()->getPoints().getView();
   auto view_v_fluid = fluid->getVariables()->v.getView();
   auto view_rho_fluid = fluid->getVariables()->rho.getView();

   auto retypeFluidToOutlet = [=] __cuda_callable__ ( int i ) mutable
   {
      const GlobalIndexType p = receivingParticleMark_view[ i ];

      view_r_buffer[ numberOfBufferParticles + i ] = view_r_fluid[ p ];
      view_rho_buffer[ numberOfBufferParticles + i ] = view_rho_fluid[ p ];
      view_v_buffer[ numberOfBufferParticles + i ] = view_v_fluid[ p ];

      view_r_fluid[ p ] = FLT_MAX;
   };
   Algorithms::parallelFor< DeviceType >( 0, fluidToBufferCount, retypeFluidToOutlet );

   openBoundary->getParticles()->setNumberOfParticles( numberOfBufferParticles + fluidToBufferCount );
   fluid->getParticles()->setNumberOfParticlesToRemove( fluid->getParticles()->getNumberOfParticlesToRemove() + fluidToBufferCount );

}

template< typename SPHConfig, typename ModelConfig >
template< typename FluidPointer,
          typename OpenBoundaryPointer >
void
OpenBoundaryConditionsBuffers< SPHConfig, ModelConfig >::applyOpenBoundary( RealType dt,
                                                                            FluidPointer& fluid,
                                                                            OpenBoundaryPointer& openBoundary,
                                                                            OpenBoundaryConfig& openBoundaryParams )
{
   if( openBoundaryParams.type == WCSPH_BCTypes::OpenBoundaryConditionsType::Inlet )
      applyInletBoundaryCondition( dt, fluid, openBoundary, openBoundaryParams );
   else if( openBoundaryParams.type == WCSPH_BCTypes::OpenBoundaryConditionsType::Outlet )
      applyOuletBoundaryCondition( dt, fluid, openBoundary, openBoundaryParams );
   else
   {
      std::cerr << "Invalid open boundary type." << std::endl;
   }
}

template< typename SPHConfig, typename ModelConfig >
template< typename FluidPointer,
          typename OpenBoundaryPointer >
void
OpenBoundaryConditionsBuffers< SPHConfig, ModelConfig >::applyInletBoundaryCondition( RealType dt,
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

template< typename SPHConfig, typename ModelConfig >
template< typename FluidPointer,
          typename OpenBoundaryPointer >
void
OpenBoundaryConditionsBuffers< SPHConfig, ModelConfig >::applyOuletBoundaryCondition( RealType dt,
                                                                                      FluidPointer& fluid,
                                                                                      OpenBoundaryPointer& openBoundary,
                                                                                      OpenBoundaryConfig& openBoundaryParams )
{
   //Remove leaving buffer particles:
   const GlobalIndexType bufferToVoidCount = moveOutletBufferParticles( dt, openBoundary );

   sortBufferParticlesByMark( openBoundary );
   openBoundary->getParticles()->setNumberOfParticles( openBoundary->getParticles()->getNumberOfParticles() - bufferToVoidCount );

   //Convert fluid to buffer;
   const GlobalIndexType fluidToBufferCount = getFluidParticlesEnteringOutlet( fluid, openBoundary );
   if( fluidToBufferCount == 0 )
      return;
   convertFluidToBuffer( fluid, openBoundary, fluidToBufferCount );
}

template< typename SPHConfig, typename ModelConfig >
template< typename FluidPointer, typename PeriodicBoundaryPatch >
void
OpenBoundaryConditionsBuffers< SPHConfig, ModelConfig >::periodicityParticleTransfer( FluidPointer& fluid,
                                                                                      PeriodicBoundaryPatch& periodicPatch )
{
   const auto zoneParticleIndices_view = periodicPatch->particleZone.getParticlesInZone().getConstView();
   const GlobalIndexType numberOfZoneParticles = periodicPatch->particleZone.getNumberOfParticles();

   auto view_r_fluid = fluid->getParticles()->getPoints().getView();

   const VectorType posShift = periodicPatch->config.shift;
   const VectorType bufferPosition = periodicPatch->config.position;
   const VectorType bufferOrientation = periodicPatch->config.orientation;

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
} // TNL
