#include <SPH/shared/Measuretool.h>
#include <cmath>

namespace TNL {
namespace SPH {
namespace features {

template< typename SPHSimulation, typename FluidPointer, typename BoundaryPointer, typename OpenBoundaryPointer, typename ModelParams, typename RealType,typename VectorType >
void updateOpenBCHydrostaticProfile( FluidPointer& fluid,
                                     BoundaryPointer& boundary,
                                     OpenBoundaryPointer& openBoundary,
                                     ModelParams& modelParams,
                                     const VectorType& referentialPoint,
                                     const VectorType& direction,
                                     const RealType& waterLevelTrashold,
                                     const RealType& kineticPressure = 0.f )
{
   using DeviceType = typename ModelParams::SPHConfig::DeviceType;
   using EOS = typename ModelParams::EOS;
   using GlobalIndexType = typename ModelParams::SPHTraitsType::GlobalIndexType;

   if( fluid->getNumberOfParticles() == 0 )
      return;

   ////obtain point water level - using max z particle component
   //const auto view_points_fluid = fluid->getParticles()->getPoints().getConstView();
   //auto fetch = [=] __cuda_callable__ ( int i )
   //{
   //   return ( view_points_fluid[ i ], direction );
   //};
   //const RealType waterLevel = Algorithms::reduce< DeviceType >( 0, fluid->getNumberOfParticles(), fetch, TNL::Max() );

   ////obtain point water level - using measuretool
   ////const RealType waterLevel = ...;
   using SensorWaterLevel = SensorWaterLevel< typename ModelParams::SPHConfig, SPHSimulation >;
   SensorWaterLevel waterLevelSensor;
   std::vector< VectorType > pointsToMeasureWaterLevel = { referentialPoint };
   waterLevelSensor.init( pointsToMeasureWaterLevel,
                          1,
                          modelParams.dp,
                          direction,
                          ( referentialPoint, direction ),
                          0.5 );
   waterLevelSensor.template interpolate< typename ModelParams::KernelFunction, typename ModelParams::EOS >(
         fluid, boundary, modelParams );
   const RealType waterLevel = waterLevelSensor.getSensorData().getElement( 0, 0 );

   //update inlet profile
   const RealType speedOfSound = modelParams.speedOfSound;
   const RealType rho0 = modelParams.rho0;
   const VectorType gravity = modelParams.gravity;
   const RealType gravityMagnitude = TNL::l2Norm( gravity );
   const typename EOS::ParamsType eosParams( modelParams );
   const RealType kineticDensity = EOS::pressureToDensity( kineticPressure, eosParams ) - rho0;
   std::cout << "direction/nptcs: " << direction << "/" << fluid->getNumberOfParticles() << " waterLevel: " << waterLevel << " kineticDensity: " << kineticDensity << std::endl;

   const auto view_points_openBound = openBoundary->getParticles()->getPoints().getConstView();
   auto view_rho_openBound = openBoundary->getVariables()->rho.getView();

   // update only if water level is above defined trashold
   if( waterLevel > waterLevelTrashold ){
      auto updateParticleDensity = [=] __cuda_callable__ ( GlobalIndexType i ) mutable
      {
         const VectorType r_i = view_points_openBound[ i ];
         const RealType h_depth = waterLevel - ( r_i, direction );
         const RealType hydrostaticPressure = rho0 * gravityMagnitude * h_depth;
         const RealType hydrostaticDensity = EOS::pressureToDensity( hydrostaticPressure, eosParams );

         view_rho_openBound[ i ] = hydrostaticDensity + kineticDensity;
      };
      openBoundary->getParticles()->forAll( updateParticleDensity );
   }
   //else{
   //   view_rho_openBound = rho0 + kineticDensity;
   //}
}

} // features
} // SPH
} // TNL

