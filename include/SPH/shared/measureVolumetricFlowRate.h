#pragma once

#include <fstream>
#include "../SPHTraits.h"
#include <TNL/Particles/GhostZone.h>

namespace TNL {
namespace SPH {

template< typename ParticlesConfig, typename SPHDefs >
class MeasureVolumetricFlowRate
{
public:
   using SPHConfig = typename SPHDefs::SPHConfig;
   using DeviceType = typename SPHConfig::DeviceType;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using IndexType = typename SPHConfig::GlobalIndexType;
   using IndexVectorType = Containers::StaticVector< SPHConfig::spaceDimension, IndexType >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   using ParticlesZone = TNL::ParticleSystem::ParticleZone< ParticlesConfig, DeviceType >;

   template< typename FluidPointer >
   void
   init( const FluidPointer& fluid,
         const VectorType pointMin,
         const VectorType pointMax,
         const VectorType orientation,
         const IndexType numberOfParticlesPerCell = 75 )
   {
      // turn the pointMin and pointMax to zone points
      this->pointMin = pointMin;
      this->pointMax = pointMax;
      this->orientation = orientation;

      const VectorType gridOrigin = fluid->getParticles()->getGridOrigin();
      const IndexVectorType gridDimensions = fluid->getParticles()->getGridDimensions();
      const RealType searchRadius = fluid->getParticles()->getSearchRadius();
      const VectorType zoneFirstPoint = pointMin - orientation * searchRadius;
      const VectorType zoneSecondPoint = pointMax + orientation * searchRadius;

      // initialize the zone
      zone.setNumberOfParticlesPerCell( numberOfParticlesPerCell );
      zone.assignCells( zoneFirstPoint, zoneSecondPoint, gridDimensions, gridOrigin, searchRadius );
   }

   const RealType
   getVolumetricFlowRate() const
   {
      return volumetricFlux;
   }

   template< typename FluidPointer, typename ModelParams >
   void
   measureVolumetricFlowRate( FluidPointer& fluid, ModelParams& modelParams, const RealType dt )
   {
      const RealType m = modelParams.mass;
      const VectorType p_min = pointMax;
      const VectorType normal = orientation;

      const auto r_view = fluid->getParticles()->getPoints().getConstView();
      const auto v_view = fluid->getVariables()->v.getConstView();
      const auto rho_view = fluid->getVariables()->rho.getConstView();

      zone.updateParticlesInZone( fluid->getParticles() );
      const auto zoneParticleIndices_view = zone.getParticlesInZone().getConstView();
      const IndexType numberOfZoneParticles = zone.getNumberOfParticles();

      auto countParticles = [=] __cuda_callable__ ( int i ) -> RealType
      {
         const IndexType p = zoneParticleIndices_view[ i ];
         const VectorType r_i = r_view[ p ];
         const VectorType v_i = v_view[ p ];
         const RealType rho_j = rho_view[ p ];
         const VectorType r_i_new = r_i + dt * v_i;

         // project the positions to the axis of the plane
         const RealType r_i_n_relative = ( r_i - p_min, normal );
         const RealType r_i_new_n_relative = ( r_i_new - p_min, normal );

         if( ( r_i_n_relative < 0.f ) && ( r_i_new_n_relative > 0.f ) )
            return m / rho_j;
         else if ( ( r_i_n_relative > 0.f ) && ( r_i_new_n_relative < 0.f ) )
            return ( -1.f ) * m / rho_j;
         else
            return 0.f;
      };
      const RealType volume = Algorithms::reduce< DeviceType >(
            0, numberOfZoneParticles, countParticles, TNL::Plus() );

      // average and evaluate
      intervalVolumeSum += volume;
      intervalTime += dt;

      if( intervalTime > averagingInterval ){
         volumetricFlux = intervalVolumeSum / intervalTime;
         intervalVolumeSum = 0.f;
         intervalTime = 0.f;
         cumulativeVolumetricFlux += volumetricFlux;
      }
   }

   void
   output( const std::string& outputPath, const IndexType step, const RealType time )
   {
      std::ofstream outfile;
      outfile.open( outputPath, std::ios_base::app );
      outfile << step << " " << time << " " << volumetricFlux << " " << cumulativeVolumetricFlux << std::endl;
   }

   VectorType pointMin;
   VectorType pointMax;
   VectorType orientation;

   RealType averagingInterval = 0.1;

   RealType intervalTime = 0.f;
   RealType intervalVolumeSum = 0.f;
   RealType volumetricFlux = 0.f;
   RealType cumulativeVolumetricFlux = 0.f;

   ParticlesZone zone;
};

} // SPH
} // TNL

