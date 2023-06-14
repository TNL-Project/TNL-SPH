#include <TNL/Meshes/Writers/VTKWriter.h>
#include <TNL/Containers/NDArray.h>

#include "../SPHTraits.h"
#include "../../Particles/neighborSearchLoop.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHConfig, typename SPHSimulation >
class InterpolateToGrid
{
public:

   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using DeviceType = typename SPHConfig::DeviceType;

   using LocalIndexType = typename SPHTraitsType::LocalIndexType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using IndexVectorType = typename SPHTraitsType::IndexVectorType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;


   using FluidPointer = typename SPHSimulation::FluidPointer;
   using BoundaryPointer = typename SPHSimulation::BoundaryPointer;
   using Variables = typename SPHSimulation::FluidVariables;
   using VariablesPointer = typename Pointers::SharedPointer< Variables, DeviceType >;

   using GridType = Meshes::Grid< 2, RealType, DeviceType, GlobalIndexType >;
   using CoordinatesType = typename GridType::CoordinatesType;

   //temp
   using ParticleSystem = typename SPHSimulation::ParticleSystemType;

   template< typename InterpolationConfigType >
   InterpolateToGrid( InterpolationConfigType configInit )
   : variables( configInit.gridSize[ 0 ] * configInit.gridSize[ 1 ] ),
     gridDimension( configInit.gridSize )
   {
      interpolationGrid.setOrigin( configInit.gridOrigin );
      interpolationGrid.setDimensions( configInit.gridSize );
      interpolationGrid.setSpaceSteps( configInit.gridStep );
   }

   template< typename SPHKernelFunction, typename SPHState >
   void
   interpolate( FluidPointer& fluid, BoundaryPointer& boundary, SPHState& sphState );

   void
   save( const std::string outputFileName );

protected:

   GridType interpolationGrid;
   CoordinatesType gridDimension;
   VariablesPointer variables;


};

template< typename SPHConfig, typename SPHSimulation >
class SensorInterpolation
{
   public:
   using DeviceType = typename SPHConfig::DeviceType;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;

   using LocalIndexType = typename SPHTraitsType::LocalIndexType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using IndexVectorType = typename SPHTraitsType::IndexVectorType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;

   using FluidPointer = typename SPHSimulation::FluidPointer;
   using BoundaryPointer = typename SPHSimulation::BoundaryPointer;

   //temp
   using ParticleSystem = typename SPHSimulation::ParticleSystemType;

   using SensorsDataArray = Containers::NDArray< RealType,  // Value
                                                 Containers::SizesHolder< int, 0, 0 >,     // SizesHolder
                                                 std::index_sequence< 0, 1 >,  // Permutation
                                                 DeviceType >;         // Device

   SensorInterpolation( GlobalIndexType numberOfSavedSteps, std::vector< VectorType >& sensorsPoints )
   : sensorPositions( sensorsPoints ),
     numberOfSensors( sensorsPoints.size() ),
     numberOfSavedSteps( numberOfSavedSteps )
   {
      sensors.setSizes( numberOfSavedSteps, numberOfSensors );
   }

   template<typename SPHKernelFunction, typename EOS, typename SPHState >
   void
   interpolate( FluidPointer& fluid, BoundaryPointer& boundary, SPHState& sphState );

   void
   save( const std::string outputFileName );

protected:

   GlobalIndexType numberOfSensors;
   SensorsDataArray sensors;
   VectorArrayType sensorPositions;

   GlobalIndexType numberOfSavedSteps;
   GlobalIndexType sensorIndexer = 0;
};

template< typename SPHConfig, typename SPHSimulation >
class SensorWaterLevel
{
   public:
   using DeviceType = typename SPHConfig::DeviceType;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;

   using LocalIndexType = typename SPHTraitsType::LocalIndexType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using IndexVectorType = typename SPHTraitsType::IndexVectorType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;
   using IndexArrayType = typename SPHTraitsType::IndexArrayType;

   using FluidPointer = typename SPHSimulation::FluidPointer;
   using BoundaryPointer = typename SPHSimulation::BoundaryPointer;

   //temp
   using ParticleSystem = typename SPHSimulation::ParticleSystemType;

   using SensorsDataArray = Containers::NDArray< RealType,  // Value
                                                 Containers::SizesHolder< int, 0, 0 >,     // SizesHolder
                                                 std::index_sequence< 0, 1 >,  // Permutation
                                                 Devices::Host >;         // Device - store data on HOST

   SensorWaterLevel( GlobalIndexType numberOfSavedSteps,
                     std::vector< VectorType >& sensorsPoints,
                     RealType levelIncrement,
                     VectorType direction,
                     RealType startLevel,
                     RealType endLevel )
   : sensorPositions( sensorsPoints ), numberOfSensors( sensorsPoints.size() ), numberOfSavedSteps( numberOfSavedSteps ),
     levelIncrement( levelIncrement ), direction( direction )
   {
      sensors.setSizes( numberOfSavedSteps, numberOfSensors );

      //New
      numberOfLevels = TNL::ceil( ( endLevel - startLevel ) / levelIncrement );
      levels.setSize( numberOfLevels );
   }

   template< typename SPHKernelFunction, typename EOS, typename SPHState >
   void
   interpolate( FluidPointer& fluid, BoundaryPointer& boundary, SPHState& sphState );

   void
   save( const std::string outputFileName );

protected:

   GlobalIndexType numberOfSensors;
   SensorsDataArray sensors;
   VectorArrayType sensorPositions;

   GlobalIndexType numberOfSavedSteps;
   GlobalIndexType sensorIndexer = 0;

   //New
   IndexArrayType levels;
   GlobalIndexType numberOfLevels; //numberOfLevels = TNL::ceil( TNL::norm( endPoints - startPoint ) / levelIncrement );
   RealType levelIncrement;
   VectorType direction;

   RealType startLevel;
   RealType endLevel;
};

} // SPH
} // ParticleSystem
} // TNL

#include "Measuretool.hpp"

