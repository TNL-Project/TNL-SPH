#include <TNL/Meshes/Writers/VTKWriter.h>
#include <TNL/Containers/NDArray.h>
#include <TNL/Pointers/SharedPointer.h>

#include "../SPHTraits.h"
#include "../../Particles/neighborSearchLoop.h"
#include "TNL/Config/ParameterContainer.h"

namespace TNL {
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

   //There should be dimension only 2
   using GridType = Meshes::Grid< SPHConfig::spaceDimension, RealType, DeviceType, GlobalIndexType >;
   using CoordinatesType = typename GridType::CoordinatesType;

   //temp
   using ParticleSystem = typename SPHSimulation::ParticlesType;

   InterpolateToGrid() : variables() {}

   InterpolateToGrid( const VectorType& gridSize, const VectorType& gridOrigin, const IndexVectorType& stepSize )
   : variables( gridSize[ 0 ] * gridSize[ 1 ] ), gridDimension( gridSize ), gridStep( gridStep )
   {
      interpolationGrid.setOrigin( gridOrigin );
      interpolationGrid.setDimensions( gridSize );
      interpolationGrid.setSpaceSteps( gridStep );
   }

   void
   init( TNL::Config::ParameterContainer& parameters, const std::string& prefix )
   {
      const VectorType gridSize = parameters.getXyz< IndexVectorType >( prefix +"gridSize" );
      variables->setSize( gridSize[ 0 ] * gridSize[ 1 ] );
      interpolationGrid.setOrigin( parameters.getXyz< VectorType >( prefix + "gridOrigin" ) );
      interpolationGrid.setDimensions( gridSize );
      interpolationGrid.setSpaceSteps( parameters.getXyz< VectorType >( prefix + "gridStep" ) );

      gridDimension = gridSize;
      gridStep = parameters.getXyz< VectorType >( prefix + "gridStep" );
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

   VectorType gridStep;


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
   using ParticleSystem = typename SPHSimulation::ParticlesType;

   using SensorsDataArray = Containers::NDArray< RealType,  // Value
                                                 Containers::SizesHolder< int, 0, 0 >,     // SizesHolder
                                                 std::index_sequence< 0, 1 >,  // Permutation
                                                 DeviceType >;         // Device

   SensorInterpolation() : sensorPositions( 0 ) {}

   SensorInterpolation( GlobalIndexType numberOfSavedSteps, std::vector< VectorType >& sensorsPoints )
   : sensorPositions( sensorsPoints ), numberOfSensors( sensorsPoints.size() ), numberOfSavedSteps( numberOfSavedSteps )
   {
      sensors.setSizes( numberOfSavedSteps, numberOfSensors );
   }

   void
   init( std::vector< VectorType >& points, const int numberOfSensors, const int numberOfSavedSteps, bool includeBoundary )
   {
      this->numberOfSensors = numberOfSensors;
      sensors.setSizes( numberOfSavedSteps, numberOfSensors );
      sensorPositions.setSize( numberOfSensors );
      sensorPositions = points;
      this->numberOfSavedSteps = numberOfSavedSteps;
   }

   template<typename SPHKernelFunction, typename EOS, typename SPHState >
   void
   interpolate( FluidPointer& fluid, BoundaryPointer& boundary, SPHState& sphState, bool includeBoundary );

   void
   save( const std::string outputFileName );

protected:

   GlobalIndexType numberOfSensors;
   SensorsDataArray sensors;
   VectorArrayType sensorPositions;
   bool includeBoundary = false;

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
   using ParticleSystem = typename SPHSimulation::ParticlesType;

   using SensorsDataArray = Containers::NDArray< RealType,  // Value
                                                 Containers::SizesHolder< int, 0, 0 >,     // SizesHolder
                                                 std::index_sequence< 0, 1 >,  // Permutation
                                                 Devices::Host >;         // Device - store data on HOST

   SensorWaterLevel() : sensorPositions( 0 ) {}

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

   void
   init( std::vector< VectorType >& points,
         const int numberOfSavedSteps,
         RealType levelIncrement,
         VectorType direction,
         RealType startLevel,
         RealType endLevel )
   {
      this->numberOfSensors = points.size();
      sensors.setSizes( numberOfSavedSteps, numberOfSensors );
      sensorPositions.setSize( numberOfSensors );
      sensorPositions = points;

      numberOfLevels = TNL::ceil( ( endLevel - startLevel ) / levelIncrement );
      levels.setSize( numberOfLevels );
      this->direction = direction;
      this->numberOfSavedSteps = numberOfSavedSteps;
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

   IndexArrayType levels;
   GlobalIndexType numberOfLevels;
   RealType levelIncrement;
   VectorType direction;
   RealType startLevel;
   RealType endLevel;
};

} // SPH
} // TNL

#include "Measuretool.hpp"

