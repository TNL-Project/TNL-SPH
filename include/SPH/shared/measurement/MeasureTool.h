// SPDX-FileComment: This file is part of TNL-SPH (https://gitlab.com/tnl-project/tnl-sph)
// SPDX-License-Identifier: MIT

#pragma once

#include <TNL/Meshes/Writers/VTKWriter.h>
#include <TNL/Containers/NDArray.h>
#include <TNL/Pointers/SharedPointer.h>
#include <TNL/Config/ParameterContainer.h>
#include <TNL/Particles/CellIndexer.h>

#include <SPH/SPHTraits.h>
#include <SPH/shared/measurement/MeasurementLevelView.h>

namespace TNL::SPH {

template< typename SPHConfig, typename SPHSimulation >
class InterpolateToGrid
{
public:
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using DeviceType = typename SPHConfig::DeviceType;

   using LocalIndexType = typename SPHTraitsType::LocalIndexType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;
   using ScalarArrayType = typename SPHTraitsType::ScalarArrayType;

   using FluidPointer = typename SPHSimulation::FluidPointer;
   using BoundaryPointer = typename SPHSimulation::BoundaryPointer;
   using Variables = typename SPHSimulation::FluidVariables;
   using VariablesPointer = typename Pointers::SharedPointer< Variables, DeviceType >;
   using IndexArrayType = typename SPHTraitsType::IndexArrayType;

   //There should be dimension only 2
   using GridType = Meshes::Grid< SPHConfig::spaceDimension, RealType, DeviceType, GlobalIndexType >;
   using CoordinatesType = typename GridType::CoordinatesType;

   // typenames used in interpolateUsingGrid
   using GridCell = typename GridType::Cell;
   using GridVertex = typename GridType::Vertex;

   // typenames used in interpolateUsingParallelFor
   using CellIndexer = TNL::Particles::SimpleCellIndex< SPHConfig::spaceDimension >;

   /**
    * \brief Particle type is required to search through particle seats we want to
    * monitor with measuretool. Measuretool works with CLL search, so in case we use
    * CLLwL, for measuretool we need to load the base class.
    */
   //using ParticlesType = typename SPHSimulation::ParticlesType;
   using ParticlesType = std::conditional_t<
      SPHSimulation::ParticlesType::specifySearchedSetExplicitly(),
      typename SPHSimulation::ParticlesType::BaseType,
      typename SPHSimulation::ParticlesType >;

   InterpolateToGrid()
   : variables()
   {}

   void
   init(
      TNL::Config::ParameterContainer& parameters,
      const std::string& prefix,
      const int interpolatedGridEntity = SPHConfig::spaceDimension )
   {
      //TODO: test: interpolateGridEntity can be 0 i.e point or spaceDimension, i.e. cell
      this->interpolatedGridEntity = interpolatedGridEntity;

      interpolationGrid.setOrigin( getXyzCamelCase< VectorType >( parameters, prefix + "GridOrigin" ) );
      interpolationGrid.setDimensions( getXyzCamelCase< CoordinatesType >( parameters, prefix + "GridSize" ) );
      interpolationGrid.setSpaceSteps( getXyzCamelCase< VectorType >( parameters, prefix + "GridStep" ) );

      variables->setSize( interpolationGrid.getEntitiesCount( interpolatedGridEntity ) );
      gamma.setSize( interpolationGrid.getEntitiesCount( interpolatedGridEntity ) );
   }

   template< typename SPHKernelFunction, typename SPHState >
   void
   interpolate( FluidPointer& fluid, BoundaryPointer& boundary, SPHState& sphState );

   template< typename SPHKernelFunction, typename ModelParams, typename LevelView >
   void
   interpolate( const std::vector< LevelView >& levels, ModelParams& modelParams );

   template< typename SPHKernelFunction, typename SPHState >
   void
   interpolateUsingGrid( FluidPointer& fluid, BoundaryPointer& boundary, SPHState& sphState );

   template< typename SPHKernelFunction, typename SPHState >
   void
   interpolateUsingParallelFor( FluidPointer& fluid, BoundaryPointer& boundary, SPHState& sphState );

   void
   save( const std::string outputFileName );

protected:
   int interpolatedGridEntity;

   GridType interpolationGrid;
   CoordinatesType dimensions;
   VariablesPointer variables;
   ScalarArrayType gamma;

   VectorType gridStep;

public:
   template< typename SPHKernelFunction, typename SourcePointer >
   void
   accumulateGridSource(
      const SourcePointer& source,
      const IndexArrayType& owners,
      GlobalIndexType levelId,
      RealType h,
      RealType mass,
      const VectorType& gridOrigin,
      const VectorType& gridSpaceSteps );
};

template< typename SPHConfig, typename SPHSimulation >
class SensorInterpolation
{
public:
   using DeviceType = typename SPHConfig::DeviceType;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;

   using LocalIndexType = typename SPHTraitsType::LocalIndexType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using CoordinatesType = typename SPHTraitsType::CoordinatesType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;
   using ScalarArrayType = typename SPHTraitsType::ScalarArrayType;
   using IndexArrayType = typename SPHTraitsType::IndexArrayType;

   using FluidPointer = typename SPHSimulation::FluidPointer;
   using BoundaryPointer = typename SPHSimulation::BoundaryPointer;

   /**
    * \brief Particle type is required to search through particle seats we want to
    * monitor with measuretool. Measuretool works with CLL search, so in case we use
    * CLLwL, for measuretool we need to load the base class.
    */
   //using ParticlesType = typename SPHSimulation::ParticlesType;
   using ParticlesType = std::conditional_t<
      SPHSimulation::ParticlesType::specifySearchedSetExplicitly(),
      typename SPHSimulation::ParticlesType::BaseType,
      typename SPHSimulation::ParticlesType >;

   using SensorsDataArray = Containers::NDArray<
      RealType,                              // Value
      Containers::SizesHolder< int, 0, 0 >,  // SizesHolder
      std::index_sequence< 0, 1 >,           // Permutation
      DeviceType >;                          // Device

   SensorInterpolation()
   : sensorPositions( 0 )
   {}

   SensorInterpolation( GlobalIndexType numberOfSavedSteps, std::vector< VectorType >& sensorsPoints )
   : sensorPositions( sensorsPoints ),
     numberOfSensors( sensorsPoints.size() ),
     numberOfSavedSteps( numberOfSavedSteps )
   {
      sensors.setSizes( numberOfSavedSteps, numberOfSensors );
   }

   void
   init( std::vector< VectorType >& points, const int numberOfSensors, const int numberOfSavedSteps, bool includeBoundary )
   {
      this->numberOfSensors = numberOfSensors;
      sensors.setSizes( 2 * numberOfSavedSteps, numberOfSensors );
      sensorPositions.setSize( numberOfSensors );
      sensorPositions = points;
      this->numberOfSavedSteps = numberOfSavedSteps + 1;
      this->includeBoundary = includeBoundary;
   }

   template< typename SPHKernelFunction, typename EOS, typename SPHState >
   void
   interpolate( FluidPointer& fluid, BoundaryPointer& boundary, SPHState& sphState );

   template< typename SPHKernelFunction, typename EOS, typename ModelParams, typename LevelView >
   void
   interpolate( const std::vector< LevelView >& levels, ModelParams& modelParams );

   void
   save( const std::string outputFileName );

protected:
   GlobalIndexType numberOfSensors;
   SensorsDataArray sensors;
   VectorArrayType sensorPositions;
   bool includeBoundary = false;

   GlobalIndexType numberOfSavedSteps;
   GlobalIndexType sensorIndexer = 0;

public:
   template< typename SPHKernelFunction, typename EOS, typename SourcePointer >
   void
   accumulatePressureSource(
      const SourcePointer& source,
      const IndexArrayType& owners,
      ScalarArrayType& pressure,
      ScalarArrayType& gamma,
      GlobalIndexType levelId,
      RealType h,
      RealType mass,
      const typename EOS::ParamsType& eosParams );
};

template< typename SPHConfig, typename SPHSimulation >
class SensorWaterLevel
{
public:
   using DeviceType = typename SPHConfig::DeviceType;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;

   using LocalIndexType = typename SPHTraitsType::LocalIndexType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using CoordinatesType = typename SPHTraitsType::CoordinatesType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;
   using IndexArrayType = typename SPHTraitsType::IndexArrayType;
   using ScalarArrayType = typename SPHTraitsType::ScalarArrayType;

   using FluidPointer = typename SPHSimulation::FluidPointer;
   using BoundaryPointer = typename SPHSimulation::BoundaryPointer;

   /**
    * \brief Particle type is required to search through particle seats we want to
    * monitor with measuretool. Measuretool works with CLL search, so in case we use
    * CLLwL, for measuretool we need to load the base class.
    */
   //using ParticlesType = typename SPHSimulation::ParticlesType;
   using ParticlesType = std::conditional_t<
      SPHSimulation::ParticlesType::specifySearchedSetExplicitly(),
      typename SPHSimulation::ParticlesType::BaseType,
      typename SPHSimulation::ParticlesType >;

   using SensorsDataArray = Containers::NDArray<
      RealType,                              // Value
      Containers::SizesHolder< int, 0, 0 >,  // SizesHolder
      std::index_sequence< 0, 1 >,           // Permutation
      Devices::Host >;                       // Device - store data on HOST

   SensorWaterLevel()
   : sensorPositions( 0 )
   {}

   SensorWaterLevel(
      GlobalIndexType numberOfSavedSteps,
      std::vector< VectorType >& sensorsPoints,
      RealType levelIncrement,
      VectorType direction,
      RealType startLevel,
      RealType endLevel )
   : sensorPositions( sensorsPoints ),
     numberOfSensors( sensorsPoints.size() ),
     numberOfSavedSteps( numberOfSavedSteps ),
     levelIncrement( levelIncrement ),
     direction( direction )
   {
      sensors.setSizes( numberOfSavedSteps, numberOfSensors );

      //New
      numberOfLevels = TNL::ceil( ( endLevel - startLevel ) / levelIncrement );
      levels.setSize( numberOfLevels );
   }

   void
   init(
      std::vector< VectorType >& points,
      const int numberOfSavedSteps,
      RealType levelIncrement,
      VectorType direction,
      RealType startLevel,
      RealType endLevel )
   {
      this->numberOfSensors = points.size();
      sensors.setSizes( 2 * numberOfSavedSteps, numberOfSensors );
      sensorPositions.setSize( numberOfSensors );
      sensorPositions = points;

      numberOfLevels = TNL::ceil( ( endLevel - startLevel ) / levelIncrement );
      levels.setSize( numberOfLevels );
      this->direction = direction;
      this->levelIncrement = levelIncrement;
      this->numberOfSavedSteps = numberOfSavedSteps + 1;

      this->startLevel = startLevel;
      this->endLevel = endLevel;
   }

   [[nodiscard]] const SensorsDataArray&
   getSensorData() const;

   template< typename SPHKernelFunction, typename EOS, typename SPHState >
   void
   interpolate( FluidPointer& fluid, BoundaryPointer& boundary, SPHState& sphState );

   template< typename SPHKernelFunction, typename EOS, typename ModelParams, typename LevelView >
   void
   interpolate( const std::vector< LevelView >& levels, ModelParams& modelParams );

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

public:
   template< typename SPHKernelFunction, typename SourcePointer >
   void
   accumulateWaterSource(
      const SourcePointer& source,
      const VectorArrayType& samplePositions,
      const IndexArrayType& sampleOwners,
      ScalarArrayType& sampleGamma,
      GlobalIndexType levelId,
      RealType h,
      RealType mass );
};

}  // namespace TNL::SPH

#include <SPH/shared/measurement/MeasureTool.hpp>

#include <SPH/ParseConfigFile.h>
