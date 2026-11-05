// SPDX-FileComment: This file is part of TNL-SPH (https://gitlab.com/tnl-project/tnl-sph)
// SPDX-License-Identifier: MIT

#pragma once

#include <TNL/Logger.h>
#include <TNL/Config/ConfigDescription.h>
#include <TNL/Config/ParameterContainer.h>
#include <TNL/Functional.h>

#include <SPH/SPHTraits.h>
#include <SPH/ParseConfigFile.h>
#include <SPH/shared/measurement/MeasureTool.h>
#include <SPH/shared/measurement/measureVolumetricFlowRate.h>

namespace TNL::SPH {

template< typename SimulationType >
class SimulationMonitor
{
public:
   using SPHConfig = typename SimulationType::ModelType::SPHConfig;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;

   using LocalIndexType = typename SPHTraitsType::LocalIndexType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using CoordinatesType = typename SPHTraitsType::CoordinatesType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;

   using GridInterpolation = TNL::SPH::InterpolateToGrid< SPHConfig, SimulationType >;
   using SensorPressure = TNL::SPH::SensorInterpolation< SPHConfig, SimulationType >;
   using SensorWaterLevel = TNL::SPH::SensorWaterLevel< SPHConfig, SimulationType >;

   using FluidPointer = typename SimulationType::FluidPointer;
   using BoundaryPointer = typename SimulationType::BoundaryPointer;
   using ModelParams = typename SimulationType::ModelParams;
   using TimeStepping = typename SimulationType::TimeStepping;

   using ParticlesConfig = typename SimulationType::ParticlesType::Config;
   using SPHDefs = typename SimulationType::ModelType::ModelConfigType;
   using FlowRateMonitorType = TNL::SPH::MeasureVolumetricFlowRate< ParticlesConfig, SPHDefs >;

   // Config values stored separately because zone setup is deferred to
   // setupVolumetricFlowRateZones() which needs the fluid pointer.
   struct VolumetricFlowRateEntry
   {
      FlowRateMonitorType monitor;
      VectorType pointMin;
      VectorType pointMax;
      VectorType orientation;
      RealType averagingInterval = 0.1;
      std::string outputPath;
   };

   void
   init( TNL::Config::ParameterContainer& parameters, TimeStepping& timeStepping, TNL::Logger& logger )
   {
      TNL::Config::ConfigDescription measuretoolConfig;
      outputDirecotry = parameters.getParameter< std::string >( "outputDirectory" );

      this->numberOfInterpolationPlanes = parameters.getParameter< int >( "interpolationPlanesCount" );
      configSetupInterpolation( measuretoolConfig, parameters );
      this->numberOfPressureSensors = parameters.getParameter< int >( "pressureSensorsCount" );
      configSetupPressureSensor( measuretoolConfig, parameters );
      this->numberOfWaterLevelSensors = parameters.getParameter< int >( "waterLevelSensorsCount" );
      configWaterLevelSensor( measuretoolConfig, parameters );
      this->numberOfVolumetricFlowRatePlanes = parameters.getParameter< int >( "volumetricFlowRatePlanesCount" );
      configSetupVolumetricFlowRate( measuretoolConfig, parameters );

      // parse parameters
      TNL::Config::ParameterContainer measuretoolParameters;
      const std::string configPath = parameters.getParameter< std::string >( "measuretoolConfig" );
      if( ! configPath.empty() ) {
         logger.writeParameter( "Reading measuretool config:", configPath );
         parseMeasuretoolConfig( measuretoolConfig, measuretoolParameters, configPath );
      }
      else if( parameters.checkParameter( "measuretool" ) ) {
         const std::string inlineConfig = parameters.getParameter< std::string >( "measuretool" );
         if( ! inlineConfig.empty() ) {
            logger.writeParameter( "Reading measuretool config:", "inline (config.jsonc subsection)" );
            try {
               measuretoolParameters = TNL::SPH::parseJSONConfigString( inlineConfig, "inline measuretool", measuretoolConfig );
            }
            catch( const std::exception& e ) {
               std::cerr << "Failed to parse the inline measuretool configuration due to the following error:\n"
                         << e.what() << std::endl;
            }
         }
      }
      initializeInterpolations( parameters, measuretoolParameters, timeStepping, logger );
      initializePressureSensors( parameters, measuretoolParameters, timeStepping, logger );
      initializeWaterLevelSensors( parameters, measuretoolParameters, timeStepping, logger );
      initializeVolumetricFlowRatePlanes( parameters, measuretoolParameters, timeStepping, logger );
   }

   void
   configSetupInterpolation( TNL::Config::ConfigDescription& measuretoolConfig, TNL::Config::ParameterContainer& parameters )
   {
      for( int i = 0; i < this->numberOfInterpolationPlanes; i++ ) {
         std::string prefix = "plane" + std::to_string( i + 1 );
         measuretoolConfig.addEntry< std::string >( prefix + "Identifier", "Identifier of the interpolation plane.", prefix );
         measuretoolConfig.addEntry< RealType >(
            prefix + "OutputTime",
            "Define interval of pressure measurements. ",
            parameters.getParameter< RealType >( "snapshotPeriod" ) );

         measuretoolConfig.addRequiredEntry< RealType >(
            prefix + "GridOriginX", "Define interpolation plane origin x componenet." );
         measuretoolConfig.addRequiredEntry< RealType >(
            prefix + "GridOriginY", "Define interpolation plane origin y componenet." );
         if constexpr( SPHConfig::spaceDimension == 3 )
            measuretoolConfig.addRequiredEntry< RealType >(
               prefix + "GridOriginZ", "Define interpolation plane origin z componenet." );

         measuretoolConfig.addRequiredEntry< int >( prefix + "GridSizeX", "Define interpolation plane size in x direction." );
         measuretoolConfig.addRequiredEntry< int >( prefix + "GridSizeY", "Define interpolation plane size in y direction." );
         if constexpr( SPHConfig::spaceDimension == 3 )
            measuretoolConfig.addRequiredEntry< int >(
               prefix + "GridSizeZ", "Define interpolation plane size in z direction." );

         measuretoolConfig.addRequiredEntry< RealType >( prefix + "GridStepX", "Define grid step size in x direction." );
         measuretoolConfig.addRequiredEntry< RealType >( prefix + "GridStepY", "Define grid step size in y direction." );
         if constexpr( SPHConfig::spaceDimension == 3 )
            measuretoolConfig.addRequiredEntry< RealType >( prefix + "GridStepZ", "Define grid step size in z direction." );
      }
   }

   void
   configSetupPressureSensor( TNL::Config::ConfigDescription& measuretoolConfig, TNL::Config::ParameterContainer& parameters )
   {
      measuretoolConfig.addEntry< RealType >(
         "sensorPOutputTime",
         "Define interval of pressure measurements. ",
         parameters.getParameter< RealType >( "snapshotPeriod" ) );
      measuretoolConfig.addEntry< bool >( "sensorPIncludeBoundary", "Include boundary data to sensor measurement. ", false );
      for( int i = 0; i < this->numberOfPressureSensors; i++ ) {
         std::string prefix = "sensorPPoint" + std::to_string( i + 1 );
         measuretoolConfig.addRequiredEntry< RealType >(
            prefix + "X", "Define pressure sensor point " + std::to_string( i + 1 ) + " x component." );
         measuretoolConfig.addRequiredEntry< RealType >(
            prefix + "Y", "Define pressure sensor point " + std::to_string( i + 1 ) + " y component." );
         if constexpr( SPHConfig::spaceDimension == 3 )
            measuretoolConfig.addRequiredEntry< RealType >(
               prefix + "Z", "Define pressure sensor point " + std::to_string( i + 1 ) + " z component." );
      }
   }

   void
   configWaterLevelSensor( TNL::Config::ConfigDescription& measuretoolConfig, TNL::Config::ParameterContainer& parameters )
   {
      measuretoolConfig.addEntry< RealType >(
         "sensorWlOutputTime",
         "Define interval of pressure measurements. ",
         parameters.getParameter< RealType >( "snapshotPeriod" ) );
      measuretoolConfig.addEntry< RealType >( "sensorWlStartLevel", "Define the starting point to measure water level. ", 0 );
      measuretoolConfig.addEntry< RealType >( "sensorWlEndLevel", "Define the end point to measure water level. ", 0 );
      measuretoolConfig.addEntry< RealType >(
         "sensorWlLevelIncrement",
         "Define water level measurement increment. ",
         parameters.getParameter< RealType >( "snapshotPeriod" ) );
      measuretoolConfig.addEntry< RealType >(
         "sensorWlDirectionX", "Dirction in which the water level measurement is performed, x component. ", 0 );
      measuretoolConfig.addEntry< RealType >(
         "sensorWlDirectionY", "Dirction in which the water level measurement is performed, y component. ", 0 );
      if constexpr( SPHConfig::spaceDimension == 3 )
         measuretoolConfig.addEntry< RealType >(
            "sensorWlDirectionZ", "Dirction in which the water level measurement is performed, y component. ", 0 );
      for( int i = 0; i < this->numberOfWaterLevelSensors; i++ ) {
         std::string prefix = "sensorWlPoint" + std::to_string( i + 1 );
         measuretoolConfig.addRequiredEntry< RealType >(
            prefix + "X", "Define water level sensor point " + std::to_string( i + 1 ) + " x component." );
         measuretoolConfig.addRequiredEntry< RealType >(
            prefix + "Y", "Define water level sensor point " + std::to_string( i + 1 ) + " y component." );
         if constexpr( SPHConfig::spaceDimension == 3 )
            measuretoolConfig.addRequiredEntry< RealType >(
               prefix + "Z", "Define water level sensor point " + std::to_string( i + 1 ) + " z component." );
      }
   }

   void
   configSetupVolumetricFlowRate(
      TNL::Config::ConfigDescription& measuretoolConfig,
      TNL::Config::ParameterContainer& parameters )
   {
      for( int i = 0; i < this->numberOfVolumetricFlowRatePlanes; i++ ) {
         std::string prefix = "vfrPlane" + std::to_string( i + 1 );
         measuretoolConfig.addEntry< std::string >(
            prefix + "Identifier", "Identifier of the volumetric flow rate plane.", prefix );
         measuretoolConfig.addEntry< RealType >(
            prefix + "AveragingInterval", "Averaging interval for the volumetric flow rate.", 0.1 );
         measuretoolConfig.addRequiredEntry< RealType >(
            prefix + "PointMinX", "Define the first point of the measurement plane, x component." );
         measuretoolConfig.addRequiredEntry< RealType >(
            prefix + "PointMinY", "Define the first point of the measurement plane, y component." );
         if constexpr( SPHConfig::spaceDimension == 3 )
            measuretoolConfig.addRequiredEntry< RealType >(
               prefix + "PointMinZ", "Define the first point of the measurement plane, z component." );
         measuretoolConfig.addRequiredEntry< RealType >(
            prefix + "PointMaxX", "Define the second point of the measurement plane, x component." );
         measuretoolConfig.addRequiredEntry< RealType >(
            prefix + "PointMaxY", "Define the second point of the measurement plane, y component." );
         if constexpr( SPHConfig::spaceDimension == 3 )
            measuretoolConfig.addRequiredEntry< RealType >(
               prefix + "PointMaxZ", "Define the second point of the measurement plane, z component." );
         measuretoolConfig.addRequiredEntry< RealType >(
            prefix + "OrientationX", "Define the orientation of the measurement plane, x component." );
         measuretoolConfig.addRequiredEntry< RealType >(
            prefix + "OrientationY", "Define the orientation of the measurement plane, y component." );
         if constexpr( SPHConfig::spaceDimension == 3 )
            measuretoolConfig.addRequiredEntry< RealType >(
               prefix + "OrientationZ", "Define the orientation of the measurement plane, z component." );
      }
   }

   void
   parseMeasuretoolConfig(
      TNL::Config::ConfigDescription& measuretoolConfig,
      TNL::Config::ParameterContainer& measuretoolParameters,
      const std::string& configPath )
   {
      try {
         measuretoolParameters = TNL::SPH::parseConfigFile( configPath, measuretoolConfig );
      }
      catch( const std::exception& e ) {
         std::cerr << "Failed to parse the measuretool configuration file " << configPath << " due to the following error:\n"
                   << e.what() << std::endl;
      }
      catch( ... ) {
         std::cerr << "Failed to parse the measuretool configuration file " << configPath << " due to an unknown C++ exception."
                   << std::endl;
         throw;
      }
   }

   void
   initializeInterpolations(
      TNL::Config::ParameterContainer& parameters,
      TNL::Config::ParameterContainer& measuretoolParameters,
      TimeStepping& timeStepping,
      TNL::Logger& logger )
   {
      const RealType simulationEndTime = parameters.getParameter< RealType >( "finalTime" );

      for( int i = 0; i < this->numberOfInterpolationPlanes; i++ ) {
         std::string prefix = "plane" + std::to_string( i + 1 );
         std::string interpolationIdentifier = measuretoolParameters.getParameter< std::string >( prefix + "Identifier" );
         interpolations.insert( { interpolationIdentifier, GridInterpolation() } );
         interpolations[ interpolationIdentifier ].init( measuretoolParameters, prefix );
         const float interpolationPeriod = measuretoolParameters.getParameter< float >( prefix + "OutputTime" );
         timeStepping.addOutputTimer( "interpolate", interpolationPeriod );
         logger.writeParameter( "Initialized interpolation plane:", interpolationIdentifier );
      }
   }

   void
   initializePressureSensors(
      TNL::Config::ParameterContainer& parameters,
      TNL::Config::ParameterContainer& measuretoolParameters,
      TimeStepping& timeStepping,
      TNL::Logger& logger )
   {
      const RealType simulationEndTime = parameters.getParameter< RealType >( "finalTime" );

      const float pressureSensorPeriod = measuretoolParameters.getParameter< float >( "sensorPOutputTime" );
      const float numberOfSavedStepsPressure = simulationEndTime / pressureSensorPeriod;
      std::vector< VectorType > sensorsPoints( this->numberOfPressureSensors );
      for( int i = 0; i < this->numberOfPressureSensors; i++ ) {
         std::string prefix = "sensorPPoint" + std::to_string( i + 1 );
         sensorsPoints[ i ] = getXyzCamelCase< VectorType >( measuretoolParameters, prefix );
      }
      const bool includeBoundary = measuretoolParameters.getParameter< bool >( "sensorPIncludeBoundary" );
      pressureSensors.init( sensorsPoints, this->numberOfPressureSensors, numberOfSavedStepsPressure, includeBoundary );
      timeStepping.addOutputTimer( "sensor_pressure", pressureSensorPeriod );
      logger.writeParameter( "Number of initialized pressure sensors:", this->numberOfPressureSensors );
   }

   void
   initializeWaterLevelSensors(
      TNL::Config::ParameterContainer& parameters,
      TNL::Config::ParameterContainer& measuretoolParameters,
      TimeStepping& timeStepping,
      TNL::Logger& logger )
   {
      const RealType simulationEndTime = parameters.getParameter< RealType >( "finalTime" );

      const float waterLevelSensorPeriod = measuretoolParameters.getParameter< float >( "sensorWlOutputTime" );
      const float numberOfSavedStepsWaterLevel = simulationEndTime / waterLevelSensorPeriod;
      const RealType levelIncrement = measuretoolParameters.getParameter< RealType >( "sensorWlLevelIncrement" );
      const RealType startWaterLevel = measuretoolParameters.getParameter< RealType >( "sensorWlStartLevel" );
      const RealType endWaterLevel = measuretoolParameters.getParameter< RealType >( "sensorWlEndLevel" );
      const VectorType measureInDirection = getXyzCamelCase< VectorType >( measuretoolParameters, "sensorWlDirection" );
      std::vector< VectorType > waterLevelSensorPoints( this->numberOfWaterLevelSensors );
      for( int i = 0; i < this->numberOfWaterLevelSensors; i++ ) {
         std::string prefix = "sensorWlPoint" + std::to_string( i + 1 );
         waterLevelSensorPoints[ i ] = getXyzCamelCase< VectorType >( measuretoolParameters, prefix );
      }
      waterLevelSensors.init(
         waterLevelSensorPoints,
         numberOfSavedStepsWaterLevel,
         levelIncrement,
         measureInDirection,
         startWaterLevel,
         endWaterLevel );
      timeStepping.addOutputTimer( "sensor_waterLevel", waterLevelSensorPeriod );
      logger.writeParameter( "Number of initialized waterLevelSensors sensors:", this->numberOfWaterLevelSensors );
   }

   void
   initializeVolumetricFlowRatePlanes(
      TNL::Config::ParameterContainer& parameters,
      TNL::Config::ParameterContainer& measuretoolParameters,
      TimeStepping& timeStepping,
      TNL::Logger& logger )
   {
      for( int i = 0; i < this->numberOfVolumetricFlowRatePlanes; i++ ) {
         std::string prefix = "vfrPlane" + std::to_string( i + 1 );
         std::string identifier = measuretoolParameters.getParameter< std::string >( prefix + "Identifier" );

         VolumetricFlowRateEntry entry;
         entry.pointMin = getXyzCamelCase< VectorType >( measuretoolParameters, prefix + "PointMin" );
         entry.pointMax = getXyzCamelCase< VectorType >( measuretoolParameters, prefix + "PointMax" );
         entry.orientation = getXyzCamelCase< VectorType >( measuretoolParameters, prefix + "Orientation" );
         entry.averagingInterval = measuretoolParameters.getParameter< RealType >( prefix + "AveragingInterval" );
         entry.monitor.averagingInterval = entry.averagingInterval;
         entry.outputPath = outputDirecotry + "/volumetricFlowRate_" + identifier + ".dat";

         volumetricFlowRateMonitors.insert( { identifier, std::move( entry ) } );
         logger.writeParameter( "Initialized volumetric flow rate plane:", identifier );
      }
   }

   void
   setupVolumetricFlowRateZones( FluidPointer& fluid )
   {
      for( auto& [ key, entry ] : volumetricFlowRateMonitors ) {
         entry.monitor.init( fluid, entry.pointMin, entry.pointMax, entry.orientation );
      }
   }

   void
   measureVolumetricFlowRatesStep( FluidPointer& fluid, ModelParams& modelParams, TimeStepping& timeStepping )
   {
      const RealType dt = timeStepping.getTimeStep();
      for( auto& [ key, entry ] : volumetricFlowRateMonitors ) {
         entry.monitor.measureVolumetricFlowRate( fluid, modelParams, dt );
      }
   }

   void
   updateAndOutputVolumetricFlowRates( TimeStepping& timeStepping )
   {
      const RealType dt = timeStepping.getTimeStep();
      const typename FlowRateMonitorType::IndexType step = timeStepping.getStep();
      const RealType time = timeStepping.getTime();
      for( auto& [ key, entry ] : volumetricFlowRateMonitors ) {
         entry.monitor.updateVolumetricFlowRate( dt );
         entry.monitor.output( entry.outputPath, step, time );
      }
   }

   void
   measureVolumetricFlowRates( FluidPointer& fluid, ModelParams& modelParams, TimeStepping& timeStepping )
   {
      measureVolumetricFlowRatesStep( fluid, modelParams, timeStepping );
      updateAndOutputVolumetricFlowRates( timeStepping );
   }

   template< typename KernelFunction, typename EOS >
   void
   measure(
      FluidPointer& fluid,
      BoundaryPointer& boundary,
      ModelParams& modelParams,
      TimeStepping& timeStepping,
      TNL::Logger& logger,
      const std::string& verbose )
   {
      if( timeStepping.checkOutputTimer( "interpolate" ) ) {
         for( auto& [ key, val ] : this->interpolations ) {
            val.template interpolate< KernelFunction >( fluid, boundary, modelParams );
            std::string outputFileNameInterpolation =
               outputDirecotry + "/" + key + "_" + std::to_string( timeStepping.getTime() ) + "s_interpolation.vtk";
            val.save( outputFileNameInterpolation );
            logger.writeParameter( "Saved:", outputFileNameInterpolation );
            if( verbose == "full" )
               logger.writeParameter( "Interpolation -" + key + ":", "Done." );
         }
      }
      if( timeStepping.checkOutputTimer( "sensor_pressure" ) ) {
         if( this->numberOfPressureSensors > 0 ) {
            pressureSensors.template interpolate< KernelFunction, EOS >( fluid, boundary, modelParams );
            if( verbose == "full" )
               logger.writeParameter( "Pressure sensor measurement:", "Done." );
         }
      }

      if( timeStepping.checkOutputTimer( "sensor_waterLevel" ) ) {
         if( this->numberOfWaterLevelSensors > 0 ) {
            waterLevelSensors.template interpolate< KernelFunction, EOS >( fluid, boundary, modelParams );
            if( verbose == "full" )
               logger.writeParameter( "Water level measurement:", "Done." );
         }
      }
   }

   void
   save( TNL::Logger& logger )
   {
      if( this->numberOfPressureSensors > 0 ) {
         const std::string sensorsOutputFilename = outputDirecotry + "/sensorsPressure.dat";
         pressureSensors.save( sensorsOutputFilename );
         logger.writeParameter( "Saved:", sensorsOutputFilename );
      }
      if( this->numberOfWaterLevelSensors > 0 ) {
         const std::string sensorsOutputFilename = outputDirecotry + "/sensorsWaterLevel.dat";
         waterLevelSensors.save( sensorsOutputFilename );
         logger.writeParameter( "Saved:", sensorsOutputFilename );
      }
   }

   std::map< std::string, VolumetricFlowRateEntry > volumetricFlowRateMonitors;

protected:
   int numberOfInterpolationPlanes = 0;
   std::map< std::string, GridInterpolation > interpolations;

   //int numberOfPressureSensorsSets = 0;
   //std::map< std::string, SensorPressure > pressureSensors;
   int numberOfPressureSensors = 0;
   SensorPressure pressureSensors;

   //int numberOfWaterLevelSensorsSets = 0;
   //std::map< std::string, SensorWaterLevel > waterLevelSensors;
   int numberOfWaterLevelSensors = 0;
   SensorWaterLevel waterLevelSensors;

   //std::map< std::string, std::unique_ptr< Sensor > > sensors;
   //Sensor* sensor = nullptr;
   int numberOfVolumetricFlowRatePlanes = 0;
   //SensorsVolumetricFlowRates volumetricFlowRatePlanes;

   std::string outputDirecotry;
};

}  // namespace TNL::SPH
