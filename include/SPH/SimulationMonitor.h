#include <TNL/Config/parseINIConfigFile.h>
#include <TNL/Logger.h>


#include "SPHTraits.h"
#include "TNL/Config/ConfigDescription.h"
#include "TNL/Config/ParameterContainer.h"
#include "TNL/Functional.h"
#include "shared/Measuretool.h"

namespace TNL {
namespace SPH {

template< typename SimulationType >
class SimulationMonitor
{
public:

   using SPHConfig = typename SimulationType::ModelType::SPHConfig;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;

   using LocalIndexType = typename SPHTraitsType::LocalIndexType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using IndexVectorType = typename SPHTraitsType::IndexVectorType;
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

   //SimulationMonitor() = default;

   void
   init( TNL::Config::ParameterContainer& parameters, TimeStepping& timeStepping, TNL::Logger& logger )
   {
      TNL::Config::ConfigDescription measuretoolConfig;
      outputDirecotry = parameters.getParameter< std::string >( "output-directory" );

      this->numberOfInterpolationPlanes = parameters.getParameter< int >( "interpolation-planes-count" );
      configSetupInterpolation( measuretoolConfig, parameters );

      this->numberOfPressureSensors = parameters.getParameter< int >( "pressure-sensors-count" );
      configSetupPressureSensor( measuretoolConfig, parameters );

      this->numberOfWaterLevelSensors = parameters.getParameter< int >( "water-level-sensors-count" );
      configWaterLevelSensor( measuretoolConfig, parameters );

      // parse parameters
      TNL::Config::ParameterContainer measuretoolParameters;
      const std::string configPath = parameters.getParameter< std::string >( "measuretool-config" );
      logger.writeParameter( "Reading measuretool config:", configPath );
      parseMeasuretoolConfig( measuretoolConfig, measuretoolParameters, configPath );

      initializeInterpolations( parameters, measuretoolParameters, timeStepping, logger );

      initializePressureSensors( parameters, measuretoolParameters, timeStepping, logger );

      initializeWaterLevelSensors( parameters, measuretoolParameters, timeStepping, logger );
   }

   void
   configSetupInterpolation( TNL::Config::ConfigDescription& measuretoolConfig, TNL::Config::ParameterContainer& parameters )
   {
      for( int i = 0; i < this->numberOfInterpolationPlanes; i++ )
      {
         std::string prefix = "plane-" + std::to_string( i + 1 ) + "-";
         measuretoolConfig.addEntry< std::string >( prefix + "identifier", "Identifier of the interpolation plane.", prefix );
         measuretoolConfig.addEntry< RealType >( prefix + "outputTime", "Define interval of pressure measurements. ", parameters.getParameter< RealType >( "snapshot-period" ) );

         measuretoolConfig.addRequiredEntry< RealType >( prefix + "gridOrigin-x", "Define interpolation plane origin x componenet." );
         measuretoolConfig.addRequiredEntry< RealType >( prefix + "gridOrigin-y", "Define interpolation plane origin y componenet." );
         if constexpr( SPHConfig::spaceDimension == 3 )
            measuretoolConfig.addRequiredEntry< RealType >( prefix + "gridOrigin-z", "Define interpolation plane origin z componenet." );

         measuretoolConfig.addRequiredEntry< int >( prefix + "gridSize-x", "Define interpolation plane size in x direction." );
         measuretoolConfig.addRequiredEntry< int >( prefix + "gridSize-y", "Define interpolation plane size in y direction." );
         if constexpr( SPHConfig::spaceDimension == 3 )
            measuretoolConfig.addRequiredEntry< int >( prefix + "gridSize-z", "Define interpolation plane size in z direction." );

         measuretoolConfig.addRequiredEntry< RealType >( prefix + "gridStep-x", "Define grid step size in x direction." );
         measuretoolConfig.addRequiredEntry< RealType >( prefix + "gridStep-y", "Define grid step size in y direction." );
         if constexpr( SPHConfig::spaceDimension == 3 )
            measuretoolConfig.addRequiredEntry< RealType >( prefix + "gridStep-z", "Define grid step size in z direction." );
      }
   }

   void
   configSetupPressureSensor( TNL::Config::ConfigDescription& measuretoolConfig, TNL::Config::ParameterContainer& parameters )
   {
      measuretoolConfig.addEntry< RealType >( "sensor-p-outputTime", "Define interval of pressure measurements. ", parameters.getParameter< RealType >( "snapshot-period" ) );
      measuretoolConfig.addEntry< bool >( "sensor-p-includeBoundary", "Include boundary data to sensor measurement. ", false );
      for( int i = 0; i < this->numberOfPressureSensors; i++ )
      {
         std::string prefix = "sensor-p-point-" + std::to_string( i + 1 );
         measuretoolConfig.addRequiredEntry< RealType >( prefix + "-x", "Define pressure sensor point " + std::to_string( i + 1 ) + " x component." );
         measuretoolConfig.addRequiredEntry< RealType >( prefix + "-y", "Define pressure sensor point " + std::to_string( i + 1 ) + " y component." );
         if constexpr( SPHConfig::spaceDimension == 3 )
            measuretoolConfig.addRequiredEntry< RealType >( prefix + "-z", "Define pressure sensor point " + std::to_string( i + 1 ) + " z component." );
      }
   }

   void
   configWaterLevelSensor( TNL::Config::ConfigDescription& measuretoolConfig, TNL::Config::ParameterContainer& parameters )
   {
      measuretoolConfig.addEntry< RealType >( "sensor-wl-outputTime", "Define interval of pressure measurements. ", parameters.getParameter< RealType >( "snapshot-period" ) );
      measuretoolConfig.addEntry< RealType >( "sensor-wl-start-level", "Define the starting point to measure water level. ", 0 );
      measuretoolConfig.addEntry< RealType >( "sensor-wl-end-level", "Define the end point to measure water level. ", 0 );
      measuretoolConfig.addEntry< RealType >( "sensor-wl-level-increment", "Define water level measurement increment. ", parameters.getParameter< RealType >( "snapshot-period" ) );
      measuretoolConfig.addEntry< RealType >( "sensor-wl-direction-x", "Dirction in which the water level measurement is performed, x component. ", 0 );
      measuretoolConfig.addEntry< RealType >( "sensor-wl-direction-y", "Dirction in which the water level measurement is performed, y component. ", 0 );
      if constexpr( SPHConfig::spaceDimension == 3 )
         measuretoolConfig.addEntry< RealType >( "sensor-wl-direction-z", "Dirction in which the water level measurement is performed, y component. ", 0 );
      for( int i = 0; i < this->numberOfWaterLevelSensors; i++ )
      {
         std::string prefix = "sensor-wl-point-" + std::to_string( i + 1 );
         measuretoolConfig.addRequiredEntry< RealType >( prefix + "-x", "Define water level sensor point " + std::to_string( i + 1 ) + " x component." );
         measuretoolConfig.addRequiredEntry< RealType >( prefix + "-y", "Define water level sensor point " + std::to_string( i + 1 ) + " y component." );
         if constexpr( SPHConfig::spaceDimension == 3 )
            measuretoolConfig.addRequiredEntry< RealType >( prefix + "-z", "Define water level sensor point " + std::to_string( i + 1 ) + " z component." );
      }
   }

   void
   parseMeasuretoolConfig( TNL::Config::ConfigDescription& measuretoolConfig,
                           TNL::Config::ParameterContainer& measuretoolParameters,
                           const std::string& configPath )
   {
      try {
          measuretoolParameters = TNL::Config::parseINIConfigFile( configPath, measuretoolConfig );
      }
      catch ( const std::exception& e ) {
          std::cerr << "Failed to parse the measuretool configuration file " << configPath << " due to the following error:\n" << e.what() << std::endl;
      }
      catch (...) {
          std::cerr << "Failed to parse the measuretool configuration file " << configPath << " due to an unknown C++ exception." << std::endl;
          throw;
      }

   }

   void
   initializeInterpolations( TNL::Config::ParameterContainer& parameters,
                             TNL::Config::ParameterContainer& measuretoolParameters,
                             TimeStepping& timeStepping,
                             TNL::Logger& logger )
   {
      const RealType simulationEndTime = parameters.getParameter< RealType >( "final-time" );

      for( int i = 0; i < this->numberOfInterpolationPlanes; i++ )
      {
         std::string prefix = "plane-" + std::to_string( i + 1 ) + "-";
         std::string interpolationIdentifier = measuretoolParameters.getParameter< std::string >( prefix  + "identifier" );
         interpolations.insert( {interpolationIdentifier, GridInterpolation() });
         interpolations[ interpolationIdentifier ].init( measuretoolParameters, prefix );
         const float interpolationPeriod = measuretoolParameters.getParameter< float >( prefix + "outputTime" );
         timeStepping.addOutputTimer( "interpolate", interpolationPeriod );
         logger.writeParameter( "Initialized interpolation plane:", interpolationIdentifier );
      }

   }

   void
   initializePressureSensors( TNL::Config::ParameterContainer& parameters,
                              TNL::Config::ParameterContainer& measuretoolParameters,
                              TimeStepping& timeStepping,
                              TNL::Logger& logger )
   {
      const RealType simulationEndTime = parameters.getParameter< RealType >( "final-time" );

      const float pressureSensorPeriod = measuretoolParameters.getParameter< float >( "sensor-p-outputTime" );
      const float numberOfSavedStepsPressure = simulationEndTime / pressureSensorPeriod;
      std::vector< VectorType > sensorsPoints( this->numberOfPressureSensors );
      for( int i = 0; i < this->numberOfPressureSensors; i++ ){
         std::string prefix = "sensor-p-point-" + std::to_string( i + 1 );
         sensorsPoints[ i ] = measuretoolParameters.getXyz< VectorType >( prefix );
      }
      const bool includeBoundary = measuretoolParameters.getParameter< bool >( "sensor-p-includeBoundary" );
      pressureSensors.init( sensorsPoints, this->numberOfPressureSensors, numberOfSavedStepsPressure, includeBoundary );
      timeStepping.addOutputTimer( "sensor_pressure", pressureSensorPeriod );
      logger.writeParameter( "Number of initialized pressure sensors:", this->numberOfPressureSensors );
   }

   void
   initializeWaterLevelSensors( TNL::Config::ParameterContainer& parameters,
                                TNL::Config::ParameterContainer& measuretoolParameters,
                                TimeStepping& timeStepping,
                                TNL::Logger& logger )
   {
      const RealType simulationEndTime = parameters.getParameter< RealType >( "final-time" );

      const float waterLevelSensorPeriod = measuretoolParameters.getParameter< float >( "sensor-wl-outputTime" );
      const float numberOfSavedStepsWaterLevel = simulationEndTime / waterLevelSensorPeriod;
      const RealType levelIncrement = measuretoolParameters.getParameter< RealType >( "sensor-wl-level-increment" );
      const RealType startWaterLevel =  measuretoolParameters.getParameter< RealType >( "sensor-wl-start-level" );
      const RealType endWaterLevel =  measuretoolParameters.getParameter< RealType >( "sensor-wl-end-level" );
      const VectorType measureInDirection =  measuretoolParameters.getXyz< VectorType >( "sensor-wl-direction" );
      std::vector< VectorType > waterLevelSensorPoints( this->numberOfWaterLevelSensors );
      for( int i = 0; i < this->numberOfWaterLevelSensors; i++ ){
         std::string prefix = "sensor-wl-point-" + std::to_string( i + 1 );
         waterLevelSensorPoints[ i ] = measuretoolParameters.getXyz< VectorType >( prefix );
      }
      waterLevelSensors.init( waterLevelSensorPoints,
                              numberOfSavedStepsWaterLevel,
                              levelIncrement,
                              measureInDirection,
                              startWaterLevel,
                              endWaterLevel );
      timeStepping.addOutputTimer( "sensor_waterLevel", waterLevelSensorPeriod );
      logger.writeParameter( "Number of initialized waterLevelSensors sensors:", this->numberOfWaterLevelSensors );

   }

   template< typename KernelFunction, typename EOS >
   void
   measure( FluidPointer& fluid,
            BoundaryPointer& boundary,
            ModelParams& modelParams,
            TimeStepping& timeStepping,
            TNL::Logger& logger,
            const std::string& verbose )
   {
      if( timeStepping.checkOutputTimer( "interpolate" ) )
      {
         for ( auto& [ key, val ] : this->interpolations )
         {
            val.template interpolate< KernelFunction >( fluid, boundary, modelParams );
            std::string outputFileNameInterpolation = outputDirecotry + "/" + key + "_" +
                                                      //std::to_string( timeStepping.getStep() ) + "_interpolation.vtk";
                                                      std::to_string( timeStepping.getTime() ) + "s_interpolation.vtk";
            val.save( outputFileNameInterpolation );
            logger.writeParameter( "Saved:", outputFileNameInterpolation );
            if( verbose == "full" )
               logger.writeParameter( "Interpolation -" + key + ":", "Done." );
         }
      }
      if( timeStepping.checkOutputTimer( "sensor_pressure" ) )
      {
         if( this->numberOfPressureSensors > 0 ){
            pressureSensors.template interpolate< KernelFunction, EOS >( fluid, boundary, modelParams );
            if( verbose == "full" )
               logger.writeParameter( "Pressure sensor measurement:", "Done." );
         }
      }

      if( timeStepping.checkOutputTimer( "sensor_waterLevel" ) )
      {
         if( this->numberOfWaterLevelSensors > 0 ){
            waterLevelSensors.template interpolate< KernelFunction, EOS >( fluid, boundary, modelParams );
            if( verbose == "full" )
               logger.writeParameter( "Water level measurement:", "Done." );
         }
      }

   }

   void
   save( TNL::Logger& logger )
   {
      if( this->numberOfPressureSensors > 0 ){
         const std::string sensorsOutputFilename = outputDirecotry + "/sensorsPressure.dat";
         pressureSensors.save( sensorsOutputFilename );
         logger.writeParameter( "Saved:", sensorsOutputFilename );
      }
      if( this->numberOfWaterLevelSensors > 0 ){
         const std::string sensorsOutputFilename = outputDirecotry + "/sensorsWaterLevel.dat";
         waterLevelSensors.save( sensorsOutputFilename );
         logger.writeParameter( "Saved:", sensorsOutputFilename );
      }
   }


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

   std::string outputDirecotry;
};

} // SPH
} // TNL

