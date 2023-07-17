#include <iostream>
#include <fstream> //temp, to write output

/**
 *  Benchamrk stuff.
 */
#include <TNL/Benchmarks/Benchmarks.h>

/**
 * Particle system.
 */
#include "../../../Particles/ParticlesLinkedListFloating.h"

/**
 * Particle system reader.
 **/
#include "../../../Readers/VTKReader.h"
#include "../../../Writers/VTKWriter.h"
#include "../../../Readers/readSPHSimulation.h"

/**
 * Case configuration
 * One configuration for particle system, one for SPH.
 */
#include "sources/ParticlesConfig.h"
#include "sources/SPHCaseConfig.h"
#include "sources/MeasuretoolConfig.h"
#include "sources/SimulationControlConfig.h"


/**
 * SPH general toolds.
 */
#include "../../../SPH/SPH.h"

/**
 * SPH model.
 */
#include "../../../SPH/Models/WCSPH_DBC/Variables.h"
#include "../../../SPH/Models/WCSPH_DBC/Interactions.h"

#include "../../../SPH/Kernels.h" //TODO: Move to another.

/**
 * Time step control.
 */
#include "../../../SPH/TimeStep.h"

/**
 * Measuretool draft.
 */
#include "../../../SPH/shared/Measuretool.h"

using namespace TNL::ParticleSystem;

int main( int argc, char* argv[] )
{
   /**
    * Load simulation configs.
    * - Particle system config:
    *   config for definition of particle system (datatypes, dimension,...)
    *   config with parameters of particle system (domain size, search radius,...)
    *
    * - Configuration of particle system.
    *   config with initial parameters of the particle system
    *
    * - SPH method config:
    *   config with parameteres and constants of the SPH method
    *
    * - Simulation control:
    *   config with path to initial condition, path to store results, end time etc.
    */
   using SimulationControl = SPH::SimulationControlConfiguration::SPHSimulationControl;

   using SPHConfig = SPH::SPHConfig::SPHConfig< SimulationControl::DeviceType >;
   using SPHParams = SPH::SPHConfig::SPHParamsConfig< SPHConfig >;

   using ParticlesConfig = ParticleSystemConfig::ParticleSystemConfig< SimulationControl::DeviceType >;
   using ParticlesParams = ParticleSystemConfig::ParticleInitialSetup< ParticlesConfig >;

   /**
    * Particle and neighbor search model.
    */
   using ParticleSystem = ParticlesLinkedList< ParticlesConfig, SimulationControl::DeviceType >;

   /**
    * Define simulation SPH model and SPH formulation.
    *
    * - SPHModel: is the model of used SPH method (WCSPH_DBC, WCSPH_BI, RSPH, etc.)
    *   IMPORTANT: Constants and parameters of used model have to be defined in the SPHConfig.
    *
    * - SPHSimulation: defines the type of problem (simple fluid, problem with open or
    *   moving boundaries or multiphase flows). For the chosen type of simulation,
    *   appropriate SPH scheme is required!
    */
   using SPHModel = SPH::WCSPH_DBC< ParticleSystem, SPHConfig >;
   using SPHSimulation = SPH::SPHSimpleFluid< SPHModel >;

   /**
    * Define time step control.
    * There is const time step option and variable time step option.
    */
   //using TimeStepping = SPH::ConstantTimeStep< SPHConfig >;
   using TimeStepping = SPH::VariableTimeStep< SPHConfig >;

   /**
    * Define readers and writers to read and write initial geometry and results.
    */
   using Reader = Readers::VTKReader;
   using Writer = Writers::VTKWriter< ParticleSystem >;
   using SimulationReaderType = ReadParticles< ParticlesConfig, Reader >;

   /**
    * Create instance of SPHParams class, which is object holding all the
    * necessary SPH constants, informations about terms in particular scheme etc.
    */
   SPHParams sphParams;

   /**
    * Create instance of Simulation control class, which is object holding all the
    * information about end time, results saving times, paths to the input files
    * and paths to store results.
    */
   SimulationControl simulationControl;

   /**
    * Create instance of class with neccessary initial information to create particle system
    * and thus to initialize SPH simulation.
    */
   ParticlesParams particlesParams;

   /**
    * Create the main object - SPH simulation itself. The constructor requires
    * struct containing information to create and allocate particle system and neighbor search,
    * which includes number of particles for fluid and boundary, background grid size and its
    * origin and search radius.
    */
   SPHSimulation sph( particlesParams );
   std::cout << sph << std::endl;

   /**
    * Create instance of timeStepper, which is a class controling the time step,
    * duration of the simulation etc.
    *
    * Add output timer to control saving to files.
    */
   TimeStepping timeStepping( sphParams.dtInit, simulationControl.endTime );
   timeStepping.addOutputTimer( "save_results", simulationControl.outputTime );

   /**
    * Read the particle file.
    *
    * Read particle file with fluid and read/set initial particle variables.
    * Read particle file with boundary and read/set initial particle variables.
    */
   sph.fluid->template readParticlesAndVariables< SimulationReaderType >(
         simulationControl.inputParticleFile );
   sph.boundary->template readParticlesAndVariables< SimulationReaderType >(
         simulationControl.inputParticleFile_bound );

   /**
    * User defined measuretool sensors. Load and initialize configuration
    * for given type of measurement/sensors.
    *
    * - Load measuretool configurations:
    *   config for pressure measurement
    *    - from particles to grid interpolation.
    *
    *   config for water level measurement
    *    - sensors to measure water levels in given points
    *
    *   config for grid interpolation
    *    - sensors to measure pressure in given points
    */
   using GridInterpolation = SPH::InterpolateToGrid< SPHConfig, SPHSimulation >;
   using MeasuretoolInitGridInterpolation = SPH::MeasuretoolConfiguration::GridInterpolationConfig< SPHConfig >;
   MeasuretoolInitGridInterpolation interpolateGridParams;
   GridInterpolation interpolator( interpolateGridParams );

   using SensorInterpolation = SPH::SensorInterpolation< SPHConfig, SPHSimulation >;
   using MeasuretoolInitParametersPressure = SPH::MeasuretoolConfiguration::MeasuretoolConfigForPressure< SPHConfig >;
   MeasuretoolInitParametersPressure measuretoolPressure;
   timeStepping.addOutputTimer( "sensor_pressure", measuretoolPressure.outputTime );
   SensorInterpolation sensorInterpolation( TNL::ceil( simulationControl.endTime / measuretoolPressure.outputTime ),
         measuretoolPressure.points );

   using MeasuretoolInitParametersWaterLevel = SPH::MeasuretoolConfiguration::MeasuretoolConfigForWaterLevel< SPHConfig >;
   using SensorWaterLevel = SPH::SensorWaterLevel< SPHConfig, SPHSimulation >;
   MeasuretoolInitParametersWaterLevel measuretoolWaterLevel;
   timeStepping.addOutputTimer( "sensor_waterLevel", measuretoolWaterLevel.outputTime );
   SensorWaterLevel sensorWaterLevel( TNL::ceil( simulationControl.endTime / measuretoolWaterLevel.outputTime ), measuretoolWaterLevel.points,
         sphParams.h, measuretoolWaterLevel.direction, measuretoolWaterLevel.startMeasureAtLevel, measuretoolWaterLevel.stopMeasureAtLevel );

   /**
    * Define timers to measure computation time.
    */
   TNL::Timer timer_search, timer_interact, timer_computeTimeStep, timer_integrate, timer_pressure;
   TNL::Timer timer_search_reset, timer_search_cellIndices, timer_search_sort, timer_search_toCells;

   while( timeStepping.runTheSimulation() )
   {
      std::cout << "Time: " << timeStepping.getTime() << std::endl;

      /**
       * Find neighbors within the SPH simulation.
       */
      timer_search.start();
      sph.PerformNeighborSearch(
            timeStepping.getStep(), timer_search_reset, timer_search_cellIndices, timer_search_sort, timer_search_toCells );
      timer_search.stop();
      std::cout << "Search... done. " << std::endl;

      /**
       * Perform interaction with given model.
       */
      timer_interact.start();
      sph.template interact< SPH::WendlandKernel2D, SPHParams::DiffusiveTerm, SPHParams::ViscousTerm, SPHParams::EOS >( sphParams );
      timer_interact.stop();
      std::cout << "Interact... done. " << std::endl;

      /**
       * Compute new time step.
       */
      timer_computeTimeStep.start();
      timeStepping.computeTimeStep( sph.fluid, sphParams );
      timer_computeTimeStep.stop();
      std::cout << "New time step is: " << timeStepping.getTimeStep() << " s." << std::endl;
      std::cout << "Time step... done." << std::endl;

      /**
       * Perform time integration, i.e. update particle positions.
       */
      timer_integrate.start();
      sph.integrator->integratStepVerlet( sph.fluid, sph.boundary, timeStepping );
      timer_integrate.stop();
      std::cout << "Integrate... done. " << std::endl;

      /**
       * Output particle data
       */
      if( timeStepping.checkOutputTimer( "save_results" ) )
      {
         /**
          * Compute pressure from density.
          * This is not necessary since we do this localy, if pressure is needed.
          * Its useful for output anyway
          */
         timer_pressure.start();
         sph.model->template ComputePressureFromDensity< SPHParams::EOS >( sph.fluid, sphParams );
         timer_pressure.stop();
         std::cout << "Compute pressure... done. " << std::endl;

         timer_pressure.start();
         sph.model->template ComputePressureFromDensity< SPHParams::EOS >( sph.boundary, sphParams ); //TODO: FIX.
         timer_pressure.stop();
         std::cout << "Compute pressure... done. " << std::endl;

         sph.template save< Writer >( simulationControl.outputFileName, timeStepping.getStep() );

         /**
          * Interpolate on the grid.
          */
         std::string outputFileNameInterpolation = simulationControl.outputFileName + std::to_string( timeStepping.getStep() ) + "_interpolation.vtk";
         interpolator.template interpolate< SPH::WendlandKernel2D >( sph.fluid, sph.boundary, sphParams );
         interpolator.save( outputFileNameInterpolation );

      }

      if( timeStepping.checkOutputTimer( "sensor_pressure" ) )
      {
         sensorInterpolation.template interpolate< SPH::WendlandKernel2D, SPHParams::EOS >(
               sph.fluid, sph.boundary, sphParams, measuretoolPressure.includeBoundary );
      }

      if( timeStepping.checkOutputTimer( "sensor_waterLevel" ) )
      {
         sensorWaterLevel.template interpolate< SPH::WendlandKernel2D, SPHParams::EOS >(
               sph.fluid, sph.boundary, sphParams );
      }

      timeStepping.updateTimeStep();
   }

   std::string outputFileNameInterpolation = simulationControl.outputFileName + "_sensors.dat";
   sensorInterpolation.save( outputFileNameInterpolation );

   std::string outputFileNameWaterLevel = simulationControl.outputFileName + "_sensorsWaterLevel.dat";
   sensorWaterLevel.save( outputFileNameWaterLevel );

   /**
    * Output simulation stats.
    */
   float totalTime = ( timer_search.getRealTime() + \
   + timer_interact.getRealTime() + timer_integrate.getRealTime() + timer_pressure.getRealTime() );

   int steps = timeStepping.getStep();
   float totalTimePerStep = totalTime / steps;

   std::cout << std::endl << "COMPUTATION TIME:" << std::endl;
   std::cout << "Search........................................ " << timer_search.getRealTime() << " sec." << std::endl;
   std::cout << "Search (average time per step)................ " << timer_search.getRealTime() / steps << " sec." << std::endl;
   std::cout << "Search (percentage)........................... " << timer_search.getRealTime() / totalTime * 100 << " %." << std::endl;
   std::cout << " - Reset ..................................... " << timer_search_reset.getRealTime() << " sec." << std::endl;
   std::cout << " - Reset (average time per step).............. " << timer_search_reset.getRealTime() / steps << " sec." << std::endl;
   std::cout << " - Reset (percentage)......................... " << timer_search_reset.getRealTime() / totalTime * 100 << " %." << std::endl;
   std::cout << " - Index by cell ............................. " << timer_search_cellIndices.getRealTime() << " sec." << std::endl;
   std::cout << " - Index by cell (average time per step)...... " << timer_search_cellIndices.getRealTime() / steps << " sec." << std::endl;
   std::cout << " - Index by cell (percentage)................. " << timer_search_cellIndices.getRealTime() / totalTime * 100 << " %." << std::endl;
   std::cout << " - Sort ...................................... " << timer_search_sort.getRealTime() << " sec." << std::endl;
   std::cout << " - Sort (average time per step)............... " << timer_search_sort.getRealTime() / steps << " sec." << std::endl;
   std::cout << " - Sort (percentage).......................... " << timer_search_sort.getRealTime() / totalTime * 100 << " %." << std::endl;
   std::cout << " - Particle to cell .......................... " << timer_search_toCells.getRealTime() << " sec." << std::endl;
   std::cout << " - Particle to cell (average time per step)... " << timer_search_toCells.getRealTime() / steps << " sec." << std::endl;
   std::cout << " - Particle to cell (percentage).............. " << timer_search_toCells.getRealTime() / totalTime * 100 << " %." << std::endl;
   std::cout << "Interaction................................... " << timer_interact.getRealTime() << " sec." << std::endl;
   std::cout << "Interaction (average time per step)........... " << timer_interact.getRealTime() / steps << " sec." << std::endl;
   std::cout << "Interaction (percentage)...................... " << timer_interact.getRealTime() / totalTime * 100 << " %." << std::endl;
   std::cout << "Integrate..................................... " << timer_integrate.getRealTime() << " sec." << std::endl;
   std::cout << "Integrate (average time per step)............. " << timer_integrate.getRealTime() / steps << " sec." << std::endl;
   std::cout << "Integrate (percentage)........................ " << timer_integrate.getRealTime() / totalTime * 100 << " %." << std::endl;
   std::cout << "Pressure update............................... " << timer_pressure.getRealTime() << " sec." << std::endl;
   std::cout << "Pressure update (average time per step)....... " << timer_pressure.getRealTime() / steps << " sec." << std::endl;
   std::cout << "Pressure update (percentage).................. " << timer_pressure.getRealTime() / totalTime * 100 << " %." << std::endl;
   std::cout << "Total......................................... " << ( timer_search.getRealTime() + \
   + timer_interact.getRealTime() + timer_integrate.getRealTime() + timer_pressure.getRealTime() ) << " sec." << std::endl;
   std::cout << "Total (average time per step)................. " << ( timer_search.getRealTime() + \
   + timer_interact.getRealTime() + timer_integrate.getRealTime() + timer_pressure.getRealTime() ) / steps << " sec." << std::endl;

   //JsonMap
   std::map< std::string, std::string > timeResults;

   timeResults.insert({ "search",                              std::to_string( timer_search.getRealTime()                                ) } );
   timeResults.insert({ "search-average",                      std::to_string( timer_search.getRealTime() / steps                        ) } );
   timeResults.insert({ "search-percentage",                   std::to_string( timer_search.getRealTime() / totalTime * 100              ) } );
   timeResults.insert({ "search-reset",                        std::to_string( timer_search_reset.getRealTime()                          ) } );
   timeResults.insert({ "search-reset-average",                std::to_string( timer_search_reset.getRealTime() / steps                  ) } );
   timeResults.insert({ "search-reset-percentage",             std::to_string( timer_search_reset.getRealTime() / totalTime * 100        ) } );
   timeResults.insert({ "search-index-by-cell ",               std::to_string( timer_search_cellIndices.getRealTime()                    ) } );
   timeResults.insert({ "search-index-by-cell-average",        std::to_string( timer_search_cellIndices.getRealTime() / steps            ) } );
   timeResults.insert({ "search-index-by-cell-percentage",     std::to_string( timer_search_cellIndices.getRealTime() / totalTime * 100  ) } );
   timeResults.insert({ "search-sort",                         std::to_string( timer_search_sort.getRealTime()                           ) } );
   timeResults.insert({ "search-sort-average",                 std::to_string( timer_search_sort.getRealTime() / steps                   ) } );
   timeResults.insert({ "search-sort-percentage",              std::to_string( timer_search_sort.getRealTime() / totalTime * 100         ) } );
   timeResults.insert({ "search-particles-to-cell ",           std::to_string( timer_search_toCells.getRealTime()                        ) } );
   timeResults.insert({ "search-particles-to-cell-average",    std::to_string( timer_search_toCells.getRealTime() / steps                ) } );
   timeResults.insert({ "search-particles-to-cell-percentage", std::to_string( timer_search_toCells.getRealTime() / totalTime * 100      ) } );
   timeResults.insert({ "interaction",                         std::to_string( timer_interact.getRealTime()                              ) } );
   timeResults.insert({ "interaction-average",                 std::to_string( timer_interact.getRealTime() / steps                      ) } );
   timeResults.insert({ "interaction-percentage",              std::to_string( timer_interact.getRealTime() / totalTime * 100            ) } );
   timeResults.insert({ "integrate",                           std::to_string( timer_integrate.getRealTime()                             ) } );
   timeResults.insert({ "integrate-average",                   std::to_string( timer_integrate.getRealTime() / steps                     ) } );
   timeResults.insert({ "integrate-percentage",                std::to_string( timer_integrate.getRealTime() / totalTime * 100           ) } );
   timeResults.insert({ "pressure-update",                     std::to_string( timer_pressure.getRealTime()                              ) } );
   timeResults.insert({ "pressure-update-average",             std::to_string( timer_pressure.getRealTime() / steps                      ) } );
   timeResults.insert({ "pressure-update-percentage",          std::to_string( timer_pressure.getRealTime() / totalTime * 100            ) } );
   timeResults.insert({ "total ",                              std::to_string( totalTime                                                 ) } );
   timeResults.insert({ "total-average",                       std::to_string( totalTime / steps                                         ) } );

   TNL::Benchmarks::writeMapAsJson( timeResults, "results/time_measurements", ".json" );
   std::cout << "\nDone ... " << std::endl;

}

