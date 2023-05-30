#include <iostream>
#include <fstream> //temp, to write output

#include <TNL/Devices/Cuda.h>
#include <string>
#include <sys/types.h>

/**
 *  Benchamrk stuff.
 */
#include <TNL/Benchmarks/Benchmarks.h>

/**
 * Particle system.
 */
#include "../../../Particles/Particles.h"
#include "../../../Particles/neighbourSearch.h"

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
#include "ParticlesConfig.h"
#include "SPHCaseConfig.h"
#include "MeasuretoolConfig.h"
#include "SimulationControlConfig.h"

/**
 * SPH general toolds.
 */
#include "../../../SPH/SPH.h"

/**
 * SPH model.
 */
#include "../../../SPH/Models/WCSPH_DBC/Variables.h"
#include "../../../SPH/Models/WCSPH_DBC/Interactions.h"
#include "../../../SPH/Models/EquationOfState.h"

#include "../../../SPH/Models/EquationOfState.h"
#include "../../../SPH/Models/DiffusiveTerms.h"
#include "../../../SPH/Kernels.h"

/**
 * Time step control.
 */
#include "../../../SPH/TimeStep.h"

/**
 * Measuretool draft.
 */
#include "../../../SPH/Models/WCSPH_DBC/measuretool/Measuretool.h"

using namespace TNL;

int main( int argc, char* argv[] )
{
   /**
    * Number of particles
    */
   using Device = Devices::Cuda;

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
    *
    * - Measuretool:
    *   config for pressure measurement
    *   config for water level measurement
    *   config for grid interpolation
    */
   using ParticlesConfig = ParticleSystemConfig::ParticleSystemConfig< Device >;
   using SPHSimulationConfig = TNL::ParticleSystem::SPH::SPHSimpleFluidConfig< ParticlesConfig >;

   using ParticlesInitParameters = ParticleSystem::ParticleSystemConfig::ParticleInitialSetup;

   using SPHConfig = SPH::SPHCaseConfig< Device >;

   using SimulationControl = TNL::ParticleSystem::SPH::SimulationConstrolConfiguration::SPHSimulationControl;

   using MeasuretoolInitGridInterpolation = TNL::ParticleSystem::SPH::MeasuretoolConfiguration::GridInterpolationConfig< SPHConfig >;
   using MeasuretoolInitParametersPressure = TNL::ParticleSystem::SPH::MeasuretoolConfiguration::MeasuretoolConfigForPressure< SPHConfig >;
   using MeasuretoolInitParametersWaterLevel = TNL::ParticleSystem::SPH::MeasuretoolConfiguration::MeasuretoolConfigForWaterLevel< SPHConfig >;

   /**
    * Particle and neighbor search model.
    */
   using ParticleSystem = typename ParticleSystem::Particles< ParticlesConfig, Device >;
   using NeighborSearch = typename TNL::ParticleSystem::NeighborSearch< ParticlesConfig, ParticleSystem >;

   /**
    * Define simulation SPH model and SPH formulation.
    *
    * - SPHModel: is the model of used SPH method (WCSPH_DBC, WCSPH_BI, RSPH, etc.)
    *   IMPORATANT: Constants and parameters of used model have to be defined in the SPHConfig.
    *
    * - SPHSimulation: defines the type of problem (simple fluid, problem with open or
    *   moving boundaries or multiphase flows). For the chosen type of simulation,
    *   appropriate SPH scheme is required!
    */
   using SPHModel = TNL::ParticleSystem::SPH::WCSPH_DBC< NeighborSearch, SPHConfig >;
   using SPHSimulation = TNL::ParticleSystem::SPH::SPHSimpleFluid< SPHModel, ParticleSystem, NeighborSearch >;

   /**
    * Define particular schemes of the used SPH model.
    * Here in the case of WCSPH, we work with DiffusiveTerm, ViscousTerm and EquationOfState.
    *
    * IMPORATANT: Constants and parameters of used schemes have to be defined in the SPHConfig.
    */
   using DiffusiveTerm = TNL::ParticleSystem::SPH::MolteniDiffusiveTerm< SPHConfig >;
   using ViscousTerm = TNL::ParticleSystem::SPH::ArtificialViscosity< SPHConfig >;
   using EOS = TNL::ParticleSystem::SPH::TaitWeaklyCompressibleEOS< SPHConfig >;

   /**
    * Define tools of measuretool to evaluate certain values from the simulation.
    *
    * - Grid interpolation: In given plane, interpolate data back to continuum.
    *
    * - Sensors: Get time series of variables in given point:
    *   pressure sensors
    *   water level sensors
    */
   using GridInterpolation = TNL::ParticleSystem::SPH::GridInterpolation< SPHConfig, SPHSimulation >;

   using SensorInterpolation = TNL::ParticleSystem::SPH::SensorInterpolation< SPHConfig, SPHSimulation >;
   using SensorWaterLevel = TNL::ParticleSystem::SPH::SensorWaterLevel< SPHConfig, SPHSimulation >;

   /**
    * Define readers and writers to read and write initial geometry and results.
    */
   using Reader = TNL::ParticleSystem::Readers::VTKReader;
   using Writer = TNL::ParticleSystem::Writers::VTKWriter< ParticleSystem >;

   /**
    * Define time step control.
    * There is const time step option and variable time step option.
    */
   using TimeStepping = TNL::ParticleSystem::SPH::ConstantTimeStep< SPHConfig >;

   /**
    * Load simulation parameters.
    */
   SPHSimulationConfig mySPHSimulationConfig;
   mySPHSimulationConfig.template loadParameters< ParticlesInitParameters >();

   /**
    * Load simulation control (file names, time steps,...)
    */
   SimulationControl mySimulationControl;

   /**
    * Create the simulation.
    */
   SPHSimulation mySPHSimulation( mySPHSimulationConfig );
   std::cout << mySPHSimulation << std::endl;

   /**
    * TEMP: Determine number of interation for constant timestep.
    * Perform simulation main loop.
    */
   TimeStepping myTimeStepping( SPHConfig::dtInit, mySimulationControl.endTime );

   /**
    * Read the particle file.
    *
    * Read particle file with fluid and read/set initial particle variables.
    * Read particle file with boundary and read/set initial particle variables.
    */
   TNL::ParticleSystem::ReadParticles< ParticlesConfig, Reader > myFluidReader(
         mySimulationControl.inputParticleFile, mySPHSimulationConfig.sizeFluid, mySPHSimulationConfig.sizeAllocatedFluid );
   myFluidReader.template readParticles< ParticleSystem::PointArrayType >( mySPHSimulation.fluid->particles->getPoints() ) ;

   myFluidReader.template readParticleVariable< SPHModel::ScalarArrayType, float >(
         mySPHSimulation.fluid->getFluidVariables()->rho, "Density" );
   myFluidReader.template readParticleVariable< SPHModel::ScalarArrayType, float >(
         mySPHSimulation.fluid->getFluidVariables()->p, "Pressure" );
   myFluidReader.template readParticleVariable< SPHModel::VectorArrayType, float >(
         mySPHSimulation.fluid->getFluidVariables()->v, "Velocity" );

   TNL::ParticleSystem::ReadParticles< ParticlesConfig, Reader > myBoundaryReader(
         mySimulationControl.inputParticleFile_bound, mySPHSimulationConfig.sizeBoundary, mySPHSimulationConfig.sizeAllocatedBoundary );
   myBoundaryReader.template readParticles< ParticleSystem::PointArrayType >( mySPHSimulation.boundary->particles->getPoints() ) ;

   myBoundaryReader.template readParticleVariable< SPHModel::ScalarArrayType, float >(
         mySPHSimulation.boundary->getBoundaryVariables()->rho, "Density" );
   myBoundaryReader.template readParticleVariable< SPHModel::ScalarArrayType, float >(
         mySPHSimulation.boundary->getBoundaryVariables()->p, "Pressure" );
   myBoundaryReader.template readParticleVariable< SPHModel::VectorArrayType, float >(
         mySPHSimulation.boundary->getBoundaryVariables()->v, "Velocity" );

   /**
    * Define measuretool sensors.
    * - from particles to grid interpolation.
    * - sensors to measure water levels in given points
    * - sensors to measure pressure in given points
    */
   float saveResultsTimer = 0.f;
   MeasuretoolInitGridInterpolation measuretoolInterpolation;
   GridInterpolation myInterpolation( measuretoolInterpolation.gridOrigin, measuretoolInterpolation.gridSize, measuretoolInterpolation.gridStep );

   MeasuretoolInitParametersPressure measuretoolPressure;
   float measuretoolPressureTimer = 0.f;
   SensorInterpolation mySensorInterpolation( TNL::ceil( mySimulationControl.endTime / measuretoolPressure.outputTime ), measuretoolPressure.points );

   MeasuretoolInitParametersWaterLevel measuretoolWaterLevel;
   float measuretoolWaterLevelTimer = 0.f;
   SensorWaterLevel mySensorWaterLevel( TNL::ceil( mySimulationControl.endTime / measuretoolWaterLevel.outputTime ), measuretoolWaterLevel.points, SPHConfig::h, measuretoolWaterLevel.direction, measuretoolWaterLevel.startMeasureAtLevel, measuretoolWaterLevel.stopMeasureAtLevel );

   /**
    * Define timers to measure computation time.
    */
   TNL::Timer timer_search, timer_interact, timer_integrate, timer_pressure;
   TNL::Timer timer_search_reset, timer_search_cellIndices, timer_search_sort, timer_search_toCells;


   //for( unsigned int iteration = 0; iteration < steps; iteration ++ )
   while( myTimeStepping.runTheSimulation() )
   {
      std::cout << "Time: " << myTimeStepping.getTime() << std::endl;

      /**
       * Find neighbors within the SPH simulation.
       */
      timer_search.start();
      mySPHSimulation.PerformNeighborSearch(
            myTimeStepping.getStep(), timer_search_reset, timer_search_cellIndices, timer_search_sort, timer_search_toCells );
      timer_search.stop();
      std::cout << "Search... done. " << std::endl;

      /**
       * Perform interaction with given model.
       */
      timer_interact.start();
      mySPHSimulation.template Interact< SPH::WendlandKernel2D, DiffusiveTerm, ViscousTerm, EOS >();
      timer_interact.stop();
      std::cout << "Interact... done. " << std::endl;

      /**
       * Perform time integration, i.e. update particle positions.
       */
      timer_integrate.start();
      if( myTimeStepping.getStep() % 20 == 0 ) {
         mySPHSimulation.integrator->IntegrateEuler< typename SPHSimulation::FluidPointer >( SPHConfig::dtInit, mySPHSimulation.fluid );
         mySPHSimulation.integrator->IntegrateEulerBoundary< typename SPHSimulation::BoundaryPointer >( SPHConfig::dtInit, mySPHSimulation.boundary );
      }
      else {
         mySPHSimulation.integrator->IntegrateVerlet< typename SPHSimulation::FluidPointer >( SPHConfig::dtInit, mySPHSimulation.fluid );
         mySPHSimulation.integrator->IntegrateVerletBoundary< typename SPHSimulation::BoundaryPointer >( SPHConfig::dtInit, mySPHSimulation.boundary );
      }
      timer_integrate.stop();

      /**
       * Output particle data
       */
      if( myTimeStepping.getTime() > saveResultsTimer )
      {
         saveResultsTimer += mySimulationControl.outputTime;

         /**
          * Compute pressure from density.
          * This is not necessary since we do this localy, if pressure is needed.
          * Its useful for output anywal.
          */
         timer_pressure.start();
         mySPHSimulation.model->template ComputePressureFromDensity< EOS >( mySPHSimulation.fluid->variables, mySPHSimulation.fluid->particles->getNumberOfParticles() ); //TODO: FIX.
         timer_pressure.stop();
         std::cout << "Compute pressure... done. " << std::endl;

         std::string outputFileNameFluid = mySimulationControl.outputFileName + std::to_string( myTimeStepping.getStep() ) + "_fluid.vtk";
         std::ofstream outputFileFluid ( outputFileNameFluid, std::ofstream::out );
         Writer myWriter( outputFileFluid );
         myWriter.writeParticles( *mySPHSimulation.fluid->particles );
         myWriter.template writePointData< SPHModel::ScalarArrayType >(
               mySPHSimulation.fluid->getFluidVariables()->p, "Pressure", mySPHSimulation.fluid->particles->getNumberOfParticles(), 1 );
         myWriter.template writeVector< SPHModel::VectorArrayType, SPHConfig::RealType >(
               mySPHSimulation.fluid->getFluidVariables()->v, "Velocity", 3, mySPHSimulation.fluid->particles->getNumberOfParticles() );

         timer_pressure.start();
         mySPHSimulation.model->template ComputePressureFromDensity< EOS >( mySPHSimulation.boundary->variables, mySPHSimulation.boundary->particles->getNumberOfParticles() ); //TODO: FIX.
         timer_pressure.stop();
         std::cout << "Compute pressure... done. " << std::endl;

         std::string outputFileNameBound = mySimulationControl.outputFileName + std::to_string( myTimeStepping.getStep() ) + "_boundary.vtk";
         std::ofstream outputFileBound ( outputFileNameBound, std::ofstream::out );
         Writer myWriterBoundary( outputFileBound );
         myWriterBoundary.writeParticles( *mySPHSimulation.boundary->particles );
         myWriterBoundary.template writePointData< SPHModel::ScalarArrayType >(
               mySPHSimulation.boundary->getBoundaryVariables()->p, "Pressure", mySPHSimulation.boundary->particles->getNumberOfParticles(), 1 );
         myWriterBoundary.template writeVector< SPHModel::VectorArrayType, SPHConfig::RealType >(
               mySPHSimulation.boundary->getBoundaryVariables()->v, "Velocity", 3, mySPHSimulation.boundary->particles->getNumberOfParticles() );

         /**
          * Interpolate on the grid.
          */
         std::string outputFileNameInterpolation = mySimulationControl.outputFileName + std::to_string( myTimeStepping.getStep() ) + "_interpolation.vtk";
         myInterpolation.template InterpolateGrid< SPH::WendlandKernel2D >( mySPHSimulation.fluid, mySPHSimulation.boundary );
         myInterpolation.saveInterpolation( outputFileNameInterpolation );

      }

      if( myTimeStepping.getTime() > measuretoolPressureTimer )
      {
         measuretoolPressureTimer += measuretoolPressure.outputTime;
         mySensorInterpolation.template interpolateSensors< SPH::WendlandKernel2D,
                                                            EOS >( mySPHSimulation.fluid, mySPHSimulation.boundary );
      }

      if( myTimeStepping.getTime() > measuretoolWaterLevelTimer )
      {
         measuretoolWaterLevelTimer += measuretoolWaterLevel.outputTime;
         mySensorWaterLevel.template interpolateSensors< SPH::WendlandKernel2D,
                                                         EOS >( mySPHSimulation.fluid, mySPHSimulation.boundary );
      }

      myTimeStepping.updateTimeStep();
   }

   std::string outputFileNameInterpolation = mySimulationControl.outputFileName + "_sensors.dat";
   mySensorInterpolation.saveSensors( outputFileNameInterpolation );

   std::string outputFileNameWaterLevel = mySimulationControl.outputFileName + "_sensorsWaterLevel.dat";
   mySensorWaterLevel.saveSensors( outputFileNameWaterLevel );

   /**
    * Output simulation stats.
    */
   float totalTime = ( timer_search.getRealTime() + \
   + timer_interact.getRealTime() + timer_integrate.getRealTime() + timer_pressure.getRealTime() );

   int steps = myTimeStepping.getStep();
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

   TNL::Benchmarks::writeMapAsJson( timeResults, "time_measurements", ".json" );
   std::cout << "\nDone ... " << std::endl;

}

