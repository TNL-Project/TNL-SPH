#include <iostream>
#include <fstream> //temp, to write output

#include <TNL/Devices/Cuda.h>
#include <string>
#include <sys/types.h>

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
//#include "ParticlesConfigNew.h"
//#include "SPHCaseConfig.h"
//#include "MeasuretoolConfig.h"
//#include "SimulationControlConfig.h"
//#include "OpenBoundaryConfigNew.h"

#include "sources/ParticlesConfig.h"
#include "sources/SPHCaseConfig.h"
//#include "sources/MeasuretoolConfig.h"
#include "sources/SimulationControlConfig.h"
#include "sources/OpenBoundaryConfig.h"

/**
 * SPH general toolds.
 */
#include "../../../SPH/SPHOpen.h"

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

using namespace TNL::ParticleSystem;

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
    *   IMPORATANT: Constants and parameters of used model have to be defined in the SPHConfig.
    *
    * - SPHSimulation: defines the type of problem (simple fluid, problem with open or
    *   moving boundaries or multiphase flows). For the chosen type of simulation,
    *   appropriate SPH scheme is required!
    */
   using SPHModel = SPH::WCSPH_DBC< ParticleSystem, SPHParams >;
   using SPHSimulation = SPH::SPHOpenSystem< SPHModel >;

   /**
    * Define time step control.
    * There is const time step option and variable time step option.
    */
   using TimeStepping = TNL::ParticleSystem::SPH::ConstantTimeStep< SPHConfig >;

   /**
    * Particle reader and writer.
    */
   using Reader = TNL::ParticleSystem::Readers::VTKReader;
   using Writer = TNL::ParticleSystem::Writers::VTKWriter< ParticleSystem >;
   using SimulationReaderType = ReadParticles< ParticlesConfig, Reader >;

   /**
    * Create instance of SPHParams class, which is object holding all the
    * necessary SPH constants, informations about terms in particular scheme etc.
    */
   SPHParams sphParams;
   TNL::ParticleSystem::SPH::InletBuffer< SPHConfig > inletBufferParams;
   TNL::ParticleSystem::SPH::OutletBuffer< SPHConfig > outletBufferParams;

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

   //TODO: Resolve this
   sph.addOpenBoundaryPatch( particlesParams );

   /**
    * Read the particle file and setup the inlets and outlets.
    */
   sph.fluid->template readParticlesAndVariables< SimulationReaderType >(
         simulationControl.inputParticleFile );
   sph.fluid->variables->v = inletBufferParams.velocity;
   sph.fluid->variables->v_swap = inletBufferParams.velocity;

   sph.boundary->template readParticlesAndVariables< SimulationReaderType >(
         simulationControl.inputParticleFile_bound );

   sph.openBoundaryPatches[ 0 ]->template readParticlesAndVariables< SimulationReaderType >(
         simulationControl.inputParticleFile_inlet );
   sph.openBoundaryPatches[ 0 ]->readOpenBoundaryParameters( inletBufferParams );
   //FIXME:
   sph.openBoundaryPatches[ 0 ]->variables->v = inletBufferParams.velocity;
   sph.openBoundaryPatches[ 0 ]->variables->v_swap = inletBufferParams.velocity;

   sph.openBoundaryPatches[ 1 ]->template readParticlesAndVariables< SimulationReaderType >(
         simulationControl.inputParticleFile_outlet );
   sph.openBoundaryPatches[ 1 ]->readOpenBoundaryParameters( outletBufferParams );
   sph.openBoundaryPatches[ 1 ]->variables->v = inletBufferParams.velocity;
   sph.openBoundaryPatches[ 1 ]->variables->v_swap = inletBufferParams.velocity;


   std::cout << "Inlet parameters: " << std::endl;
   std::cout << "Orientation ................. " << sph.openBoundaryPatches[ 0 ]->parameters.orientation << std::endl;
   std::cout << "Velocity .................... " << sph.openBoundaryPatches[ 0 ]->parameters.velocity << std::endl;
   std::cout << "BufferWidth ................. " << sph.openBoundaryPatches[ 0 ]->parameters.bufferWidth << std::endl;

  std::cout << "Inlet2 parameters: " << std::endl;
  std::cout << "Orientation ................. " << sph.openBoundaryPatches[ 1 ]->parameters.orientation << std::endl;
  std::cout << "Velocity .................... " << sph.openBoundaryPatches[ 1 ]->parameters.velocity << std::endl;
  std::cout << "BufferWidth ................. " << sph.openBoundaryPatches[ 1 ]->parameters.bufferWidth << std::endl;

   /**
    * Define timers to measure computation time.
    */
   TNL::Timer timer_search, timer_interact, timer_integrate, timer_inlet, timer_outlet, timer_pressure;
   TNL::Timer timer_search_reset, timer_search_cellIndices, timer_search_sort, timer_search_toCells;

   while( timeStepping.runTheSimulation() )
   {
      std::cout << "Time: " << timeStepping.getTime() << " step: " <<  timeStepping.getStep() << std::endl;
      std::cout << "Number of fluid particles: " << sph.fluid->particles->getNumberOfParticles() << std::endl;

      /**
       * Find neighbors within the SPH simulation.
       */
      timer_search.start();
      sph.PerformNeighborSearch(
            timeStepping.getStep(), timer_search_reset, timer_search_cellIndices, timer_search_sort, timer_search_toCells );
      timer_search.stop();
      std::cout << "Search... done. " << std::endl;

      //TODO:
      sph.model->extrapolateOpenBoundaryData< SPHSimulation::FluidPointer,
                                              SPHSimulation::OpenBoundaryPointer,
                                              SPH::WendlandKernel2D,
                                              SPHParams::EOS >( sph.fluid, sph.openBoundaryPatches[ 1 ], sphParams );

      /**
       * Perform interaction with given model.
       */
      timer_interact.start();
      sph.template interact< SPH::WendlandKernel2D, SPHParams::DiffusiveTerm, SPHParams::ViscousTerm, SPHParams::EOS >( sphParams );
      timer_interact.stop();
      std::cout << "Interact... done. " << std::endl;

      //#include "outputForDebug.h"

      /**
       * Perform time integration, i.e. update particle positions.
       */
      timer_integrate.start();
      sph.integrator->integratStepVerlet( sph.fluid, sph.boundary, timeStepping );
      timer_integrate.stop();
      std::cout << "Integration... done. " << std::endl;

      timer_inlet.start();
      sph.integrator->updateBuffer< typename SPHSimulation::FluidPointer, typename SPHSimulation::OpenBoundaryPointer >( timeStepping.getTimeStep(), sph.fluid, sph.openBoundaryPatches[ 0 ] );
      timer_inlet.stop();
      timer_outlet.start();
      sph.integrator->updateOutletBuffer< typename SPHSimulation::FluidPointer, typename SPHSimulation::OpenBoundaryPointer >( timeStepping.getTimeStep(), sph.fluid, sph.openBoundaryPatches[ 1 ] );
      timer_outlet.stop();
      std::cout << "Open boundary... done. " << std::endl;


      /**
       * Output particle data
       */
      if( timeStepping.checkOutputTimer( "save_results" ) )
      {
         /**
          * Compute pressure from density.
          * This is not necessary since we do this localy, if pressure is needed.
          * It's useful for output anyway.
          */
         timer_pressure.start();
         sph.model->template computePressureFromDensity< SPHParams::EOS >( sph.fluid, sphParams );
         timer_pressure.stop();
         std::cout << "Compute pressure... done. " << std::endl;

         timer_pressure.start();
         sph.model->template computePressureFromDensity< SPHParams::EOS >( sph.boundary, sphParams );
         timer_pressure.stop();
         std::cout << "Compute pressure... done. " << std::endl;

         sph.template save< Writer >( simulationControl.outputFileName, timeStepping.getStep() );
      }

      timeStepping.updateTimeStep();
   }

   /**
    * Output simulation stats.
    */
   float openBoundaryTime = timer_inlet.getRealTime() + timer_outlet.getRealTime();

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
   std::cout << "Open boundary................................. " << openBoundaryTime << " sec." << std::endl;
   std::cout << "Open boundary (average time per step)......... " << openBoundaryTime / steps << " sec." << std::endl;
   std::cout << "Open boundary (percentage).................... " << openBoundaryTime / totalTime * 100 << " %." << std::endl;
   std::cout << " - Update inlet............................... " << timer_inlet.getRealTime() << " sec." << std::endl;
   std::cout << " - Update inlet (average time per step)....... " << timer_inlet.getRealTime() / steps << " sec." << std::endl;
   std::cout << " - Update inlet (percentage).................. " << timer_inlet.getRealTime() / totalTime * 100 << " %." << std::endl;
   std::cout << " - Update outlet.............................. " << timer_outlet.getRealTime() << " sec." << std::endl;
   std::cout << " - Update outlet (average time per step)...... " << timer_outlet.getRealTime() / steps << " sec." << std::endl;
   std::cout << " - Update outlet (percentage)................. " << timer_outlet.getRealTime() / totalTime * 100 << " %." << std::endl;
   std::cout << "Pressure update............................... " << timer_pressure.getRealTime() << " sec." << std::endl;
   std::cout << "Pressure update (average time per step)....... " << timer_pressure.getRealTime() / steps << " sec." << std::endl;
   std::cout << "Pressure update (percentage).................. " << timer_pressure.getRealTime() / totalTime * 100 << " %." << std::endl;
   std::cout << "Total......................................... " << ( timer_search.getRealTime() + \
   + timer_interact.getRealTime() + timer_integrate.getRealTime() + timer_pressure.getRealTime() ) << " sec." << std::endl;
   std::cout << "Total (average time per step)................. " << ( timer_search.getRealTime() + \
   + timer_interact.getRealTime() + timer_integrate.getRealTime() + timer_pressure.getRealTime() ) / steps << " sec." << std::endl;

   std::cout << "\nDone ... " << std::endl;
}

