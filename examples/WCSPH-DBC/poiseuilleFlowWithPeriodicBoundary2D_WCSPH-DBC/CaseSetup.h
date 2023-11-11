/**
 * Include configuration files containing data for case definition.
 * - "SimulationControlConfig.h" contains core informations to control the simulation
 * - "ParticleConfig.h" contains information about domain and sizes of problem
 * - "SPHCaseConfig.h" contains parameter of SPH method
 * - "OpenBoundaryConfig.h" contains parameter of SPH method
 */
#include "sources/SimulationControlConfig.h"
using SimulationControl = TNL::ParticleSystem::SPH::SimulationControlConfiguration::SPHSimulationControl;

#include "sources/ParticlesConfig.h"
using ParticlesParams = TNL::ParticleSystem::ParticleSystemConfig::ParticleInitialSetup< SimulationControl::DeviceType >;

#include "sources/SPHCaseConfig.h"
using SPHParams = TNL::ParticleSystem::SPH::SPHConfig::SPHParamsConfig< SimulationControl::DeviceType >;

#include "sources/PeriodicBoundaryConfig.h"
using LeftPeriodicityParams = TNL::ParticleSystem::SPH::PeriodicityLeftBuffer< SPHParams::SPHConfig >;
using RightPeriodicityParams = TNL::ParticleSystem::SPH::PeriodicityRightBuffer< SPHParams::SPHConfig >;

/**
 * Include type of particle system.
 */
#include <Particles/ParticlesLinkedListFloating.h>
using ParticlesSys = TNL::ParticleSystem::ParticlesLinkedList< ParticlesParams::ParticlesConfig, SimulationControl::DeviceType >;

/**
 * Include particular formulation of SPH method.
 */
#include <SPH/Models/WCSPH_DBC/Interactions.h>
using SPHModel = TNL::ParticleSystem::SPH::WCSPH_DBC< ParticlesSys, SPHParams >;

/**
 * Include type of SPH simulation.
 */
#include <SPH/SPHOpen.h>
using SPHSimulation = TNL::ParticleSystem::SPH::SPHOpenSystem< SPHModel >;

/**
 * Particle system reader.
 */
#include <Readers/VTKReader.h>
#include <Writers/VTKWriter.h>
#include <Readers/readSPHSimulation.h>
using Reader = TNL::ParticleSystem::Readers::VTKReader;
using Writer = TNL::ParticleSystem::Writers::VTKWriter< ParticlesSys >;
using SimulationReaderType = TNL::ParticleSystem::ReadParticles< ParticlesParams::ParticlesConfig, Reader >;

/**
 *  Used to write computation time to json format.
 */
#include <TNL/Benchmarks/Benchmarks.h>


int main( int argc, char* argv[] )
{

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
    * Create instance of classes with parameters for inlet and outlet buffers. This
    * contains configurations for open boundary patches.
    */
   LeftPeriodicityParams leftPeriodicityParams;
   RightPeriodicityParams rightPeriodicityParams;

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
   SPHParams::TimeStepping timeStepping( sphParams.dtInit, simulationControl.endTime );
   timeStepping.addOutputTimer( "save_results", simulationControl.outputTime );

   //TODO: Resolve this
   sph.addOpenBoundaryPatch( particlesParams );

   /**
    * Read the particle file and setup the inlets and outlets.
    */
   sph.fluid->template readParticlesAndVariables< SimulationReaderType >(
         simulationControl.inputParticleFile );

   sph.boundary->template readParticlesAndVariables< SimulationReaderType >(
         simulationControl.inputParticleFile_bound );

   sph.openBoundaryPatches[ 0 ]->readOpenBoundaryParameters( leftPeriodicityParams, particlesParams );
   sph.openBoundaryPatches[ 1 ]->readOpenBoundaryParameters( rightPeriodicityParams, particlesParams );

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

      sph.integrator->applyPeriodicBoundary( sph.fluid, sph.openBoundaryPatches[ 0 ], sph.openBoundaryPatches[ 1 ], leftPeriodicityParams.shift );
      std::cout << "Update periodic (copy data)... done. " << std::endl;

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
      sph.template interact< SPHParams::KernelFunction, SPHParams::DiffusiveTerm, SPHParams::ViscousTerm, SPHParams::EOS >( sphParams );
      timer_interact.stop();
      std::cout << "Interact... done. " << std::endl;

      /**0
       * Perform time integration, i.e. update particle positions.
       */
      timer_integrate.start();
      sph.integrator->integratStepVerlet( sph.fluid, sph.boundary, timeStepping );
      timer_integrate.stop();
      std::cout << "Integration... done. " << std::endl;

      sph.integrator->periodicityParticleTransfer( sph.fluid, sph.openBoundaryPatches[ 0 ], leftPeriodicityParams.shift );
      sph.integrator->periodicityParticleTransfer( sph.fluid, sph.openBoundaryPatches[ 1 ], rightPeriodicityParams.shift );
      std::cout << "Update periodic (movement)... done. " << std::endl;

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
    * Write out simulation computation time.
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

