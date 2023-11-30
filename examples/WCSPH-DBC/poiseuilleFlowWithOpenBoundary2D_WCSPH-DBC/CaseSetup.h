/**
 * Include configuration files containing data for case definition.
 * - "SimulationControlConfig.h" contains core informations to control the simulation
 * - "ParticleConfig.h" contains information about domain and sizes of problem
 * - "SPHCaseConfig.h" contains parameter of SPH method
 * - "OpenBoundaryConfig.h" contains parameter of SPH method
 */
#include "sources/SimulationControlConfig.h"
#include <string>
using SimulationControl = TNL::ParticleSystem::SPH::SimulationControlConfiguration::SPHSimulationControl;

#include "sources/ParticlesConfig.h"
using ParticlesParams = TNL::ParticleSystem::ParticleSystemConfig::ParticleInitialSetup< SimulationControl::DeviceType >;

#include "sources/SPHCaseConfig.h"
using SPHParams = TNL::ParticleSystem::SPH::SPHConfig::SPHParamsConfig< SimulationControl::DeviceType >;

#include "sources/OpenBoundaryConfig.h"
using InletParams = TNL::ParticleSystem::SPH::InletBuffer< SPHParams::SPHConfig >;
using OutletParams = TNL::ParticleSystem::SPH::OutletBuffer< SPHParams::SPHConfig >;

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
#include <SPH/TimeMeasurement.h>

#include <SPH/Models/WCSPH_DBC/OpenBoundaryConfig.h>

#include "sources/OpenBoundaryConfig.test.h"

/**
 * Test
 */
#include <SPH/Models/WCSPH_DBC/OpenBoundaryConfig.h>

int main( int argc, char* argv[] )
{

   /**
    * Create instance of SPHParams class, which is object holding all the
    * necessary SPH constants, informations about terms in particular scheme etc.
    */
   SPHParams sphParams;

   std::vector< SPHModel::OpenBoundaryConfig > openBoundaryConfigs( 2 );
   #include "sources/OpenBoundaryConfig.test2.h"

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
   InletParams inletParams;
   OutletParams outletParams;

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
   //sph.addOpenBoundaryPatch( particlesParams );
   sph.addOpenBoundaryPatch( particlesParams, openBoundaryConfigs );

   /**
    * Read the particle file and setup the inlets and outlets.
    */
   sph.fluid->template readParticlesAndVariables< SimulationReaderType >(
         simulationControl.inputParticleFile );

   sph.boundary->template readParticlesAndVariables< SimulationReaderType >(
         simulationControl.inputParticleFile_bound );

   sph.openBoundaryPatches[ 0 ]->template readParticlesAndVariables< SimulationReaderType >(
         simulationControl.inputParticleFile_inlet );
   sph.openBoundaryPatches[ 0 ]->readOpenBoundaryParameters( inletParams );

   sph.openBoundaryPatches[ 1 ]->template readParticlesAndVariables< SimulationReaderType >(
         simulationControl.inputParticleFile_outlet );
   sph.openBoundaryPatches[ 1 ]->readOpenBoundaryParameters( outletParams );

   //sph.template addOpenBoundaryPatch< TNL::ParticleSystem::SPH::WCSPH_BCTypes::Inlet >( "inlet", particlesParams, 0 );
   //sph.template addOpenBoundaryPatch< TNL::ParticleSystem::SPH::WCSPH_BCTypes::Outlet >( "outlet", particlesParams, 1 );


   //FIXME - problem with 2D vectors:
   sph.fluid->variables->v = inletParams.velocity;
   sph.fluid->variables->v_swap = inletParams.velocity;
   sph.openBoundaryPatches[ 0 ]->variables->v = inletParams.velocity;
   sph.openBoundaryPatches[ 0 ]->variables->v_swap = inletParams.velocity;
   sph.openBoundaryPatches[ 1 ]->variables->v = inletParams.velocity;
   sph.openBoundaryPatches[ 1 ]->variables->v_swap = inletParams.velocity;

   /**
    * Define timers to measure computation time.
    */
   TNL::Timer timer_search, timer_interact, timer_integrate, timer_inlet, timer_outlet, timer_pressure;
   TNL::Timer timer_search_reset, timer_search_cellIndices, timer_search_sort, timer_search_toCells;

   TNL::ParticleSystem::SPH::TimerMeasurement timeMeasurement;

   TNL::Logger log( 100, std::cout );

   //while( timeStepping.runTheSimulation() )
   while( timeStepping.getStep() < 3 )
   {
      log.writeHeader( "Time: " + std::to_string( timeStepping.getTime() ) + ", simulation step: " \
                       + std::to_string( timeStepping.getStep() ) );

      /**
       * Find neighbors within the SPH simulation.
       */
      timeMeasurement.start( "search" );
      sph.performNeighborSearch( timeStepping.getStep(), timeMeasurement, log );
      timeMeasurement.stop( "search" );
      log.writeParameter( "Search...", "Done." );

      sph.template extrapolateOpenBC< SPHParams::KernelFunction, SPHParams::EOS >( sphParams, openBoundaryConfigs );
      log.writeParameter( "Extrapolate open BC...", "Done." );

      /**
       * Perform interaction with given model.
       */
      timeMeasurement.start( "interact" );
      sph.template interact< SPHParams::KernelFunction, SPHParams::DiffusiveTerm, SPHParams::ViscousTerm, SPHParams::EOS >( sphParams );
      timeMeasurement.stop( "interact" );
      log.writeParameter( "Interact...", "Done." );

      /**
       * Perform time integration, i.e. update particle positions.
       */
      timeMeasurement.start( "integrate" );
      sph.integrator->integratStepVerlet( sph.fluid, sph.boundary, timeStepping );
      timeMeasurement.stop( "integrate" );
      log.writeParameter( "Integrate...", "Done." );

      timer_inlet.start();
      sph.integrator->applyOpenBoundary( timeStepping.getTimeStep(), sph.fluid, sph.openBoundaryPatches[ 0 ], openBoundaryConfigs[ 0 ] );
      timer_inlet.stop();
      timer_outlet.start();
      sph.integrator->applyOpenBoundary( timeStepping.getTimeStep(), sph.fluid, sph.openBoundaryPatches[ 1 ], openBoundaryConfigs[ 1 ] );
      timer_outlet.stop();
      log.writeParameter( "Update open BC...", "Done." );

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

         timer_pressure.start();
         sph.model->template computePressureFromDensity< SPHParams::EOS >( sph.boundary, sphParams );
         timer_pressure.stop();

         sph.template save< Writer >( simulationControl.outputFileName, timeStepping.getStep() );
      }

      timeStepping.updateTimeStep();
   }

   timeMeasurement.print( timeStepping.getStep() );
   std::cout << "\nDone ... " << std::endl;
}

