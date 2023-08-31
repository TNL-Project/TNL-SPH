#include <iostream>
#include <fstream> //temp, to write output

#include <TNL/Devices/Cuda.h>
#include <ostream>
#include <string>
#include <sys/types.h>

#include <TNL/MPI/ScopedInitializer.h>

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
#include "sources/SimulationControlConfig.h"

/**
 * SPH general toolds.
 */
#include "../../../SPH/DistributedSPH.h"

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

   TNL::MPI::ScopedInitializer mpi(argc, argv);
#ifdef HAVE_MPI
      std::cout << "Running with MPI." << std::endl;
#endif

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
   using ParticlesParams = ParticleSystemConfig::DistributedInitialParticleSetup< ParticlesConfig >;

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
   using SPHSimulation = SPH::SPHSimpleFluid< SPHModel >;

   /**
    * Define readers and writers to read and write initial geometry and results.
    */
   using Reader = Readers::VTKReader;
   using Writer = Writers::VTKWriter< ParticleSystem >;
   using SimulationReaderType = ReadParticles< ParticlesConfig, Reader >;

   /**
    * Distributed simulation
    */
   using DistributedSPHSimulation = SPH::DistributedSPHSimpleFluid< SPHSimulation >;

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
   ParticlesParams allParticleParams;

   /**
    * Create the simulation.
    */
   SPHSimulation localSPHSimulation( allParticleParams.particlesParams[ TNL::MPI::GetRank() ] );
   std::cout << localSPHSimulation << std::endl;

   /**
    * TEMP: Determine number of interation for constant timestep.
    * Perform simulation main loop.
    */
   SPHParams::TimeStepping timeStepping( sphParams.dtInit, simulationControl.endTime );
   timeStepping.addOutputTimer( "save_results", simulationControl.outputTime );

   /**
    * Read the particle file.
    *
    * Read particle file with fluid and read/set initial particle variables.
    * Read particle file with boundary and read/set initial particle variables.
    */
   localSPHSimulation.fluid->template readParticlesAndVariables< SimulationReaderType >(
         simulationControl.inputParticleFile[ TNL::MPI::GetRank() ] );
   localSPHSimulation.boundary->template readParticlesAndVariables< SimulationReaderType >(
         simulationControl.inputParticleFile_bound[ TNL::MPI::GetRank() ] );

   /**
    *
    */
   DistributedSPHSimulation distributedSph( std::move( localSPHSimulation ) );

   distributedSph.localSimulation.fluid->subdomainInfo.loadParameters( allParticleParams.subdomainParams[ TNL::MPI::GetRank() ] );
   distributedSph.localSimulation.boundary->subdomainInfo.loadParameters( allParticleParams.subdomainParams[ TNL::MPI::GetRank() ] );

   /**
    * Define timers to measure computation time.
    */
   TNL::Timer timer_search, timer_interact, timer_integrate, timer_pressure;
   TNL::Timer timer_search_reset, timer_search_cellIndices, timer_search_sort, timer_search_toCells;
   TNL::Timer timer_synchronize, timer_synchronize_updateInfo, timer_synchronize_transfer, timer_synchronize_arrange;

   distributedSph.localSimulation.fluid->centerObjectArraysInMemory();
   distributedSph.localSimulation.boundary->centerObjectArraysInMemory();


   TNL::MPI::Barrier( distributedSph.communicator );

   while( timeStepping.runTheSimulation() )
   {
      std::cout << "Time: " << timeStepping.getTime() << " step: " << timeStepping.getStep() << std::endl;

      TNL::MPI::Barrier( distributedSph.communicator );

      /**
       * Resize the domains based on the computation time
       * and numbers of particles.
       */
      if( ( timeStepping.getStep() > 0 ) && (  timeStepping.getStep() % 500 == 0 ) )
         distributedSph.performLoadBalancing();

      TNL::MPI::Barrier( distributedSph.communicator );

      /**
       * Find neighbors within the SPH simulation.
       */
      timer_search.start();
      distributedSph.localSimulation.PerformNeighborSearch(
            0, timer_search_reset, timer_search_cellIndices, timer_search_sort, timer_search_toCells );
      timer_search.stop();
      std::cout << "Search... done. " << std::endl;

      TNL::MPI::Barrier( distributedSph.communicator );

      /**
       * Update informations about subdomaints.
       */
      timer_synchronize_updateInfo.start();
      distributedSph.updateLocalSubdomain();
      timer_synchronize_updateInfo.stop();
      std::cout << "Update local simulation info... done. " << std::endl;

      TNL::MPI::Barrier( distributedSph.communicator );

      /**
       * Perform interaction with given model.
       */
      timer_interact.start();
      distributedSph.template interact< SPHParams::KernelFunction,
                                        SPHParams::DiffusiveTerm,
                                        SPHParams::ViscousTerm,
                                        SPHParams::EOS >( sphParams );
      timer_interact.stop();
      std::cout << "Interact... done. " << std::endl;

      /**
       * Perform time integration, i.e. update particle positions.
       */
      timer_integrate.start();
      distributedSph.localSimulation.integrator->integratStepVerlet(
            distributedSph.localSimulation.fluid,
            distributedSph.localSimulation.boundary,
            timeStepping );
      timer_integrate.stop();
      std::cout << "Integrate... done. " << std::endl;

      TNL::MPI::Barrier( distributedSph.communicator );

      /**
       * Transfer the data between domaints.
       */
      timer_synchronize_transfer.start();
      distributedSph.synchronize();
      timer_synchronize_transfer.stop();
      std::cout << "Synchronization... done. " << std::endl;

      TNL::MPI::Barrier( distributedSph.communicator );

      /**
       * Output particle data
       */
      if( timeStepping.checkOutputTimer( "save_results" ) )
      {
         /**
          * Compute pressure from density.
          * This is not necessary since we do this localy, if pressure is needed.
          * Its useful for output anywal.
          */
         timer_pressure.start();
         distributedSph.localSimulation.model->template computePressureFromDensity< SPHParams::EOS >(
               distributedSph.localSimulation.fluid, sphParams );
         timer_pressure.stop();
         std::cout << "Compute pressure... done. " << std::endl;

         distributedSph.template save< Writer >( simulationControl.outputFileName, timeStepping.getStep() );

      }

      timeStepping.updateTimeStep();
   }

   /**
    * Output simulation stats.
    */
   float totalTime = ( timer_search.getRealTime() + \
   + timer_interact.getRealTime() + timer_integrate.getRealTime() + timer_pressure.getRealTime() \
   + timer_synchronize_updateInfo.getRealTime() + timer_synchronize_transfer.getRealTime() + timer_synchronize.getRealTime() );

   float totalTimeSynchronize = timer_synchronize_updateInfo.getRealTime() \
   + timer_synchronize_transfer.getRealTime() + timer_synchronize.getRealTime();

   float totalTimeWithoutMPI = ( timer_search.getRealTime() + \
   + timer_interact.getRealTime() + timer_integrate.getRealTime() + timer_pressure.getRealTime() );

   int steps = timeStepping.getStep();
   float totalTimePerStep = totalTime / steps;


   TNL::MPI::Barrier( distributedSph.communicator );

   if( TNL::MPI::GetRank() == 0 )
   {
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
      std::cout << "MPI: Synchronize ............................. " << totalTimeSynchronize << " sec." << std::endl;
      std::cout << "MPI: Synchronize (average time per step)...... " << totalTimeSynchronize / steps << " sec." << std::endl;
      std::cout << "MPI: Synchronize (percentage)................. " << totalTimeSynchronize / totalTime * 100 << " %." << std::endl;
      std::cout << " - MPI: Update info .......................... " << timer_synchronize_updateInfo.getRealTime() << " sec." << std::endl;
      std::cout << " - MPI: Update info (average time per step)... " << timer_synchronize_updateInfo.getRealTime() / steps << " sec." << std::endl;
      std::cout << " - MPI: Update info (percentage).............. " << timer_synchronize_updateInfo.getRealTime() / totalTime * 100 << " %." << std::endl;
      std::cout << " - MPI: Transfer ............................. " << timer_synchronize_transfer.getRealTime() << " sec." << std::endl;
      std::cout << " - MPI: Transfer (average time per step)...... " << timer_synchronize_transfer.getRealTime() / steps << " sec." << std::endl;
      std::cout << " - MPI: Transfer (percentage)................. " << timer_synchronize_transfer.getRealTime() / totalTime * 100 << " %." << std::endl;
      std::cout << " - MPI: Arrange ............................. " << timer_synchronize_arrange.getRealTime() << " sec." << std::endl;
      std::cout << " - MPI: Arrange (average time per step)...... " << timer_synchronize_arrange.getRealTime() / steps << " sec." << std::endl;
      std::cout << " - MPI: Arrange (percentage)................. " << timer_synchronize_arrange.getRealTime() / totalTime * 100 << " %." << std::endl;
      std::cout << "Total without MPI............................. " << totalTimeWithoutMPI << " sec." << std::endl;
      std::cout << "Total without MPI (average time per step)..... " << totalTimeWithoutMPI / steps << " sec." << std::endl;
      std::cout << "Total......................................... " << totalTime << " sec." << std::endl;
      std::cout << "Total (average time per step)................. " << totalTime / steps << " sec." << std::endl;
   }

   TNL::MPI::Barrier( distributedSph.communicator );

   if( TNL::MPI::GetRank() == 1 )
   {
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
      std::cout << "MPI: Synchronize ............................. " << totalTimeSynchronize << " sec." << std::endl;
      std::cout << "MPI: Synchronize (average time per step)...... " << totalTimeSynchronize / steps << " sec." << std::endl;
      std::cout << "MPI: Synchronize (percentage)................. " << totalTimeSynchronize / totalTime * 100 << " %." << std::endl;
      std::cout << " - MPI: Update info .......................... " << timer_synchronize_updateInfo.getRealTime() << " sec." << std::endl;
      std::cout << " - MPI: Update info (average time per step)... " << timer_synchronize_updateInfo.getRealTime() / steps << " sec." << std::endl;
      std::cout << " - MPI: Update info (percentage).............. " << timer_synchronize_updateInfo.getRealTime() / totalTime * 100 << " %." << std::endl;
      std::cout << " - MPI: Transfer ............................. " << timer_synchronize_transfer.getRealTime() << " sec." << std::endl;
      std::cout << " - MPI: Transfer (average time per step)...... " << timer_synchronize_transfer.getRealTime() / steps << " sec." << std::endl;
      std::cout << " - MPI: Transfer (percentage)................. " << timer_synchronize_transfer.getRealTime() / totalTime * 100 << " %." << std::endl;
      std::cout << " - MPI: Arrange ............................. " << timer_synchronize_arrange.getRealTime() << " sec." << std::endl;
      std::cout << " - MPI: Arrange (average time per step)...... " << timer_synchronize_arrange.getRealTime() / steps << " sec." << std::endl;
      std::cout << " - MPI: Arrange (percentage)................. " << timer_synchronize_arrange.getRealTime() / totalTime * 100 << " %." << std::endl;
      std::cout << "Total without MPI............................. " << totalTimeWithoutMPI << " sec." << std::endl;
      std::cout << "Total without MPI (average time per step)..... " << totalTimeWithoutMPI / steps << " sec." << std::endl;
      std::cout << "Total......................................... " << totalTime << " sec." << std::endl;
      std::cout << "Total (average time per step)................. " << totalTime / steps << " sec." << std::endl;
   }

   TNL::MPI::Barrier( distributedSph.communicator );

   std::cout << "\nDone ... " << std::endl;
}

