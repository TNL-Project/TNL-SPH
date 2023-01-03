#include <iostream>
#include <fstream> //temp, to write output

#include <TNL/Devices/Cuda.h>
#include <string>
#include <sys/types.h>

/**
 * Particle system.
 */
#include "../../Particles/Particles.h"
#include "../../Particles/neighbourSearch.h"

/**
 * Particle system reader.
 **/
#include "../../Readers/VTKReader.h"
#include "../../Writers/VTKWriter.h"

/**
 * Case configuration
 * One configuration for particle system, one for SPH.
 */
//#include "../SPHCaseSetup/damBreak2769particles/ParticlesConfig.h"
//#include "../SPHCaseSetup/damBreak2769particles/SPHCaseConfig.h"
//const std::string inputParticleFile = "../SPHCaseSetup/damBreak2769particles/dambreak.vtk";

//:#include "../SPHCaseSetup/damBreak49821particles/ParticlesConfig.h"
//:#include "../SPHCaseSetup/damBreak49821particles/SPHCaseConfig.h"
//:const std::string inputParticleFile = "../SPHCaseSetup/damBreak49821particles/dambreak.vtk";

#include "../SPHCaseSetup/damBreak49821particles/ParticlesConfig.h"
#include "../SPHCaseSetup/damBreak49821particles/SPHCaseConfig.h"
const std::string inputParticleFile = "../SPHCaseSetup/damBreak49821particles/dambreak_fluid.vtk";
const std::string inputParticleFile_bound = "../SPHCaseSetup/damBreak49821particles/dambreak_boundary.vtk";

//#include "../SPHCaseSetup/damBreakN49206particles/ParticlesConfig.h"
//#include "../SPHCaseSetup/damBreakN49206particles/SPHCaseConfig.h"
//const std::string inputParticleFile = "../SPHCaseSetup/damBreakN49206particles/dambreak.vtk";

//#include "../SPHCaseSetup/damBreak189636particles/ParticlesConfig.h"
//#include "../SPHCaseSetup/damBreak189636particles/SPHCaseConfig.h"
//const std::string inputParticleFile = "../SPHCaseSetup/damBreak189636particles/dambreak.vtk";

//#include "../SPHCaseSetup/damBreak189636particles/ParticlesConfig.h"
//#include "../SPHCaseSetup/damBreak189636particles/SPHCaseConfig.h"
//const std::string inputParticleFile = "../SPHCaseSetup/damBreak189636particles/dambreak_fluid.vtk";
//const std::string inputParticleFile_bound = "../SPHCaseSetup/damBreak189636particles/dambreak_boundary.vtk";

//#include "../SPHCaseSetup/damBreakN736806particles/ParticlesConfig.h"
//#include "../SPHCaseSetup/damBreakN736806particles/SPHCaseConfig.h"
//const std::string inputParticleFile = "../SPHCaseSetup/damBreakN736806particles/dambreak.vtk";

//#include "../SPHCaseSetup/damBreakN2913606particles/ParticlesConfig.h"
//#include "../SPHCaseSetup/damBreakN2913606particles/SPHCaseConfig.h"
//const std::string inputParticleFile = "../SPHCaseSetup/damBreakN2913606particles/dambreak.vtk";

//const float endTime = 0.0002;
const float endTime = 0.05;
const int outputStep = 2500;

std::string outputFileName = "results/particles";

/**
 * SPH general toolds.
 */
#include "../../SPH/SPH.h"

/**
 * SPH model.
 */
#include "../../SPH/Models/WCSPH_DBC/Variables.h"
#include "../../SPH/Models/WCSPH_DBC/Interactions.h"
#include "../../SPH/Models/EquationOfState.h"

#include "../../SPH/Models/EquationOfState.h"
#include "../../SPH/Models/DiffusiveTerms.h"
#include "../../SPH/Kernels.h"

using namespace TNL;

int main( int argc, char* argv[] )
{
   /**
    * Number of particles
    */
   using Device = Devices::Cuda;
   using ParticlesConfig = ParticleSystemConfig< Device >;
   using ParticlesConfig_bound = ParticleSystemConfig_boundary< Device >;
   using SPHConfig = SPH::SPHCaseConfig< Device >;

   /**
    * Particle and neighbor search model.
    */
   using ParticleSystem = typename ParticleSystem::Particles< ParticlesConfig, Device >;
   using NeighborSearch = typename TNL::ParticleSystem::NeighborSearch< ParticlesConfig, ParticleSystem >;

   /**
    * SPH model.
    */
   using SPHModel = typename TNL::ParticleSystem::SPH::WCSPH_DBC< ParticleSystem, SPHConfig >;
   using SPHSimulation = typename TNL::ParticleSystem::SPH::SPHSimpleFluid< SPHModel, ParticleSystem, NeighborSearch >;

   /**
    * SPH schemes.
    */
   using DiffusiveTerm = TNL::ParticleSystem::SPH::MolteniDiffusiveTerm< SPHConfig >;
   using ViscousTerm = TNL::ParticleSystem::SPH::ArtificialViscosity< SPHConfig >;
   using EOS = TNL::ParticleSystem::SPH::TaitWeaklyCompressibleEOS< SPHConfig >;

   /**
    * Particle reader and writer.
    */
   using Reader = TNL::ParticleSystem::Readers::VTKReader;
   using Writer = TNL::ParticleSystem::Writers::VTKWriter< ParticleSystem >;

   /**
    * Create the simulation.
    */
   SPHSimulation mySPHSimulation(
         ParticlesConfig::numberOfParticles, ParticlesConfig_bound::numberOfParticles,
         ParticlesConfig::searchRadius, ParticlesConfig::gridXsize * ParticlesConfig::gridYsize );

   /**
     * TEMP.
     */
   mySPHSimulation.model->h = SPHConfig::h;
   mySPHSimulation.model->m = SPHConfig::mass;
   mySPHSimulation.model->speedOfSound = SPHConfig::speedOfSound;
   mySPHSimulation.model->coefB = SPHConfig::coefB;
   mySPHSimulation.model->rho0 = SPHConfig::rho0;

   /**
    * Read the particle file.
    */
   #include "readParticleData.h"

   /**
    * Define timers to measure computation time.
    */
   TNL::Timer timer_search, timer_interact, timer_integrate, timer_pressure;
   TNL::Timer timer_search_reset, timer_search_cellIndices, timer_search_sort, timer_search_toCells;

   /**
    * TEMP: Determine number of interation for constant timestep.
    * Perform simulation main loop.
    */
   int steps = endTime / SPHConfig::dtInit;
   std::cout << "Number of steps: " << steps << std::endl;

   for( unsigned int iteration = 0; iteration < steps; iteration ++ )
   {
      std::cout << "STEP: " << iteration << std::endl;

      /**
       * Find neighbors within the SPH simulation.
       */
      timer_search.start();
      mySPHSimulation.PerformNeighborSearch(
            iteration, timer_search_reset, timer_search_cellIndices, timer_search_sort, timer_search_toCells );
      timer_search.stop();
      std::cout << "Search... done. " << std::endl;

      /**
       * Perform interaction with given model.
       */
      timer_interact.start();
      mySPHSimulation.template InteractModel< SPH::WendlandKernel, DiffusiveTerm, ViscousTerm, EOS >();
      timer_interact.stop();
      std::cout << "Interact... done. " << std::endl;

      //#include "outputForDebug.h"

      /**
       * Perform time integration, i.e. update particle positions.
       */
      timer_integrate.start();
      if( iteration % 20 == 0 ) {
         mySPHSimulation.integrator->IntegrateEuler( SPHConfig::dtInit );
         mySPHSimulation.integrator->IntegrateEulerBoundary( SPHConfig::dtInit );
      }
      else {
         mySPHSimulation.integrator->IntegrateVerlet( SPHConfig::dtInit );
         mySPHSimulation.integrator->IntegrateVerletBoundary( SPHConfig::dtInit );
      }
      timer_integrate.stop();
      std::cout << "integrate... done. " << std::endl;


      /**
       * Output particle data
       */
      if( ( iteration % outputStep ==  0) && (iteration > 0) )
      {
         /**
          * Compute pressure from density.
          * This is not necessary since we do this localy, if pressure is needed.
          * Its useful for output anywal.
          */
         timer_pressure.start();
         mySPHSimulation.model->template ComputePressureFromDensity< EOS >();
         timer_pressure.stop();
         std::cout << "Compute pressure... done. " << std::endl;

         outputFileName += std::to_string( iteration ) + ".vtk";
         std::ofstream outputFile (outputFileName, std::ofstream::out);
         Writer myWriter( outputFile, VTK::FileFormat::binary );
         myWriter.writeParticles( *mySPHSimulation.particles );
         myWriter.template writePointData< SPHModel::ScalarArrayType >(
               mySPHSimulation.model->FluidVariables.p, "Pressure" );
         myWriter.template writeVector< SPHModel::VectorArrayType, SPHConfig::RealType >(
               mySPHSimulation.model->FluidVariables.v, "Velocity", 3 );
      }
   }

   /**
    * Output simulation stats.
    */
   std::cout << std::endl << "COMPUTATION TIME:" << std::endl;
   std::cout << "Search........................................ " << timer_search.getRealTime() << " sec." << std::endl;
   std::cout << "Search (average time per step)................ " << timer_search.getRealTime() / steps << " sec." << std::endl;
   std::cout << " - Reset ..................................... " << timer_search_reset.getRealTime() << " sec." << std::endl;
   std::cout << " - Reset (average time per step).............. " << timer_search_reset.getRealTime() / steps << " sec." << std::endl;
   std::cout << " - Index by cell ............................. " << timer_search_cellIndices.getRealTime() << " sec." << std::endl;
   std::cout << " - Index by cell (average time per step)...... " << timer_search_cellIndices.getRealTime() / steps << " sec." << std::endl;
   std::cout << " - Sort ...................................... " << timer_search_sort.getRealTime() << " sec." << std::endl;
   std::cout << " - Sort (average time per step)............... " << timer_search_sort.getRealTime() / steps << " sec." << std::endl;
   std::cout << " - Particle to cell .......................... " << timer_search_toCells.getRealTime() << " sec." << std::endl;
   std::cout << " - Particle to cell (average time per step)... " << timer_search_toCells.getRealTime() / steps << " sec." << std::endl;
   std::cout << "Interaction................................... " << timer_interact.getRealTime() << " sec." << std::endl;
   std::cout << "Interaction (average time per step)........... " << timer_interact.getRealTime() / steps << " sec." << std::endl;
   std::cout << "Integrate..................................... " << timer_integrate.getRealTime() << " sec." << std::endl;
   std::cout << "Integrate (average time per step)............. " << timer_integrate.getRealTime() / steps << " sec." << std::endl;
   std::cout << "Pressure update............................... " << timer_pressure.getRealTime() << " sec." << std::endl;
   std::cout << "Pressure update (average time per step)....... " << timer_pressure.getRealTime() / steps << " sec." << std::endl;
   std::cout << "Total......................................... " << ( timer_search.getRealTime() + \
   + timer_interact.getRealTime() + timer_integrate.getRealTime() + timer_pressure.getRealTime() ) << " sec." << std::endl;
   std::cout << "Total (average time per step)................. " << ( timer_search.getRealTime() + \
   + timer_interact.getRealTime() + timer_integrate.getRealTime() + timer_pressure.getRealTime() ) / steps << " sec." << std::endl;


   std::cout << "\nDone ... " << std::endl;
}
