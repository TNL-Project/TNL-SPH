#include <iostream>
#include <fstream> //temp, to write output

#include <TNL/Devices/Cuda.h>
#include <string>

/**
 * Particle system.
 */
#include "../../Particles/Particles.h"
#include "../../Particles/neighbourSearch.h"

/**
 * Particle system reader.
 **/
#include "../../Readers/VTKReader.h"
//#include "../../Writers/VTKWriter.h"

/**
 * Case configuration
 * One configuration for particle system, one for SPH.
 */
//#include "../SPHCaseSetup/damBreak2769particles/ParticlesConfig.h"
//#include "../SPHCaseSetup/damBreak2769particles/SPHCaseConfig.h"
//const std::string inputParticleFile = "../SPHCaseSetup/damBreak2769particles/dambreak.vtk";

#include "../SPHCaseSetup/damBreak49821particles/ParticlesConfig.h"
#include "../SPHCaseSetup/damBreak49821particles/SPHCaseConfig.h"
const std::string inputParticleFile = "../SPHCaseSetup/damBreak49821particles/dambreak.vtk";

//#include "../SPHCaseSetup/damBreakN49206particles/ParticlesConfig.h"
//#include "../SPHCaseSetup/damBreakN49206particles/SPHCaseConfig.h"
//const std::string inputParticleFile = "../SPHCaseSetup/damBreakN49206particles/dambreak.vtk";

//#include "../SPHCaseSetup/damBreak189636particles/ParticlesConfig.h"
//#include "../SPHCaseSetup/damBreak189636particles/SPHCaseConfig.h"
//const std::string inputParticleFile = "../SPHCaseSetup/damBreak189636particles/dambreak.vtk";

//#include "../SPHCaseSetup/damBreakN736806particles/ParticlesConfig.h"
//#include "../SPHCaseSetup/damBreakN736806particles/SPHCaseConfig.h"
//const std::string inputParticleFile = "../SPHCaseSetup/damBreakN736806particles/dambreak.vtk";

//#include "../SPHCaseSetup/damBreakN2913606particles/ParticlesConfig.h"
//#include "../SPHCaseSetup/damBreakN2913606particles/SPHCaseConfig.h"
//const std::string inputParticleFile = "../SPHCaseSetup/damBreakN2913606particles/dambreak.vtk";

//const float endTime = 0.0002;
const float endTime = 0.05;

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
   using SPHConfig = SPH::SPHCaseConfig< Device >;

   using ParticleSystem = typename ParticleSystem::Particles< ParticlesConfig, Device >;
   using Variables = typename TNL::ParticleSystem::SPH::SPHFluidVariables< SPHConfig >;
   using NeighborSearch = typename TNL::ParticleSystem::NeighborSearch< ParticlesConfig, ParticleSystem >;

   using SPHModel = typename TNL::ParticleSystem::SPH::WCSPH_DBC< ParticleSystem, SPHConfig >;
   using SPHSimulation = typename TNL::ParticleSystem::SPH::SPHSimulation< SPHModel, ParticleSystem, NeighborSearch >;

   using Reader = TNL::ParticleSystem::Readers::VTKReader;
   //using Writer = TNL::ParticleSystem::Writers::VTKWriter< ParticleSystem >;

   /**
    * SPH solver models.
    */
   using DiffusiveTerm = TNL::ParticleSystem::SPH::MolteniDiffusiveTerm< SPHConfig >; //-> template
   using ViscousTerm = TNL::ParticleSystem::SPH::ArtificialViscosity< SPHConfig >; //-> template
   using EOS = TNL::ParticleSystem::SPH::TaitWeaklyCompressibleEOS< SPHConfig >; //move this inside model

   /**
    * Create the simulation.
    */
   SPHSimulation mySPHSimulation( ParticlesConfig::numberOfParticles, ParticlesConfig::searchRadius, ParticlesConfig::gridXsize * ParticlesConfig::gridYsize );

   /**
     * Temporary.
     */
   mySPHSimulation.model->h = SPHConfig::h;
   mySPHSimulation.model->m = SPHConfig::mass;
   mySPHSimulation.model->speedOfSound = SPHConfig::speedOfSound;
   mySPHSimulation.model->coefB = SPHConfig::coefB;
   mySPHSimulation.model->rho0 = SPHConfig::rho0;

   /**
    * Read the particle file.
    */
   //Create temp particle system to read data in:
   using ParticleSystemToReadData = typename ParticleSystem::Particles< ParticlesConfig, Devices::Host >;
   ParticleSystemToReadData particlesToRead( ParticlesConfig::numberOfParticles, ParticlesConfig::searchRadius );

   const std::string inputFileName = inputParticleFile;
   Reader myReader( inputFileName );
   myReader.detectParticleSystem();
   //myReader.template loadParticle< ParticleSystem >( *mySPHSimulation.particles );
   myReader.template loadParticle< ParticleSystemToReadData >( particlesToRead );

   /**
    * Setup type for boundary particles and initial condition.
    */
   #include "asignIinitialCondition.h"

   //const std::string outputFileName = "output.vtk";
   //std::ofstream outputFile (outputFileName, std::ofstream::out);
   //Writer myWriter( outputFile, VTK::FileFormat::ascii );
   ////myWriter.template loadParticle< ParticleSystem >( *mySPHSimulation.particles );
   //myWriter.writeParticles( *mySPHSimulation.particles );
   ////myWriter.template writeMetadata( 1, 0 );
   //myWriter.template writePointData< SPHModel::ScalarArrayType >( mySPHSimulation.model->getRho(), "Density" );
   ////myWriter.template writePointData< SPHModel::VectorArrayType >( mySPHSimulation.model->getVel(), "Velocity" );

   TNL::Timer timer_search, timer_interact, timer_integrate, timer_pressure;
   TNL::Timer timer_search_reset, timer_search_cellIndices, timer_search_sort, timer_search_toCells;

   int steps = endTime / SPHConfig::dtInit;
   //int steps = 1800;
   std::cout << "Number of steps: " << steps << std::endl;
   for( unsigned int time = 0; time < steps; time ++ ) //2500
   {
      std::cout << "STEP: " << time << std::endl;

      /**
       * Find neighbors within the SPH simulation.
       */
      timer_search.start();
      mySPHSimulation.PerformNeighborSearch( time, timer_search_reset, timer_search_cellIndices, timer_search_sort, timer_search_toCells );
      timer_search.stop();
      std::cout << "Search... done. " << std::endl;

      timer_interact.start();
      mySPHSimulation.template InteractModel< SPH::WendlandKernel, DiffusiveTerm, ViscousTerm, EOS >();
      timer_interact.stop();
      std::cout << "Interact... done. " << std::endl;

      //#include "outputForDebug.h"

      timer_integrate.start();
      if( time % 20 == 0 ) {
         mySPHSimulation.model->IntegrateEuler( SPHConfig::dtInit ); //0.00005/0.00002/0.00001
      }
      else {
         mySPHSimulation.model->IntegrateVerlet( SPHConfig::dtInit );
      }
      timer_integrate.stop();
      std::cout << "integrate... done. " << std::endl;

      timer_pressure.start();
      mySPHSimulation.model->template ComputePressureFromDensity< EOS >();
      timer_pressure.stop();
      std::cout << "Compute pressure... done. " << std::endl;

      if( ( time % 2500 ==  0) && (time > 1) )
      {
         std::string outputFileName = "results/particles";
         outputFileName += std::to_string(time) + ".ptcs";
         #include "writeParticleData.h"
      }
   }

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


   //:std::string outputFileName = "particles.ptcs";
   //:#include "writeParticleData.h"

   //:std::cout << "\nDone ... " << std::endl;
}
