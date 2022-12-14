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

/**
 * Case configuration
 * One configuration for particle system, one for SPH.
 */
//#include "../testSetups/2769particles/ParticlesConfig.h"
//#include "../testSetups/2769particles/SPHCaseConfig.h"
//const std::string inputParticleFile = "../testSetups/2769particles/dambreak.vtk";

#include "../testSetups/49821particles/ParticlesConfig.h"
#include "../testSetups/49821particles/SPHCaseConfig.h"
const std::string inputParticleFile = "../testSetups/49821particles/dambreak.vtk";

//#include "../testSetups/189636particles/ParticlesConfig.h"
//#include "../testSetups/189636particles/SPHCaseConfig.h"
//const std::string inputParticleFile = "../testSetups/189636particles/dambreak.vtk";

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
  using ParticlesConfig = ParticleSystemConfig;
  using Device = ParticlesConfig::DeviceType;
  using SPHConfig = SPH::SPHCaseConfig< Device >;

  using ParticleSystem = typename ParticleSystem::Particles< ParticlesConfig, Device >;
  using Variables = typename TNL::ParticleSystem::SPH::SPHFluidVariables< SPHConfig >;
  using NeighborSearch = typename TNL::ParticleSystem::NeighborSearch< ParticlesConfig, ParticleSystem >;

  using SPHModel = typename TNL::ParticleSystem::SPH::WCSPH_DBC< ParticleSystem, SPHConfig >;
  using SPHSimulation = typename TNL::ParticleSystem::SPH::SPHSimulation< SPHModel, ParticleSystem, NeighborSearch >;

  using Reader = TNL::ParticleSystem::Readers::VTKReader;

  using DiffusiveTerm = TNL::ParticleSystem::SPH::MolteniDiffusiveTerm< SPHConfig >; //-> template
  using ViscousTerm = TNL::ParticleSystem::SPH::ArtificialViscosity< SPHConfig >; //-> template
  using EOS = TNL::ParticleSystem::SPH::TaitWeaklyCompressibleEOS< SPHConfig >; //move this inside model

	/**
	 * SPH solver models.
	 */

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
  double eps = 0.0001;
  #include "asignIinitialCondition.h"

  TNL::Timer timer;

  for( unsigned int time = 0; time < 2500; time ++ ) //2500
  {
    std::cout << "STEP: " << time << std::endl;

  	/**
  	 * Find neighbors within the SPH simulation.
  	 */
  	timer.reset();
		timer.start();
    mySPHSimulation.PerformNeighborSearch( time );
		timer.stop();
    std::cout << "Search... done. ";
		std::cout << "Total time: " << timer.getRealTime() << " sec." << std::endl;

  	timer.reset();
		timer.start();
  	mySPHSimulation.template Interact< SPH::WendlandKernel, DiffusiveTerm, ViscousTerm >();
		timer.stop();
    std::cout << "Interact... done. ";
		std::cout << "Total time: " << timer.getRealTime() << " sec." << std::endl;

    //#include "outputForDebug.h"

  	timer.reset();
		timer.start();
		if( time % 20 == 0 ) {
    	 mySPHSimulation.model->IntegrateEuler( 0.00002 ); //0.00005/0.000025/0.00001
		}
		else {
    	 mySPHSimulation.model->IntegrateVerlet( 0.00002 );
		}
		timer.stop();
    std::cout << "integrate... done. " ;
		std::cout << "Total time: " << timer.getRealTime() << " sec." << std::endl;

  	timer.reset();
		timer.start();
    mySPHSimulation.model->template ComputePressureFromDensity< EOS >();
		timer.stop();
    std::cout << "Compute pressure... done. ";
		std::cout << "Total time: " << timer.getRealTime() << " sec." << std::endl;

		if( ( time % 2500 ==  0) && (time > 1) )
		{
			std::string outputFileName = "results/particles";
			outputFileName += std::to_string(time) + ".ptcs";
			#include "writeParticleData.h"
		}
  }


	std::string outputFileName = "particles.ptcs";
	#include "writeParticleData.h"

  std::cout << "\nDone ... " << std::endl;
}
