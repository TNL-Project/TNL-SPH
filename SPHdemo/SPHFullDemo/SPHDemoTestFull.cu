#include <iostream>
#include <fstream> //temp, to write output
#include <TNL/Meshes/Writers/VTIWriter.h> //temp,t o write grid

//#include <TNL/Devices/Host.h>
#include <TNL/Devices/Cuda.h>

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
//#include "testSetups/2769particles/ParticlesConfig.h"
//#include "testSetups/2769particles/SPHCaseConfig.h"
//const std::string inputParticleFile = "testSetups/2769particles/dambreak.vtk";

#include "testSetups/49821particles/ParticlesConfig.h"
#include "testSetups/49821particles/SPHCaseConfig.h"
const std::string inputParticleFile = "testSetups/49821particles/dambreak.vtk";

//#include "testSetups/189636particles/ParticlesConfig.h"
//#include "testSetups/189636particles/SPHCaseConfig.h"
//const std::string inputParticleFile = "testSetups/189636particles/dambreak.vtk";

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

  /**
   * Create the simulation.
   */
  SPHSimulation mySPHSimulation( ParticlesConfig::numberOfParticles, ParticlesConfig::searchRadius, ParticlesConfig::gridXsize * ParticlesConfig::gridYsize );

  /*** TEMP ***/
  mySPHSimulation.model->h = SPHConfig::h;
  mySPHSimulation.model->m = SPHConfig::mass;
  mySPHSimulation.model->speedOfSound = SPHConfig::speedOfSound;
  mySPHSimulation.model->coefB = SPHConfig::coefB;
  mySPHSimulation.model->rho0 = SPHConfig::rho0;

  /***
   * Read the particle file.
   */
  const std::string inputFileName = inputParticleFile;
  Reader myReader( inputFileName );
  myReader.detectParticleSystem();
  myReader.template loadParticle< ParticleSystem >( *mySPHSimulation.particles );

  /***
   * Setup type for boundary particles and initial condition.
   */
  double eps = 0.0001;
  #include "asignIinitialCondition.h"

  //std::cout << "mySPHSimulation particle positions: " << mySPHSimulation.model.points << std::endl;
  //:info: std::cout << "Grid informations: " << std::endl;
  //:info: mySPHSimulation.particles.GetParticlesInformations();

  //:debug: /**
  //:debug:  * Find neighbors within the SPH simulation - step by step.
  //:debug:  */
  //:debug: mySPHSimulation.particles.computeGridCellIndices();
  //:debug: mySPHSimulation. particles.computeParticleCellIndices();
  //:debug: //model.sortParticlesAndVariables();
  //:debug: mySPHSimulation.particles.sortParticles();
  //:debug: std::cout << "Grid cell indices: " << mySPHSimulation.particles.getGridCellIndices() << std::endl;
  //:debug std::cout << "Particle cell indices: " << mySPHSimulation.particles.getParticleCellIndices() << std::endl;

  using DiffusiveTerm = TNL::ParticleSystem::SPH::MolteniDiffusiveTerm< SPHConfig >; //-> template
  using ViscousTerm = TNL::ParticleSystem::SPH::ArtificialViscosity< SPHConfig >; //-> template

  //test: /**
  //test:  * Find neighbors within the SPH simulation.
  //test:  */
  //test: std::cout << "\nFind neighbros within the SPH simulation." << std::endl;
  //test: mySPHSimulation.PerformNeighborSearch();

  //test: /**
  //test:  * Test the loop over particle neighbors.
  //test:  */
  //test: std::cout << "\nTest the loop over particle neighbros." << std::endl;
  //test: mySPHSimulation.Interact();
  //test: std::cout << std::endl;
  //test: //std::cout << "mySPHSimulation interaction values: " << mySPHSimulation.model.vars.DrhoDv << std::endl;


  //test: double eps = 0.0001;
  //test: #include "outputForDebug.h"


  //test: std::cout << "\nTest one step of integration." << std::endl;
  //test: mySPHSimulation.model.integrator.IntegrateVerlet( ParticlesConfig::numberOfParticles, 0.00005 );
  //test: //mySPHSimulation.model.integrator.IntegrateVerlet( ParticlesConfig::numberOfParticles, 0.0001 );

  //std::cout << "mySPHSimulation points after integration: " << mySPHSimulation.particles.getPoints() << std::endl;
  //std::cout << "mySPHSimulation interaction values after integration step: " << mySPHSimulation.model.vars.rho << std::endl;


  using EOS = TNL::ParticleSystem::SPH::TaitWeaklyCompressibleEOS< SPHConfig >; //move this inside model

  TNL::Timer timer;

  for( unsigned int time = 0; time < 200; time ++ )
  //for( unsigned int time = 0; time < 200; time ++ )
  {

    std::cout << "STEP: " << time << std::endl;
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
    	mySPHSimulation.model->IntegrateEuler( 0.00005 );
		}
		else {
    	mySPHSimulation.model->IntegrateVerlet( 0.00005 );
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

		//#include "outputForDebug.h"
  }

	#include "writeParticleData.h"

  std::cout << "\nDone ... " << std::endl;
}
