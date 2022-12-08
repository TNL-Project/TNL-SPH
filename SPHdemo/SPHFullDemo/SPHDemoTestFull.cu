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
#include "ParticlesConfig.h"
#include "SPHCaseConfig.h"

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
  const std::string inputFileName="./dambreak.vtk";
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

  for( unsigned int time = 0; time < 2; time ++ )
  //for( unsigned int time = 0; time < 200; time ++ )
  {

    std::cout << "STEP: " << time << std::endl;
    mySPHSimulation.PerformNeighborSearch( time );
    std::cout << "search... done." << std::endl;

  	mySPHSimulation.template Interact< SPH::WendlandKernel, DiffusiveTerm, ViscousTerm >();
    std::cout << "interact... done." << std::endl;

    //#include "outputForDebug.h"

		if( time % 20 == 0 ) {
    	//mySPHSimulation.model.integrator.IntegrateEuler( ParticlesConfig::numberOfParticles, 0.00005 );
    	mySPHSimulation.model->IntegrateEuler( 0.00005 );
		}
		else {
    	//mySPHSimulation.model.integrator.IntegrateVerlet( ParticlesConfig::numberOfParticles, 0.00005 );
    	mySPHSimulation.model->IntegrateVerlet( 0.00005 );
		}
		//if(time < 5 )
   	//mySPHSimulation.model.integrator.IntegrateVerlet( ParticlesConfig::numberOfParticles, 0.00005 );

		//: mySPHSimulation.model.integrator.IntegrateEuler( ParticlesConfig::numberOfParticles, 0.00005 );
    std::cout << "integrate... done." << std::endl;

    mySPHSimulation.model->template ComputePressureFromDensity< EOS >();
    std::cout << "compute pressure... done." << std::endl;

		//#include "outputForDebugNbs.h"

    //std::cout << "mySPHSimulation POINTS: " << time << std::endl << mySPHSimulation.particles.getPoints() << std::endl;
    //std::cout << "mySPHSimulation DERIVATIVES: " << time << std::endl << mySPHSimulation.model.vars.DrhoDv << std::endl;
    //std::cout << "mySPHSimulation PRESSURE: " << mySPHSimulation.model->getPress() << std::endl;
    std::cout << "mySPHSimulation DENSIY: " << mySPHSimulation.model->getDrho() << std::endl;
    //std::cout << "mySPHSimulation RHO: " << mySPHSimulation.model->getAcc() << std::endl;
  }

	//#include "writeBoundaryParticleData.h"

  //#include "printNeighborList.h"
	//#include "writeFluidParticleData.h"
	#include "writeParticleData.h"

	//mySPHSimulation.particles.GetParticlesInformations();


	/*
  std::ofstream file( "grid.vti" );
	using Grid = typename ParticleSystem::GridType;
  using Writer = TNL::Meshes::Writers::VTUWriter< Grid >;
  Writer writer( file );
  writer.writeImageData( *mySPHSimulation.particles.grid );
	*/

  //std::cout << "mySPHSimulation points after integration step: " << time << std::endl << mySPHSimulation.particles.getPoints() << std::endl;
  //std::cout << "mySPHSimulation interaction values after integration step: " << mySPHSimulation.model.vars.rho << std::endl;

  //std::cout << "mySPHSimulation points after integration: " << mySPHSimulation.particles.getPoints() << std::endl;

  //:: /**
  //::  * Test particle positions -> model points bridge.
  //::  */
  //:: std::cout << "\nTest particle positions -> model points bridge." << std::endl;
  //:: std::cout << mySPHSimulation.model.points << std::endl; //For some reason this works,...
  //:: std::cout << "Particle points bridge - point with idx 1: " << mySPHSimulation.model.points[1] << std::endl;


  std::cout << "\nDone ... " << std::endl;


}
