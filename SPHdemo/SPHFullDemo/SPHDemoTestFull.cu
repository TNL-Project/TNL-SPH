#include <iostream>

#include <TNL/Devices/Host.h>
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

using namespace TNL;

int main( int argc, char* argv[] )
{

  /**
   * Number of particles
   */
  using ParticlesConfig = ParticleSystemConfig;
  using Device = Devices::Host;
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
  mySPHSimulation.model.vars.h = ParticlesConfig::searchRadius;
  mySPHSimulation.model.vars.m = SPHConfig::mass;
  mySPHSimulation.model.vars.speedOfSound = SPHConfig::speedOfSound;
  mySPHSimulation.model.vars.coefB = SPHConfig::coefB;
  mySPHSimulation.model.vars.rho0 = SPHConfig::rho0;

  /***
   * Read the particle file.
   */
  const std::string inputFileName="./dambreak.vtk";
  Reader myReader( inputFileName );
  myReader.detectParticleSystem();
  myReader.template loadParticle< ParticleSystem >( mySPHSimulation.particles );

  /***
   * Setup type for boundary particles and initial condition.
   */
  for( unsigned int p = 0; p < mySPHSimulation.particles.getNumberOfParticles(); p ++ )
  {

    if( ( mySPHSimulation.particles.getPoint( p )[ 0 ] == 0. ) ||
        ( mySPHSimulation.particles.getPoint( p )[ 0 ] == -0.025 ) ||
        ( mySPHSimulation.particles.getPoint( p )[ 0 ] == -0.05 ) ||
        ( mySPHSimulation.particles.getPoint( p )[ 1 ] == 0. ) ||
        ( mySPHSimulation.particles.getPoint( p )[ 1 ] == -0.025 ) ||
        ( mySPHSimulation.particles.getPoint( p )[ 1 ] == -0.05 ) ||
        ( mySPHSimulation.particles.getPoint( p )[ 1 ] == 3.9 ) ||
        ( mySPHSimulation.particles.getPoint( p )[ 1 ] == 3.875 ) ||
        ( mySPHSimulation.particles.getPoint( p )[ 1 ] == 3.85 ) )
    {
      mySPHSimulation.model.vars.type[ p ] = 1.;
    }
    else
    {
      mySPHSimulation.model.vars.type[ p ] = 0.;
    }

      mySPHSimulation.model.vars.rho[ p ] = 1000.;
      mySPHSimulation.model.vars.rho[ p ] = 1000.;
      mySPHSimulation.model.vars.v[ p ] = 0.;

      //fill in integrator arrays
      mySPHSimulation.model.integrator.rhoO[ p ] = 1000.;
      mySPHSimulation.model.integrator.rhoOO[ p ] = 1000.;

      mySPHSimulation.model.integrator.vO[ p ] = 0.;
      mySPHSimulation.model.integrator.vOO[ p ] = 0.;

  }

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


  /**
   * Find neighbors within the SPH simulation.
   */
  std::cout << "\nFind neighbros within the SPH simulation." << std::endl;
  mySPHSimulation.PerformNeighborSearch();

  /**
   * Test the loop over particle neighbors.
   */
  std::cout << "\nTest the loop over particle neighbros." << std::endl;
  mySPHSimulation.Interact();
  std::cout << std::endl;
  std::cout << "mySPHSimulation interaction values: " << mySPHSimulation.model.vars.DrhoDv << std::endl;

  std::cout << "\nTest one step of integration." << std::endl;
  mySPHSimulation.model.integrator.IntegrateVerlet( ParticlesConfig::numberOfParticles, 0.0001 );
  //mySPHSimulation.model.integrator.IntegrateVerlet( ParticlesConfig::numberOfParticles, 0.0001 );

  std::cout << "mySPHSimulation points after integration: " << mySPHSimulation.particles.getPoints() << std::endl;
  std::cout << "mySPHSimulation interaction values after integration step: " << mySPHSimulation.model.vars.rho << std::endl;


  using EOS = TNL::ParticleSystem::SPH::TaitWeaklyCompressibleEOS; //move this inside model

  for( unsigned int time = 0; time < 500; time ++ )
  {

    mySPHSimulation.PerformNeighborSearch();
    std::cout << "search... done." << std::endl;

    mySPHSimulation.Interact();
    std::cout << "interact... done." << std::endl;

    mySPHSimulation.model.integrator.IntegrateVerlet( ParticlesConfig::numberOfParticles, 0.00001 );
    std::cout << "integrate... done." << std::endl;

    mySPHSimulation.model.template ComputePressureFromDensity< EOS >();
    std::cout << "compute pressure... done." << std::endl;


    std::cout << "mySPHSimulation POINTS: " << time << std::endl << mySPHSimulation.particles.getPoints() << std::endl;
    std::cout << "mySPHSimulation DERIVATIVES: " << time << std::endl << mySPHSimulation.model.vars.DrhoDv << std::endl;
    std::cout << "mySPHSimulation DENSIY: " << mySPHSimulation.model.vars.p << std::endl;
    std::cout << "mySPHSimulation RHO: " << mySPHSimulation.model.vars.v << std::endl;
  }

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

