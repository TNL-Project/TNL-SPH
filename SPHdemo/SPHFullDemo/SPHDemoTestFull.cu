#include <iostream>

#include <TNL/Devices/Host.h>
#include <TNL/Devices/Cuda.h>

#include "../../Particles/Particles.h"

//#include "../../ParticlesConfig.h"
#include "ParticlesConfig.h"

#include "../../Particles/neighbourSearch.h"

#include "../../SPH/SPH.h"
#include "../../SPH/SPHFluidVariables.h"
#include "../../SPH/SPHConfig.h"
#include "../../SPH/SPHInteractions.h"

#include "../../Readers/VTKReader.h"

using namespace TNL;

int main( int argc, char* argv[] )
{

  /**
   * Number of particles
   */
  using ParticlesConfig = ParticleSystemConfig;
  using Device = Devices::Host;

  using ParticleSystem = typename ParticleSystem::Particles< ParticlesConfig, Device >;
  using Variables = typename TNL::ParticleSystem::SPH::SPHFluidVariables< TNL::ParticleSystem::SPH::SPHFluidConfig < Device > >;
  using NeighborSearch = typename TNL::ParticleSystem::NeighborSearch< ParticlesConfig, ParticleSystem >;

  using SPHModel = typename TNL::ParticleSystem::SPH::WCSPH_DBC< ParticleSystem, TNL::ParticleSystem::SPH::SPHFluidConfig < Device >>;
  using SPHSimulation = typename TNL::ParticleSystem::SPH::SPHSimulation< SPHModel, ParticleSystem, NeighborSearch >;

  using Reader = TNL::ParticleSystem::Readers::VTKReader;

  /**
   * Create the simulation.
   */
  SPHSimulation mySPHSimulation( ParticlesConfig::numberOfParticles, ParticlesConfig::searchRadius, ParticlesConfig::gridXsize * ParticlesConfig::gridYsize );

  /***
   * Read the particle file.
   */
  const std::string inputFileName="./CaseDamBreak_DEMO.vtk";
  Reader myReader( inputFileName );
  myReader.detectParticleSystem();
  myReader.template loadParticle< ParticleSystem >( mySPHSimulation.particles );

  /***
   * Setup type for boundary particles and initial condition.
   */
  for( unsigned int p = 0; p < mySPHSimulation.particles.getNumberOfParticles(); p ++ )
  {

    if( ( mySPHSimulation.particles.getPoint( p )[ 0 ] == 0. ) ||
        ( mySPHSimulation.particles.getPoint( p )[ 0 ] == 4. ) ||
        ( mySPHSimulation.particles.getPoint( p )[ 1 ] == 0. ) )
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

  //:: /**
  //::  * Test particle positions -> model points bridge.
  //::  */
  //:: std::cout << "\nTest particle positions -> model points bridge." << std::endl;
  //:: std::cout << mySPHSimulation.model.points << std::endl; //For some reason this works,...
  //:: std::cout << "Particle points bridge - point with idx 1: " << mySPHSimulation.model.points[1] << std::endl;


  std::cout << "\nDone ... " << std::endl;


}

