#include <iostream>

#include <TNL/Devices/Host.h>
#include <TNL/Devices/Cuda.h>

#include "Particles/Particles.h"
#include "ParticlesConfig.h"

#include "Particles/neighbourSearch.h"

#include "./Readers/VTKReader.h"


using namespace TNL;

int main( int argc, char* argv[] )
{

  /**
   * Number of particles
   */
  //unsigned int nptcs = 2447; //particleSample_binary.vtk
  unsigned int nptcs = 16; //particleSample_simpleAscii.vtk

  /**
   * Create particle system
   * @number of particles
   * @search radius
   */
  using ParticleType = ParticleSystem::Particles< ParticleSystem::ParticleSystemConfig, Devices::Host >;

  ParticleSystem::Particles< ParticleSystem::ParticleSystemConfig, Devices::Host > myParticleSystem(nptcs, 1);

  std::cout << "Particle system dimension: " << myParticleSystem.getParticleDimension() << std::endl;
  std::cout << "Particle system points: " << myParticleSystem.getPoints() << std::endl;
  std::cout << "Particle system points - point with idx 1: " << myParticleSystem.getPoint(1) << std::endl;

  /**
   * Define VTKReader
   */
  std::cout << "\nRead particles." << std::endl;
  //const std::string inputFileName="/home/tomas/Documents/codes/tnl/myTNLProject/playground/particleSample_binary.vtk";
  const std::string inputFileName="/home/tomas/Documents/codes/tnl/myTNLProject/playground/particleSample_simpleAscii.vtk";
  ParticleSystem::Readers::VTKReader myVTKReader(inputFileName);

  myVTKReader.detectParticleSystem();
  myVTKReader.template loadParticle< ParticleType >( myParticleSystem );

  std::cout << "Particle system dimension: " << myParticleSystem.getParticleDimension() << std::endl;
  std::cout << "Particle system points: " << myParticleSystem.getPoints() << std::endl;
  std::cout << "Particle system points - point with idx 1: " << myParticleSystem.getPoint(1) << std::endl;

  std::cout << "\nreadParticleDEmoTest.cu: Done ... " << std::endl;


}

