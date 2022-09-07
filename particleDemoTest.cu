#include <iostream>

#include <TNL/Devices/Host.h>
#include <TNL/Devices/Cuda.h>

#include "Particles/Particles.h"
#include "ParticlesConfig.h"


using namespace TNL;

int main( int argc, char* argv[] )
{

  /**
   * Number of particles
   */
  unsigned int nptcs = 32;

  /**
   * Create particle system
   * @number of particles
   * @search radius
   */
  ParticleSystem::Particles< ParticleSystem::ParticleSystemConfig, Devices::Host > myParticleSystem(nptcs, 1);

  std::cout << "Particle system dimension: " << myParticleSystem.getParticleDimension() << std::endl;
  std::cout << "Particle system points: " << myParticleSystem.getPoints() << std::endl;
  std::cout << "Particle system points - point with idx 1: " << myParticleSystem.getPoint(1) << std::endl;

  /**
   * Generate random particles.
   */
  std::cout << "\nGenerate random particle positions." << std::endl;
  myParticleSystem.generateRandomParticles();
  std::cout << "Particle points: " << myParticleSystem.getPoints() << std::endl;

  /**
   * Show initial (empty) gird and particel cell indices.
   */
  std::cout << "\nDefault (empty) grid and particle cell indices." << std::endl;
  std::cout << "Grid cell indices: " << myParticleSystem.getGridCellIndices() << std::endl;
  std::cout << "Grid cell indices - cell index with idx 1: " << myParticleSystem.getGridCellIndex(1) << std::endl;
  std::cout << "Particle cell indices: " << myParticleSystem.getParticleCellIndices() << std::endl;
  std::cout << "Particle cell - cell index with idx 1: " << myParticleSystem.getParticleCellIndex(1) << std::endl;

  /**
   * Compute gird nad partice cell indices and show them.
   */
  std::cout << "\nCompute grid cell indices." << std::endl;
  myParticleSystem.computeGridCellIndices();
  std::cout << "Grid cell indices: " << myParticleSystem.getGridCellIndices() << std::endl;

  std::cout << "\nCompute particle cell index." << std::endl;
  myParticleSystem.computeParticleCellIndices();
  std::cout << "Particle cell indices: " << myParticleSystem.getParticleCellIndices() << std::endl;

  /**
   * Print information about particle system.
   */
  std::cout << "\nInformation about particle system:" << std::endl;
  myParticleSystem.GetParticlesInformations();

  /**
   * Sort particles by its cell indices.
   * Show sorted indices and points sortet by indices.
   */
  std::cout << "\nSort particles. " << std::endl;
  myParticleSystem.sortParticles();
  std::cout << "SORTED particle cell indices: " << myParticleSystem.getParticleCellIndices() << std::endl;
  std::cout << "SORTED particle system points: " << myParticleSystem.getPoints() << std::endl;

  std::cout << "\nDone ... " << std::endl;


}

