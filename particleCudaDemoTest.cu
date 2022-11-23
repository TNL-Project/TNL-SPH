#include <iostream>

#include <type_traits>  // std::add_const_t
#include <memory>       // std::shared_ptr

#include <TNL/Devices/Cuda.h>
#include <TNL/Pointers/DevicePointer.h>

#include "Particles/Particles.h"
#include "ParticlesConfig.h"

#include "Particles/neighbourSearch.h"


#ifdef HAVE_CUDA
template< typename ParticlePointer >
__global__ void testSetPoints( ParticlePointer gpu_particles )
{
	 if( threadIdx.x < gpu_particles->getNumberOfParticles() )
 	 {
 	     gpu_particles->setPoint( threadIdx.x, {1.5, 2.6} );
 	 }
}

template< typename ParticlePointer >
__global__ void testGetPoints( ParticlePointer gpu_particles )
{
	 if( threadIdx.x < gpu_particles->getNumberOfParticles() )
 	 {
 	     printf(" x: %f ", gpu_particles->getPoint( threadIdx.x )[ 0 ]);
 	 }
}

template< typename ParticlePointer, typename DeviceType >
void testParallelFor( ParticlePointer gpu_particles )
{
   auto init = [=] __cuda_callable__ ( int i ) mutable
   {
		 	printf(" nps: %d", gpu_particles->getNumberOfParticles() );
   };
   Algorithms::ParallelFor< DeviceType >::exec( 1, gpu_particles->getNumberOfParticles() - 1, init );
}
#endif /* HAVE_CUDA */


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
	using ParticleType = ParticleSystem::Particles< ParticleSystem::ParticleSystemConfig, Devices::Cuda >;
	using ParticlePointer = typename Pointers::SharedPointer< ParticleType, Devices::Cuda >;

  ParticlePointer myParticleSystem(nptcs, 1);


	std::cout << "Particle system dimension: " << myParticleSystem->getParticleDimension() << std::endl;
	std::cout << "Particle system points: " << myParticleSystem->getPoints() << std::endl;
	//you cant do this outside of Kernel: std::cout << "Particle system points - point with idx 1: " << myParticleSystem->getPoint(1) << std::endl;

#ifdef HAVE_CUDA

	testSetPoints< ParticlePointer ><<< 1, nptcs >>>( myParticleSystem );
  std::cout << "\nParticle system point modified: " << myParticleSystem->getPoints() << std::endl;

#endif /* HAVE_CUDA */

	/**
	 * Generate random particles.
	 */
	std::cout << "\nGenerate random particle positions." << std::endl;
	myParticleSystem->generateRandomParticles();
	//myParticleSystem.setPoint(1, {1, 1});
	std::cout << "Particle points: " << myParticleSystem->getPoints() << std::endl;

//test: #ifdef HAVE_CUDA
//test:
//test: 	std::cout << "\nTest getPoint in kernel: " << std::endl;
//test: 	testGetPoints< ParticlePointer ><<< 1, nptcs >>>( myParticleSystem );
//test: 	std::cout << std::endl;
//test:
//test: #endif /* HAVE_CUDA */

//test:	std::cout << "\nTest parallel for acces to the data:" << std::endl;
//test:	testParallelFor< ParticlePointer, Devices::Cuda >( myParticleSystem );
//test:	std::cout << std::endl;

	/**
	 * Show initial (empty) gird and particel cell indices.
	 */
	std::cout << "\nDefault (empty) grid and particle cell indices." << std::endl;
	std::cout << "Grid cell indices: " << myParticleSystem->getGridCellIndices() << std::endl;
	//you cant do this on Cuda: std::cout << "Grid cell indices - cell index with idx 1: " << myParticleSystem.getGridCellIndex(1) << std::endl;
	std::cout << "Particle cell indices: " << myParticleSystem->getParticleCellIndices() << std::endl;
	//you cant do this on Cuda: std::cout << "Particle cell - cell index with idx 1: " << myParticleSystem.getParticleCellIndex(1) << std::endl;

	/**
	 * Compute gird nad partice cell indices and show them.
	 */
	std::cout << "\nCompute grid cell indices." << std::endl;
	myParticleSystem->computeGridCellIndices();
	std::cout << "Grid cell indices: " << myParticleSystem->getGridCellIndices() << std::endl;

  std::cout << "\nCompute particle cell index." << std::endl;
  myParticleSystem->computeParticleCellIndices();
  std::cout << "Particle cell indices: " << myParticleSystem->getParticleCellIndices() << std::endl;

  //:/**
  //: * Print information about particle system.
  //: */
  //:std::cout << "\nInformation about particle system:" << std::endl;
  //:myParticleSystem.GetParticlesInformations();

  /**
   * Sort particles by its cell indices.
   * Show sorted indices and points sortet by indices.
   */
  std::cout << "\nSort particles. " << std::endl;
  myParticleSystem->sortParticles();
  std::cout << "SORTED particle cell indices: " << myParticleSystem->getParticleCellIndices() << std::endl;
  std::cout << "SORTED particle system points: " << myParticleSystem->getPoints() << std::endl;

  /**
   * Init neighbor search.
   * Show initial empty negihborSearch array with first and last particle in each cell.
   */
  std::cout << "\nNeighbor search." << std::endl;
	using NeighborSearchType = typename ParticleSystem::NeighborSearch< ParticleSystem::ParticleSystemConfig, ParticleType >;
	using NeighborSearchPointer = typename Pointers::SharedPointer< NeighborSearchType, Devices::Cuda >;

  //ParticleSystem::NeighborSearch< ParticleSystem::ParticleSystemConfig,
  //                                ParticleSystem::Particles< ParticleSystem::ParticleSystemConfig,
  //                                                           Devices::Cuda >> myNeighborSearch(myParticleSystem, 100);
	NeighborSearchPointer myNeighborSearch( myParticleSystem, 100 );

  std::cout << "First cell index: " << myNeighborSearch->getCellFirstParticleList() << std::endl;
  std::cout << "Last cell idnex: " << myNeighborSearch->getCellLastParticleList() << std::endl;

	/**
	 * Assign particles to neighborSearch arrays.
	 * Show filled negihborSearch array with first and last particle in each cell.
	 */
	myNeighborSearch->particlesToCells();
	std::cout << "First cell index after update: " << myNeighborSearch->getCellFirstParticleList() << std::endl;
	std::cout << "Last cell idnex after update: " << myNeighborSearch->getCellLastParticleList() << std::endl;


  /*** TEMP TEST ***/
  std::cout << "\nCycle over internal cells." << std::endl;
  myNeighborSearch->runCycleOverGrid();

///""
///""  //:/**
///""  //: * Print neihbor list.
///""  //: */
///""  //:#include "printNeighborList.h"


  std::cout << "\nDone ... " << std::endl;


}
