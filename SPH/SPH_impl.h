#include "SPH.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Variables, typename ParticleSystem, typename NeighborSearch >
void
SPHSimulation< Variables, ParticleSystem, NeighborSearch >::PerformNeighborSearch( GlobalIndexType step )
{
   /**
    * Compute gird nad partice cell indices.
    */
   particles.resetNeighborList();
	 if( step == 0 )
   particles.computeGridCellIndices(); //I DONT NEED TO REPEAT THIS!

   std::cout << "SPHSimulation::PerformNeighborSearch(): ... OK" << std::endl; //debug
   particles.computeParticleCellIndices();
   std::cout << "SPHSimulation::computeParticleCellIndices(): ... OK" << std::endl; //debug

	 if( step % 1 == 0 )
	 {
    	model.sortParticlesAndVariables(); //I DONT NEED TO DO THIS IN EACH STEP!
    	std::cout << "SPHSimulation::sortParticlesAndVariables(): ... OK" << std::endl; //debug
    	//particles.sortParticles();
	 }

   //std::cout << particles.getParticleCellIndices() << std::endl; //debug

   /**
    * Find neigbors.
    */
   neighborSearch.searchForNeighbors();
   std::cout << "neighborSearch.searchForNeighbors() ... OK" << std::endl; //debug
}

template< typename Variables, typename ParticleSystem, typename NeighborSearch >
void
SPHSimulation< Variables, ParticleSystem, NeighborSearch >::Interact()
{

	 /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
	 GlobalIndexType numberOfParticles = particles->getNumberOfParticles();
   static constexpr GlobalIndexType _numberOfCells = ParticleSystem::ParticleConfig::gridXsize; //FIXIT
	 const auto view_firstCellParticle = neighborSearch.getCellFirstParticleList().getView();
	 const auto view_particleCellIndex = particles.getParticleCellIndices().getView();
	 const auto view_points = particles.getPoints().getView();
	 RealType searchRadius = this->particles->getSearchRadius();

	 /* VARIABLES AND FIELD ARRAYS */

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i  ) mutable
	 {
		 const unsigned int activeCell = view_particleCellIndex[ i ];

		 /* LOAD OTHER PARTICLE DATA */

		 for( int ci = -1; ci <= 1; ci++ ){
			 for( int cj = -1; cj <= 1; cj++ ){

				 const unsigned int neighborCell = activeCell + cj * _numberOfCells + ci;
				 int j = view_firstCellParticle[ neighborCell ]; //USE INT MAX

				 while( ( j < numberOfParticles ) && ( j >= 0 ) && ( view_particleCellIndex[ j ] == neighborCell ) ){

       			//if( ( l2Norm( view_points[ i ] - view_points[ j ] ) < searchRadius ) && ( i != j ) )
			 			//{	}

					 	/* START OF LOOP OVER NEIGHBROS */


					 	/* END OF LOOP OVER NEIGHBROS */

					 j++;

				 } //while over particle in cell
			 } //for cells in y direction
		 } //for cells in x direction

	 /* SAVE INTERACTION RESULTS */

	 };
	 Algorithms::ParallelFor< DeviceType >::exec( 0, numberOfParticles, particleLoop );
}

} // SPH
} // ParticleSystem
} // TNL

