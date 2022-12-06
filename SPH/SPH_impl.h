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
   //particles->resetNeighborList();
	 if( step == 0 )
   particles->computeGridCellIndices(); //I DONT NEED TO REPEAT THIS!

   std::cout << "SPHSimulation::PerformNeighborSearch(): ... OK" << std::endl; //debug
   particles->computeParticleCellIndices();
   std::cout << "SPHSimulation::computeParticleCellIndices(): ... OK" << std::endl; //debug

	 if( step % 1 == 0 )
	 {
    	model->sortParticlesAndVariables(); //I DONT NEED TO DO THIS IN EACH STEP!
    	std::cout << "SPHSimulation::sortParticlesAndVariables(): ... OK" << std::endl; //debug
    	//particles.sortParticles();
	 }

   //debug: std::cout << particles->getParticleCellIndices() << std::endl; //debug

   //i dont need to compose nblist, interaction are calculated directly: /**
   //i dont need to compose nblist, interaction are calculated directly:  * Find neigbors.
   //i dont need to compose nblist, interaction are calculated directly:  */
   //i dont need to compose nblist, interaction are calculated directly: neighborSearch->searchForNeighbors();
   //i dont need to compose nblist, interaction are calculated directly: std::cout << "neighborSearch.searchForNeighbors() ... OK" << std::endl; //debug
}

template< typename Variables, typename ParticleSystem, typename NeighborSearch >
template< typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm >
//template< typename SPHKernelFunction = WendlandKernel, typename DiffusiveTerm = DiffusiveTerm_MT, typename VisousTerm = ViscousTerm_AV >
void
SPHSimulation< Variables, ParticleSystem, NeighborSearch >::Interact()
{

	 /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
	 GlobalIndexType numberOfParticles = particles->getNumberOfParticles();
   static constexpr GlobalIndexType _numberOfCells = ParticleSystem::Config::gridXsize; //FIXIT
	 const auto view_firstCellParticle = neighborSearch->getCellFirstParticleList().getView();
	 const auto view_particleCellIndex = particles->getParticleCellIndices().getView();
	 const auto view_points = particles->getPoints().getView();
	 RealType searchRadius = this->particles->getSearchRadius();

	 /* CONSTANT VARIABLES */
	 RealType h = model->h;
	 RealType m = model->m;
	 RealType speedOfSound = model->speedOfSound;
	 RealType coefB = model->coefB;
	 RealType rho0 = model->rho0;
	 RealType delta = model->delta;
	 RealType alpha = model->alpha;


	 /* VARIABLES AND FIELD ARRAYS */
	 const auto view_particleType = model->getParticleType().getView();
	 const auto view_rho = model->getRho().getView();
	 auto view_Drho = model->getDrho().getView();
	 const auto view_p = model->getPress().getView();
	 const auto view_v = model->getVel().getView();
	 auto view_a = model->getAcc().getView();

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i  ) mutable
	 {
		 const unsigned int activeCell = view_particleCellIndex[ i ];

		 /* This should be some interaction structure  - properties of particle A:*/
		 const PointType r_i = view_points[ i ];
		 const PointType v_i = view_v[ i ];
		 const RealType rho_i = view_rho[ i ];
		 const RealType p_i = view_p[ i ];

		 PointType a_i = {0., 0.};
		 RealType drho_i = 0.;

		 /* LOAD OTHER PARTICLE DATA */
		 /* Process fluid particle */
		 if( view_particleType[ i ] == 0 )
		 {
		 		for( int ci = -1; ci <= 1; ci++ ){
		 		  for( int cj = -1; cj <= 1; cj++ ){

		 		 	 const unsigned int neighborCell = activeCell + cj * _numberOfCells + ci;
		 		 	 int j = view_firstCellParticle[ neighborCell ]; //USE INT MAX

		 		 	 while( ( j < numberOfParticles ) && ( j >= 0 ) && ( view_particleCellIndex[ j ] == neighborCell ) ){

     		  		//if( ( l2Norm( view_points[ i ] - view_points[ j ] ) < searchRadius ) && ( i != j ) )
		 		  		//{	}

						 	// If im right, this is same for boundary as well as for fluid particles, sice we use DBC only for now.

		 		 		 	/* START OF LOOP OVER NEIGHBROS */
						 	/* This should be some interaction structure, mby. - properties of particle B: */
		 					const PointType r_j = view_points[ j ];
		 					const PointType v_j = view_v[ j ];
		 					const RealType rho_j = view_rho[ j ];
		 					const RealType p_j = view_p[ j ];

							/* Interaction: */

  						const PointType dr = r_i - r_j;
  						const PointType dv = v_i - v_j;

  						const RealType drs = l2Norm( dr );
  						const RealType F = SPHKernelFunction::F( drs, h );
  						const PointType gradW = dr * F;

							const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, drs );
							const RealType diffTerm =  psi * ( dr, gradW ) * m / rho_j;
  						drho_i += ( dv, gradW ) * m - diffTerm;
  						//const RealType drho = ( dv, gradW ) * vars.m + DiffusiveTerm::Psi( vars.rho[ i ], vars.rho[ j ], drs );

  						const RealType p_term = ( p_i + p_j ) / ( rho_i * rho_j );
  						const RealType visco =  ViscousTerm::Pi( rho_i, rho_j, drs, ( dr, dv ) );
  						a_i += ( -1.0 ) * ( p_term + visco )* gradW * m;


		 		 		 	/* END OF LOOP OVER NEIGHBROS */

		 		 		 j++;

		 		 	 } //while over particle in cell
		 		  } //for cells in y direction
		 		} //for cells in x direction

	 		/* SAVE INTERACTION RESULTS */

			view_Drho[ i ] = drho_i;
			a_i[1] -= 9.81;
			view_a[ i ] = a_i;


		 } // if - process fluid particle
		/* Process boundary particle */
		else if( view_particleType[ i ] == 1 )
		{
		 		for( int ci = -1; ci <= 1; ci++ ){
		 		  for( int cj = -1; cj <= 1; cj++ ){

		 		 	 const unsigned int neighborCell = activeCell + cj * _numberOfCells + ci;
		 		 	 int j = view_firstCellParticle[ neighborCell ]; //USE INT MAX

		 		 	 while( ( j < numberOfParticles ) && ( j >= 0 ) && ( view_particleCellIndex[ j ] == neighborCell ) ){

     		  			//if( ( l2Norm( view_points[ i ] - view_points[ j ] ) < searchRadius ) && ( i != j ) )
		 		  			//{	}

		 		 		 	/* START OF LOOP OVER NEIGHBROS */

						 	/* This should be some interaction structure, mby. - properties of particle B: */
		 					const PointType r_j = view_points[ j ];
		 					const PointType v_j = view_v[ j ];
		 					const RealType rho_j = view_rho[ j ];
		 					const RealType p_j = view_p[ j ];

							/* INteraction */

  						const PointType dr = r_i - r_j;
  						const PointType dv = v_i - v_j;

  						const RealType drs = l2Norm( dr );
  						const RealType F = SPHKernelFunction::F( drs, h );
  						const PointType gradW = dr*F;

							const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, drs );
							const RealType diffTerm =  psi * ( dr, gradW ) * m / rho_j;
  						const RealType drho = ( dv, gradW ) * m - diffTerm;
  						//const RealType drho = ( dv, gradW )*vars.m;
  						const PointType a = { 0., 0. };

		 		 		 	/* END OF LOOP OVER NEIGHBROS */

		 		 		 j++;

		 		 	 } //while over particle in cell
		 		  } //for cells in y direction
		 		} //for cells in x direction

			view_Drho[ i ] = drho_i;
			a_i = { 0., 0. };
			view_a[ i ] = a_i;

		} //else if - process boundary particle
		else
		{
			//cerr
			printf(" INVALID PARTICLE TYPE!");
		}


	 };
	 Algorithms::ParallelFor< DeviceType >::exec( 0, numberOfParticles, particleLoop );
}

} // SPH
} // ParticleSystem
} // TNL

