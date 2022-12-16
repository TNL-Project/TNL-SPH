#include "SPH.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Variables, typename ParticleSystem, typename NeighborSearch >
void
SPHSimulation< Variables, ParticleSystem, NeighborSearch >::PerformNeighborSearch( GlobalIndexType step, TNL::Timer& timer_reset, TNL::Timer& timer_cellIndices, TNL::Timer& timer_sort, TNL::Timer& timer_toCells )
{
   /**
    * Compute gird nad partice cell indices.
    */
   timer_reset.start();
   neighborSearch->resetListWithIndices();
   timer_reset.stop();
   std::cout << " - neighborSearch->resetListWithIndices();... done" << std::endl;

   if( step == 0 ) //TODO: do this better
   particles->computeGridCellIndices();

   timer_cellIndices.start();
   particles->computeParticleCellIndices();
   timer_cellIndices.stop();
   std::cout << " - particles->computeParticleCellIndices();... done " << std::endl;

   timer_sort.start();
   model->sortParticlesAndVariables(); //particles.sortParticles();
   timer_sort.stop();
   std::cout << " - model->sortParticlesAndVariables();... done " << std::endl;

   timer_toCells.start();
   neighborSearch->particlesToCells();
   timer_toCells.stop();
   std::cout << " - neighborSearch->particlesToCells();... done " << std::endl;
}

template< typename Variables, typename ParticleSystem, typename NeighborSearch >
template< typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm >
void
SPHSimulation< Variables, ParticleSystem, NeighborSearch >::Interact()
{

   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   GlobalIndexType numberOfParticles = particles->getNumberOfParticles();
   static constexpr GlobalIndexType _numberOfCells = ParticleSystem::Config::gridXsize; //FIXIT
   const auto view_firstCellParticle = neighborSearch->getCellFirstParticleList().getView();
   const auto view_lastCellParticle = neighborSearch->getCellLastParticleList().getView(); // DEBUG
   const auto view_particleCellIndex = particles->getParticleCellIndices().getView();
   const auto view_points = particles->getPoints().getView();
   const RealType searchRadius = this->particles->getSearchRadius();

   /* CONSTANT VARIABLES */
   const RealType h = model->h;
   const RealType m = model->m;
   const RealType speedOfSound = model->speedOfSound;
   const RealType coefB = model->coefB;
   const RealType rho0 = model->rho0;
   const RealType delta = model->delta;
   const RealType alpha = model->alpha;

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

      /*TODO: This should be some interaction structure  - properties of particle A:*/
      const PointType r_i = view_points[ i ];
      const PointType v_i = view_v[ i ];
      const RealType rho_i = view_rho[ i ];
      const RealType p_i = view_p[ i ];

      PointType a_i = {0., 0.};
      RealType drho_i = 0.;

      /* LOAD OTHER PARTICLE DATA */
      // Process fluid particle
      if( view_particleType[ i ] == 0 )
      {
         for( int ci = -1; ci <= 1; ci++ ){
            for( int cj = -1; cj <= 1; cj++ ){
               const unsigned int neighborCell = activeCell + cj * _numberOfCells + ci;
               int j = view_firstCellParticle[ neighborCell ];

               //also works: int j_end = view_lastCellParticle[ neighborCell ]; //test;
               //also works: if( j_end >= numberOfParticles )
               //also works:  	j_end = -1;
               //also works: while( ( j <= j_end ) ){ //test
               while( ( j < numberOfParticles ) && ( view_particleCellIndex[ j ] == neighborCell ) ){
                  if( i == j ){ j++; continue; }

                  /* START OF LOOP OVER NEIGHBROS */

                  /* This should be some interaction structure - properties of particle B:*/
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
         a_i[ 1 ] -= 9.81 ;
         view_a[ i ] = a_i;

      } // if - process fluid particle
      // Process boundary particle
      else if( view_particleType[ i ] == 1 )
      {
         for( int ci = -1; ci <= 1; ci++ ){
            for( int cj = -1; cj <= 1; cj++ ){
               const unsigned int neighborCell = activeCell + cj * _numberOfCells + ci;
               int j = view_firstCellParticle[ neighborCell ];

               //also works: int j_end = view_lastCellParticle[ neighborCell ];//test;
               //also works: if( j_end >= numberOfParticles )
               //also works:  	j_end = -1;
               //also works: while( ( j <= j_end ) ){ //test
               while( ( j < numberOfParticles ) && ( view_particleCellIndex[ j ] == neighborCell ) ){
                  if( i == j ){ j++; continue; }

                  /* START OF LOOP OVER NEIGHBROS */

                  if( view_particleType[ j ] == 0 ){

                  /* This should be some interaction structure, mby. - properties of particle B: */
                  const PointType r_j = view_points[ j ];
                  const PointType v_j = view_v[ j ];
                  const RealType rho_j = view_rho[ j ];
                  const RealType p_j = view_p[ j ];

                  /* Interaction */
                  const PointType dr = r_i - r_j;
                  const PointType dv = v_i - v_j;

                  const RealType drs = l2Norm( dr );
                  const RealType F = SPHKernelFunction::F( drs, h );
                  const PointType gradW = dr*F;

                  const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, drs );
                  const RealType diffTerm =  psi * ( dr, gradW ) * m / rho_j;
                  drho_i += ( dv, gradW ) * m - diffTerm;
               }

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
         printf( "INVALID PARTICLE TYPE!" ); //cerr
      }
   };
   Algorithms::ParallelFor< DeviceType >::exec( 0, numberOfParticles, particleLoop );
}

template< typename Variables, typename ParticleSystem, typename NeighborSearch >
template< typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm >
void
SPHSimulation< Variables, ParticleSystem, NeighborSearch >::InteractModel()
{
   model->template Interaction< NeighborSearchPointer, SPHKernelFunction, DiffusiveTerm, ViscousTerm >( neighborSearch );
}

} // SPH
} // ParticleSystem
} // TNL

