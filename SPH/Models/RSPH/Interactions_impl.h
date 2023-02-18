#include "Interactions.h"
#include "../../customParallelFor.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Particles, typename SPHFluidConfig, typename Variables >
template< typename FluidPointer, typename BoudaryPointer, typename NeighborSearchPointer, typename RiemannSolver, typename EOS >
void
RSPHSimple< Particles, SPHFluidConfig, Variables >::Interaction( FluidPointer& fluid, BoudaryPointer& boundary )
{

   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   GlobalIndexType numberOfParticles = particles->getNumberOfParticles();
   GlobalIndexType numberOfParticles_bound = boundaryParticles->getNumberOfParticles();
   const RealType searchRadius = this->particles->getSearchRadius();

   const VectorType gridOrigin = fluid->particles->getGridOrigin();
   const IndexVectorType gridSize = fluid->particles->getGridSize();

   const auto view_firstLastCellParticle = fluid->neighborSearch->getCellFirstLastParticleList().getView();
   const auto view_particleCellIndex = fluid->particles->getParticleCellIndices().getView();
   const auto view_points = fluid->particles->getPoints().getView();

   const auto view_firstLastCellParticle_bound = boundary->neighborSearch->getCellFirstLastParticleList().getView();
   const auto view_particleCellIndex_bound = boundary->particles->getParticleCellIndices().getView();
   const auto view_points_bound = boundary->particles->getPoints().getView();

   /* CONSTANT VARIABLES */
   const RealType h = this->h;
   const RealType m = this->m;
   const RealType speedOfSound = this->speedOfSound;
   const RealType coefB = this->coefB;
   const RealType rho0 = this->rho0;
   const RealType delta = this->delta;
   const RealType alpha = this->alpha;
   const VectorType gravity = this->g;

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_rho = fluid->variables->rho.getView();
   auto view_Drho = fluid->variables->drho.getView();
   const auto view_v = fluid->variables->v.getView();
   auto view_a = fluid->variables->a.getView();

   const auto view_rho_bound = boundary->variables->rho.getView();
   auto view_Drho_bound = boundary->variables->drho.getView();
   const auto view_v_bound = boundary->variables->v.getView();

   auto FluidFluid = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j, VectorType& r_i, VectorType& v_i, RealType& rho_i, RealType& p_i, RealType* drho_i, VectorType* a_i ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if (drs <= searchRadius )
      {
         const VectorType v_j = view_v[ j ];
         const RealType rho_j = view_rho[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j );

         /* Riemann problem: */
         const VectorType e_ij = ( -1.f ) * r_ij / drs;
         const RealType vL = ( v_i, e_ij );
         const RealType vR = ( v_j, e_ij );

         const RealType v_starScalar = RiemannSolver::stateVelocity( vL, vR, p_i, p_j, 0.5f * ( rho_i + rho_j ) );
         const VectorType v_star = e_ij * v_starScalar + ( 0.5f * ( v_i + v_j ) - 0.5f * e_ij * ( vL + vR ) );
         const RealType p_star = RiemannSolver::statePressure( vL, vR, p_i, p_j, 0.5f * ( rho_i + rho_j ) );

         /* Interaction: */
         const VectorType v_is = v_i - v_star;

         const RealType F = SPHKernelFunction::F( drs, h );
         const VectorType gradW = r_ij * F;

         *drho_i += 2.f * ( v_is, gradW ) * m;

         const RealType p_term = ( p_star ) / ( rho_i * rho_j );
         *a_i += ( -2.f ) * ( p_term ) * gradW * m;
      }
   };

   auto FluidBound = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j, VectorType& r_i, VectorType& v_i, RealType& rho_i, RealType& p_i, RealType* drho_i, VectorType* a_i ) mutable
   {
      const VectorType r_j = view_points_bound[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if (drs <= searchRadius )
      {
         const VectorType v_j = view_v_bound[ j ];
         const RealType rho_j = view_rho_bound[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j );
         const VectorType n_j = view_n_bound[ j ];

         /* Riemann problem: */
         const VectorType e_ij = ( -1.f ) * r_ij / drs;
         const RealType vL = ( -1.f ) * ( v_i, n_j ); //normal points to the domain
         const RealType vR = -vL;
         const RealType pR = p_i + ( -1.f ) * rho_i * ( r_ij, gravity ); //sign

         const RealType v_starScalar = RiemannSolver::stateVelocity( vL, vR, p_i, pR, 0.5f * ( rho_i + rho_j ) );
         const VectorType v_star = e_ij * v_starScalar + ( 0.5f * ( v_i + v_j ) - 0.5f * e_ij * ( vL + vR ) );
         const RealType p_star = RiemannSolver::statePressure( vL, vR, p_i, pR, 0.5f * ( rho_i + rho_j ) );

         /* Interaction: */
         const VectorType v_is = v_i - v_star;

         const RealType F = SPHKernelFunction::F( drs, h );
         const VectorType gradW = r_ij * F;

         *drho_i += 2.f * ( v_is, gradW ) * m;

         const RealType p_term = ( p_star ) / ( rho_i * rho_j );
         *a_i += ( -2.f ) * ( p_term ) * gradW * m;
      }
   };

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i, NeighborSearchPointer& neighborSearch, NeighborSearchPointer& neighborSearch_bound ) mutable
   {
      const unsigned int activeCell = view_particleCellIndex[ i ];

      /*TODO: This should be some interaction structure  - properties of particle A:*/
      const VectorType r_i = view_points[ i ];
      const VectorType v_i = view_v[ i ];
      const RealType rho_i = view_rho[ i ];
      const RealType p_i = EOS::DensityToPressure( rho_i );

      const IndexVectorType gridIndex = TNL::floor( ( r_i - gridOrigin ) / searchRadius );

      VectorType a_i = 0.f;
      RealType drho_i = 0.f;

      neighborSearch->loopOverNeighbors( i, numberOfParticles, gridIndex, gridSize, view_firstLastCellParticle, view_particleCellIndex, FluidFluid, r_i, v_i, rho_i, p_i, &drho_i, &a_i );
      neighborSearch_bound->loopOverNeighbors( i, numberOfParticles_bound, gridIndex, gridSize, view_firstLastCellParticle_bound, view_particleCellIndex, FluidBound, r_i, v_i, rho_i, p_i, &drho_i, &a_i );

      view_Drho[ i ] = drho_i;
      a_i += gravity;
      view_a[ i ] = a_i;
   };
   SPHParallelFor::exec( 0, numberOfParticles, particleLoop, fluid->neighborSearch, boundary->neighborSearch );
}

} // SPH
} // ParticleSystem
} // TNL

