#include "Interactions.h"
#include "../../customParallelFor.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Particles, typename SPHFluidConfig, typename Variables >
template< typename FluidPointer, typename BoudaryPointer, typename OpenBoundaryPointer, typename NeighborSearchPointer, typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm, typename EOS  >
void
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::Interaction( FluidPointer& fluid, BoudaryPointer& boundary, OpenBoundaryPointer& openBoundary )
{

   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   GlobalIndexType numberOfParticles = fluid->particles->getNumberOfParticles();
   GlobalIndexType numberOfParticles_bound = boundary->particles->getNumberOfParticles();
   GlobalIndexType numberOfParticles_inlet = openBoundary->particles->getNumberOfParticles();
   const RealType searchRadius = fluid->particles->getSearchRadius();

   const VectorType gridOrigin = fluid->particles->getGridOrigin();
   const IndexVectorType gridSize = fluid->particles->getGridSize();

   const auto view_firstLastCellParticle = fluid->neighborSearch->getCellFirstLastParticleList().getView();
   const auto view_particleCellIndex = fluid->particles->getParticleCellIndices().getView();
   const auto view_points = fluid->particles->getPoints().getView();

   const auto view_firstLastCellParticle_bound = boundary->neighborSearch->getCellFirstLastParticleList().getView();
   const auto view_particleCellIndex_bound = boundary->particles->getParticleCellIndices().getView();
   const auto view_points_bound = boundary->particles->getPoints().getView();

   const auto view_firstLastCellParticle_inlet = openBoundary->neighborSearch->getCellFirstLastParticleList().getView();
   const auto view_particleCellIndex_inlet = openBoundary->particles->getParticleCellIndices().getView();
   const auto view_points_inlet = openBoundary->particles->getPoints().getView();

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

   const auto view_rho_inlet = openBoundary->variables->rho.getView();
   const auto view_v_inlet = openBoundary->variables->v.getView();

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

         /* Interaction: */
         const VectorType v_ij = v_i - v_j;

         const RealType F = SPHKernelFunction::F( drs, h );
         const VectorType gradW = r_ij * F;

         const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, drs );
         const RealType diffTerm =  psi * ( r_ij, gradW ) * m / rho_j;
         *drho_i += ( v_ij, gradW ) * m - diffTerm;

         const RealType p_term = ( p_i + p_j ) / ( rho_i * rho_j );
         const RealType visco =  ViscousTerm::Pi( rho_i, rho_j, drs, ( r_ij, v_ij ) );
         *a_i += ( -1.0f ) * ( p_term + visco )* gradW * m;
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

         /* Interaction: */
         const VectorType v_ij = v_i - v_j;

         const RealType F = SPHKernelFunction::F( drs, h );
         const VectorType gradW = r_ij * F;

         const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, drs );
         const RealType diffTerm =  psi * ( r_ij, gradW ) * m / rho_j;
         *drho_i += ( v_ij, gradW ) * m - diffTerm;

         const RealType p_term = ( p_i + p_j ) / ( rho_i * rho_j );
         const RealType visco =  ViscousTerm::Pi( rho_i, rho_j, drs, ( r_ij, v_ij ) );
         *a_i += ( -1.0f ) * ( p_term + visco )* gradW * m;
      }
   };

   auto FluidInlet = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j, VectorType& r_i, VectorType& v_i, RealType& rho_i, RealType& p_i, RealType* drho_i, VectorType* a_i ) mutable
   {
      const VectorType r_j = view_points_inlet[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if (drs <= searchRadius )
      {
         const VectorType v_j = view_v_inlet[ j ];
         const RealType rho_j = view_rho_inlet[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j );

         /* Interaction: */
         const VectorType v_ij = v_i - v_j;

         const RealType F = SPHKernelFunction::F( drs, h );
         const VectorType gradW = r_ij * F;

         const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, drs );
         const RealType diffTerm =  psi * ( r_ij, gradW ) * m / rho_j;
         *drho_i += ( v_ij, gradW ) * m - diffTerm;

         const RealType p_term = ( p_i + p_j ) / ( rho_i * rho_j );
         const RealType visco =  ViscousTerm::Pi( rho_i, rho_j, drs, ( r_ij, v_ij ) );
         *a_i += ( -1.0f ) * ( p_term + visco )* gradW * m;
      }
   };

   auto BoundFluid = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j, VectorType& r_i, VectorType& v_i, RealType& rho_i, RealType& p_i, RealType* drho_i ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         const VectorType v_j = view_v[ j ];
         const RealType rho_j = view_rho[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j );

         /* Interaction */
         const VectorType v_ij = v_i - v_j;

         const RealType F = SPHKernelFunction::F( drs, h );
         const VectorType gradW = r_ij * F;

         const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, drs );
         const RealType diffTerm =  psi * ( r_ij, gradW ) * m / rho_j;
         *drho_i += ( v_ij, gradW ) * m - diffTerm;
      }
   };

   auto BoundInlet = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j, VectorType& r_i, VectorType& v_i, RealType& rho_i, RealType& p_i, RealType* drho_i ) mutable
   {
      const VectorType r_j = view_points_inlet[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         const VectorType v_j = view_v_inlet[ j ];
         const RealType rho_j = view_rho_inlet[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j );

         /* Interaction */
         const VectorType v_ij = v_i - v_j;

         const RealType F = SPHKernelFunction::F( drs, h );
         const VectorType gradW = r_ij * F;

         const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, drs );
         const RealType diffTerm =  psi * ( r_ij, gradW ) * m / rho_j;
         *drho_i += ( v_ij, gradW ) * m - diffTerm;
      }
   };

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i, NeighborSearchPointer& neighborSearch, NeighborSearchPointer& neighborSearch_bound, NeighborSearchPointer& neighborSearch_buffer ) mutable
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
      neighborSearch_buffer->loopOverNeighbors( i, numberOfParticles_inlet, gridIndex, gridSize, view_firstLastCellParticle_inlet, view_particleCellIndex, FluidInlet, r_i, v_i, rho_i, p_i, &drho_i, &a_i );

      view_Drho[ i ] = drho_i;
      //a_i[ 1 ] -= 9.81f ;
      a_i += gravity;
      view_a[ i ] = a_i;
   };
   SPHParallelFor::exec( 0, numberOfParticles, particleLoop, fluid->neighborSearch, boundary->neighborSearch, openBoundary->neighborSearch );

   auto particleLoopBoundary = [=] __cuda_callable__ ( LocalIndexType i, NeighborSearchPointer& neighborSearch, NeighborSearchPointer& neighborSearch_buffer ) mutable
   {
      const unsigned int activeCell = view_particleCellIndex[ i ];

      /*TODO: This should be some interaction structure  - properties of particle A:*/
      const VectorType r_i = view_points_bound[ i ];
      const VectorType v_i = view_v_bound[ i ];
      const RealType rho_i = view_rho_bound[ i ];
      const RealType p_i = EOS::DensityToPressure( rho_i );

      const IndexVectorType gridIndex = TNL::floor( ( r_i - gridOrigin ) / searchRadius );

      RealType drho_i = 0.;

      neighborSearch->loopOverNeighbors( i, numberOfParticles, gridIndex, gridSize, view_firstLastCellParticle, view_particleCellIndex_bound, BoundFluid, r_i, v_i, rho_i, p_i, &drho_i );
      neighborSearch_buffer->loopOverNeighbors( i, numberOfParticles_inlet, gridIndex, gridSize, view_firstLastCellParticle_inlet, view_particleCellIndex, BoundInlet, r_i, v_i, rho_i, p_i, &drho_i );

      view_Drho_bound[ i ] = drho_i;

   };
   SPHParallelFor::exec( 0, numberOfParticles_bound, particleLoopBoundary, fluid->neighborSearch, openBoundary->neighborSearch );
}


template< typename Particles, typename SPHFluidConfig, typename Variables >
template< typename EquationOfState >
void
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::ComputePressureFromDensity( VariablesPointer& variables, GlobalIndexType numberOfParticles )
{
   auto view_rho = variables->rho.getView();
   auto view_p = variables->p.getView();

   auto init = [=] __cuda_callable__ ( int i ) mutable
   {
      view_p[ i ] = EquationOfState::DensityToPressure( view_rho[ i ] );
   };
   Algorithms::ParallelFor< DeviceType >::exec( 0, numberOfParticles, init );
}

} // SPH
} // ParticleSystem
} // TNL

