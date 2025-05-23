#include "SHTC.h"
#include "details.h"
#include <TNL/Containers/StaticVector.h>

namespace TNL {
namespace SPH {

template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename BoudaryPointer >
void
SHTC< Particles, ModelConfig >::interaction( FluidPointer& fluid, BoudaryPointer& boundary, ModelParams& modelParams )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   typename Particles::NeighborsLoopParams searchInFluid( fluid->getParticles() );
   typename Particles::NeighborsLoopParams searchInBound( boundary->getParticles() );

   /* CONSTANT VARIABLES */
   const RealType searchRadius = fluid->getParticles()->getSearchRadius();
   const RealType h = modelParams.h;
   const RealType m = modelParams.mass;
   const VectorType gravity = modelParams.gravity;

   typename Stress::ParamsType stressParams( modelParams );
   typename DiffusiveTerm::ParamsType diffusiveTermsParams( modelParams );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->getParticles()->getPoints().getView();
   const auto view_rho = fluid->getVariables()->rho.getView();
   auto view_drhodt = fluid->getVariables()->drhodt.getView();
   const auto view_v = fluid->getVariables()->v.getView();
   auto view_dvdt = fluid->getVariables()->dvdt.getView();
   const auto view_A = fluid->getVariables()->A.getView();
   auto view_dAdt = fluid->getVariables()->dAdt.getView();

   // boundary
   const auto view_points_boundary = boundary->getParticles()->getPoints().getView();
   const auto view_rho_boundary = boundary->getVariables()->rho.getView();
   const auto view_v_boundary = boundary->getVariables()->v.getView();
   const auto view_A_boundary = boundary->getVariables()->A.getView();

   auto interaction = [=] __cuda_callable__ ( LocalIndexType i,
                                              LocalIndexType j,
                                              VectorType& r_i,
                                              VectorType& v_i,
                                              RealType& rho_i,
                                              MatrixType& A_i,
                                              MatrixType& stress_i,
                                              RealType* drhodt_i,
                                              VectorType* dvdt_i,
                                              MatrixType* dAdt_i,
                                              RealType* dens ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         const VectorType v_j = view_v[ j ];
         const RealType rho_j = view_rho[ j ];
         const MatrixType stress_j = Stress::distortionToStress( view_A[ j ], rho_j, stressParams );

         // Interaction
         const RealType F = KernelFunction::F( drs, h );
         const VectorType gradW = r_ij * F;
         const VectorType v_ij = v_i - v_j;

         const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, r_ij, drs, diffusiveTermsParams );
         const RealType diffTerm =  psi * ( r_ij, gradW ) * m / rho_j;
         *drhodt_i += ( v_ij, gradW ) * m; // - diffTerm;

         *dvdt_i += ( -1.0f ) * ( stress_i / ( rho_i * rho_i ) + stress_j / ( rho_j * rho_j ) ) * gradW * m;

         *dAdt_i += Matrices::matmul( A_i, Matrices::cross( v_ij, gradW ) ) * m / rho_j;
         *dens += rho_j;
      }
   };

   auto interactionBoundary = [=] __cuda_callable__ ( LocalIndexType i,
                                                      LocalIndexType j,
                                                      VectorType& r_i,
                                                      VectorType& v_i,
                                                      RealType& rho_i,
                                                      MatrixType& A_i,
                                                      MatrixType& stress_i,
                                                      RealType* drhodt_i,
                                                      VectorType* dvdt_i,
                                                      MatrixType* dAdt_i,
                                                      RealType* dens ) mutable
   {
      const VectorType r_j = view_points_boundary[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         const VectorType v_j = view_v_boundary[ j ];
         const RealType rho_j = view_rho_boundary[ j ];
         const MatrixType stress_j = Stress::distortionToStress( view_A_boundary[ j ], rho_j, stressParams );

         // Interaction
         const RealType F = KernelFunction::F( drs, h );
         const VectorType gradW = r_ij * F;
         const VectorType v_ij = v_i - v_j;

         const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, r_ij, drs, diffusiveTermsParams );
         const RealType diffTerm =  psi * ( r_ij, gradW ) * m / rho_j;
         *drhodt_i += ( v_ij, gradW ) * m;// - diffTerm;

         *dvdt_i += ( -1.0f ) * ( stress_i / ( rho_i * rho_i ) + stress_j / ( rho_j * rho_j ) ) * gradW * m;

         *dAdt_i += Matrices::matmul( A_i, Matrices::cross( v_ij, gradW ) ) * m / rho_j;
         *dens += rho_j;
      }
   };

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points[ i ];
      const VectorType v_i = view_v[ i ];
      const RealType rho_i = view_rho[ i ];
      const MatrixType A_i = view_A[ i ];
      const MatrixType stress_i = Stress::distortionToStress( A_i, rho_i, stressParams );

      RealType drhodt_i = 0.f;
      VectorType dvdt_i = 0.f;
      MatrixType dAdt_i = 0.f;
      RealType count = 0;

      Particles::NeighborsLoop::exec( i, r_i, searchInFluid, interaction, v_i, rho_i, A_i, stress_i, &drhodt_i, &dvdt_i, &dAdt_i, &count );
      Particles::NeighborsLoop::exec( i, r_i, searchInBound, interactionBoundary, v_i, rho_i, A_i, stress_i, &drhodt_i, &dvdt_i, &dAdt_i, &count );

      view_drhodt[ i ] = drhodt_i;
      dvdt_i += gravity;
      view_dvdt[ i ] = dvdt_i;
      view_dAdt[ i ] = dAdt_i;
   };
   fluid->getParticles()->forAll( particleLoop );
}

template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename BoudaryPointer >
void
SHTC< Particles, ModelConfig >::updateSolidBoundary( FluidPointer& fluid, BoudaryPointer& boundary, ModelParams& modelParams )
{
   // neighborhood search tokens
   typename Particles::NeighborsLoopParams searchInFluid( fluid->getParticles() );
   typename Particles::NeighborsLoopParams searchInBound( boundary->getParticles() );

   // constants
   const RealType searchRadius = fluid->getParticles()->getSearchRadius();
   const RealType h = modelParams.h;
   const RealType m = modelParams.mass;

   typename Stress::ParamsType stressParams( modelParams );
   typename DiffusiveTerm::ParamsType diffusiveTermsParams( modelParams );

   // variables and arrays
   const auto view_points = fluid->getParticles()->getPoints().getView();
   const auto view_rho = fluid->getVariables()->rho.getView();
   const auto view_v = fluid->getVariables()->v.getView();
   const auto view_A = fluid->getVariables()->A.getView();

   // boundary
   const auto view_points_boundary = boundary->getParticles()->getPoints().getView();
   const auto view_rho_boundary = boundary->getVariables()->rho.getView();
   const auto view_v_boundary = boundary->getVariables()->v.getView();
   const auto view_A_boundary = boundary->getVariables()->A.getView();
   auto view_dAdt_boundary = boundary->getVariables()->dAdt.getView();

   auto interaction = [=] __cuda_callable__ ( LocalIndexType i,
                                              LocalIndexType j,
                                              VectorType& r_i,
                                              VectorType& v_i,
                                              RealType& rho_i,
                                              MatrixType& A_i,
                                              MatrixType& stress_i,
                                              RealType* drhodt_i,
                                              VectorType* dvdt_i,
                                              MatrixType* dAdt_i ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         const VectorType v_j = view_v[ j ];
         const RealType rho_j = view_rho[ j ];
         const MatrixType stress_j = Stress::distortionToStress( view_A[ j ], rho_j, stressParams );

         // Interaction
         const RealType F = KernelFunction::F( drs, h );
         const VectorType gradW = r_ij * F;
         const VectorType v_ij = v_i - v_j;

         const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, r_ij, drs, diffusiveTermsParams );
         const RealType diffTerm =  psi * ( r_ij, gradW ) * m / rho_j;
         *drhodt_i += ( v_ij, gradW ) * m; // - diffTerm;

         *dvdt_i += ( -1.0f ) * ( stress_i / ( rho_i * rho_i ) + stress_j / ( rho_j * rho_j ) ) * gradW * m;

         *dAdt_i += Matrices::matmul( A_i, Matrices::cross( v_ij, gradW ) ) * m / rho_j;
      }
   };

   auto interactionBoundary = [=] __cuda_callable__ ( LocalIndexType i,
                                                      LocalIndexType j,
                                                      VectorType& r_i,
                                                      VectorType& v_i,
                                                      RealType& rho_i,
                                                      MatrixType& A_i,
                                                      MatrixType& stress_i,
                                                      RealType* drhodt_i,
                                                      VectorType* dvdt_i,
                                                      MatrixType* dAdt_i ) mutable
   {
      const VectorType r_j = view_points_boundary[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         const VectorType v_j = view_v_boundary[ j ];
         const RealType rho_j = view_rho_boundary[ j ];
         const MatrixType stress_j = Stress::distortionToStress( view_A_boundary[ j ], rho_j, stressParams );

         // Interaction
         const RealType F = KernelFunction::F( drs, h );
         const VectorType gradW = r_ij * F;
         const VectorType v_ij = v_i - v_j;

         const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, r_ij, drs, diffusiveTermsParams );
         const RealType diffTerm =  psi * ( r_ij, gradW ) * m / rho_j;
         *drhodt_i += ( v_ij, gradW ) * m; // - diffTerm;

         *dvdt_i += ( -1.0f ) * ( stress_i / ( rho_i * rho_i ) + stress_j / ( rho_j * rho_j ) ) * gradW * m;

         *dAdt_i += Matrices::matmul( A_i, Matrices::cross( v_ij, gradW ) ) * m / rho_j;
      }
   };

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points_boundary[ i ];
      const VectorType v_i = view_v_boundary[ i ];
      const RealType rho_i = view_rho_boundary[ i ];
      const MatrixType A_i = view_A_boundary[ i ];
      const MatrixType stress_i = Stress::distortionToStress( A_i, rho_i, stressParams );

      //FIXME
      if ( v_i[ 0 ] > 0.1 ) return;

      RealType drhodt_i = 0.f;
      VectorType dvdt_i = 0.f;
      MatrixType dAdt_i = 0.f;

      Particles::NeighborsLoop::exec( i, r_i, searchInFluid, interaction, v_i, rho_i, A_i, stress_i, &drhodt_i, &dvdt_i, &dAdt_i );
      Particles::NeighborsLoop::exec( i, r_i, searchInBound, interactionBoundary, v_i, rho_i, A_i, stress_i, &drhodt_i, &dvdt_i, &dAdt_i );

      view_dAdt_boundary[ i ] = dAdt_i;
   };
   boundary->getParticles()->forAll( particleLoop );

}

template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename TimeStepping >
void
SHTC< Particles, ModelConfig >::relaxDistortion( FluidPointer& fluid, TimeStepping& timeStepping, ModelParams& modelParams )
{
   const RealType dt = timeStepping.getTimeStep();
   const RealType tau = modelParams.tau;
   auto A_view = fluid->getVariables()->A.getView();

   auto relax = [=] __cuda_callable__ ( const MatrixType& A )
   {
      const MatrixType fst = Matrices::matmul( Matrices::transpose( A ), A );
      return -3.f / tau * Matrices::matmul( A, details::deviator( fst ) );
   };

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const MatrixType A_i = A_view[ i ];
      MatrixType A_new_i = A_i;
      MatrixType K;

      K = relax( A_i );
      A_new_i += dt * K / 6.f;
      K = relax( A_i + dt * K / 2.f );
      A_new_i += dt * K / 3.f;
      K = relax( A_i + dt * K / 2.f );
      A_new_i += dt * K / 3.f;
      K = relax( A_i + dt * K );
      A_new_i += dt * K / 6.f;

      A_view[ i ] = A_new_i;
   };
   fluid->getParticles()->forAll( particleLoop );
}

} // SPH
} // TNL
