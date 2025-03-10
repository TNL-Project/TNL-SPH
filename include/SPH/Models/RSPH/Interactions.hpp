#include "Interactions.h"

namespace TNL {
namespace SPH {

template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename BoudaryPointer >
void
RSPH< Particles, ModelConfig >::interaction( FluidPointer& fluid, BoudaryPointer& boundary, ModelParams& sphState )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   typename Particles::NeighborsLoopParams searchInFluid( fluid->getParticles() );
   typename Particles::NeighborsLoopParams searchInBound( boundary->getParticles() );

   /* CONSTANT VARIABLES */
   const RealType searchRadius = fluid->getParticles()->getSearchRadius();
   const RealType h = sphState.h;
   const RealType m = sphState.mass;
   const VectorType gravity = sphState.gravity;

   typename RiemannSolver::ParamsType riemannSolverParams( sphState );
   typename EOS::ParamsType eosParams( sphState );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->getParticles()->getPoints().getView();
   const auto view_rho = fluid->getVariables()->rho.getView();
   auto view_Drho = fluid->getVariables()->drho.getView();
   const auto view_v = fluid->getVariables()->v.getView();
   auto view_a = fluid->getVariables()->a.getView();

   const auto view_points_bound = boundary->getParticles()->getPoints().getView();
   const auto view_rho_bound = boundary->getVariables()->rho.getView();
   const auto view_v_bound = boundary->getVariables()->v.getView();
   const auto view_n_bound = boundary->getVariables()->n.getView();

   auto FluidFluid = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j,
         VectorType& r_i, VectorType& v_i, RealType& rho_i, RealType& p_i, RealType* drho_i, VectorType* a_i ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if (drs <= searchRadius )
      {
         const VectorType v_j = view_v[ j ];
         const RealType rho_j = view_rho[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );

         /* Riemann problem: */
         const VectorType e_ij = ( -1.f ) * r_ij / drs;
         const RealType vL = ( v_i, e_ij );
         const RealType vR = ( v_j, e_ij );

         const RealType v_starScalar = RiemannSolver::stateVelocity( vL, vR, p_i, p_j, 0.5f * ( rho_i + rho_j ), riemannSolverParams );
         const VectorType v_star = e_ij * v_starScalar + ( 0.5f * ( v_i + v_j ) - 0.5f * e_ij * ( vL + vR ) );
         const RealType p_star = RiemannSolver::statePressure( vL, vR, p_i, p_j, 0.5f * ( rho_i + rho_j ), riemannSolverParams );

         /* Interaction: */
         const VectorType v_is = v_i - v_star;

         const RealType F = KernelFunction::F( drs, h );
         const VectorType gradW = r_ij * F;

         *drho_i += 2.f * ( v_is, gradW ) * m;

         const RealType p_term = ( p_star ) / ( rho_i * rho_j );
         *a_i += ( -2.f ) * ( p_term ) * gradW * m;
      }
   };

   auto FluidBound = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j,
         VectorType& r_i, VectorType& v_i, RealType& rho_i, RealType& p_i, RealType* drho_i, VectorType* a_i ) mutable
   {
      const VectorType r_j = view_points_bound[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if (drs <= searchRadius )
      {
         const VectorType v_j = view_v_bound[ j ];
         const RealType rho_j = view_rho_bound[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );
         const VectorType n_j = view_n_bound[ j ];

         /* Riemann problem: */
         const VectorType e_ij = ( -1.f ) * r_ij / drs;
         const RealType vL = ( -1.f ) * ( v_i, n_j ); //normal points to the domain
         const RealType vR = -vL;
         const RealType pR = p_i + ( -1.f ) * rho_i * ( r_ij, gravity ); //sign

         const RealType v_starScalar = RiemannSolver::stateVelocity( vL, vR, p_i, pR, 0.5f * ( rho_i + rho_j ), riemannSolverParams );
         const VectorType v_star = e_ij * v_starScalar + ( 0.5f * ( v_i + v_j ) - 0.5f * e_ij * ( vL + vR ) );
         const RealType p_star = RiemannSolver::statePressure( vL, vR, p_i, pR, 0.5f * ( rho_i + rho_j ), riemannSolverParams );

         /* Interaction: */
         const VectorType v_is = v_i - v_star;

         const RealType F = KernelFunction::F( drs, h );
         const VectorType gradW = r_ij * F;

         *drho_i += 2.f * ( v_is, gradW ) * m;

         const RealType p_term = ( p_star ) / ( rho_i * rho_j );
         *a_i += ( -2.f ) * ( p_term ) * gradW * m;
      }
   };

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points[ i ];
      const VectorType v_i = view_v[ i ];
      const RealType rho_i = view_rho[ i ];
      const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );

      VectorType a_i = 0.f;
      RealType drho_i = 0.f;

      Particles::NeighborsLoop::exec( i, r_i, searchInFluid, FluidFluid, v_i, rho_i, p_i, &drho_i, &a_i );
      Particles::NeighborsLoop::exec( i, r_i, searchInBound, FluidBound, v_i, rho_i, p_i, &drho_i, &a_i );

      view_Drho[ i ] = drho_i;
      a_i += gravity;
      view_a[ i ] = a_i;
   };
   fluid->getParticles()->forAll( particleLoop );
}

template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename BoudaryPointer >
void
RSPH< Particles, ModelConfig >::updateSolidBoundary( FluidPointer& fluid, BoudaryPointer& boundary, ModelParams& sphState )
{}

template< typename Particles, typename ModelConfig >
template< typename EquationOfState, typename PhysicalObjectPointer >
void
RSPH< Particles, ModelConfig >::computePressureFromDensity( PhysicalObjectPointer& physicalObject, ModelParams& sphState )
{
   auto view_rho = physicalObject->getVariables()->rho.getView();
   auto view_p = physicalObject->getVariables()->p.getView();

   typename EOS::ParamsType eosParams( sphState );

   auto evalPressure = [=] __cuda_callable__ ( int i ) mutable
   {
      view_p[ i ] = EquationOfState::DensityToPressure( view_rho[ i ], eosParams );
   };
   physicalObject->getParticles()->forAll( evalPressure ); //TODO: forloop?
}

} // SPH
} // TNL
