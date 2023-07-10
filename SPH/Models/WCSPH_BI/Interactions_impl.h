#include "Interactions.h"
#include "../../customParallelFor.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename ParticleSystem, typename SPHFluidConfig, typename Variables >
template< typename FluidPointer,
          typename BoudaryPointer,
          typename SPHKernelFunction,
          typename DiffusiveTerm,
          typename ViscousTerm,
          typename EOS,
          typename SPHState >
void
WCSPH_BI< ParticleSystem, SPHFluidConfig, Variables >::Interaction( FluidPointer& fluid,
                                                                    BoudaryPointer& boundary,
                                                                    SPHState& sphState )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   GlobalIndexType numberOfParticles = fluid->particles->getNumberOfParticles();
   GlobalIndexType numberOfParticles_bound = boundary->particles->getNumberOfParticles();
   const RealType searchRadius = fluid->particles->getSearchRadius();

   typename ParticleSystem::NeighborsLoopParams searchInFluid( fluid->particles );
   typename ParticleSystem::NeighborsLoopParams searchInBound( boundary->particles );

   /* CONSTANT VARIABLES */
   const RealType h = sphState.h;
   const RealType m = sphState.mass;
   const RealType ds = sphState.boundaryElementSize;
   const RealType rho0 = sphState.rho0;
   const VectorType gravity = sphState.gravity;

   typename DiffusiveTerm::ParamsType diffusiveTermsParams( sphState );
   typename ViscousTerm::ParamsType viscousTermsParams( sphState );
   typename EOS::ParamsType eosParams( sphState );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->particles->getPoints().getView();
   const auto view_rho = fluid->variables->rho.getView();
   auto view_Drho = fluid->variables->drho.getView();
   const auto view_v = fluid->variables->v.getView();
   auto view_a = fluid->variables->a.getView();

   const auto view_points_bound = boundary->particles->getPoints().getView();
   auto view_rho_bound = boundary->variables->rho.getView();
   const auto view_v_bound = boundary->variables->v.getView();
   const auto view_n_bound = boundary->variables->n.getView();

   auto FluidFluid = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j,
         VectorType& r_i, VectorType& v_i, RealType& rho_i, RealType& p_i, RealType* drho_i, VectorType* a_i, RealType* gamma_i ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if (drs <= searchRadius )
      {
         const VectorType v_j = view_v[ j ];
         const RealType rho_j = view_rho[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );

         /* Interaction: */
         const VectorType v_ij = v_i - v_j;

         const RealType F = SPHKernelFunction::F( drs, h );
         const RealType W = SPHKernelFunction::W( drs, h );
         const VectorType gradW = r_ij * F;

         const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, drs, diffusiveTermsParams );
         const RealType diffTerm =  psi * ( r_ij, gradW ) * m / rho_j;
         *drho_i += ( v_ij, gradW ) * m - diffTerm;

         const RealType p_term = ( p_i + p_j ) / ( rho_i * rho_j );
         const RealType visco =  ViscousTerm::Pi( rho_i, rho_j, drs, ( r_ij, v_ij ), viscousTermsParams );
         *a_i += ( -1.0f ) * ( p_term + visco )* gradW * m;

         *gamma_i += W * m / rho_j;
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

         /* Interaction: */
         const VectorType v_ij = v_i - v_j;

         const RealType W = SPHKernelFunction::W( drs, h );

         *drho_i += ( -1.f ) * ( v_ij, n_j ) * W * rho_j * ds;

         const RealType p_term = ( p_i + p_j ) / ( rho_i * rho_j );
         const RealType visco =  ViscousTerm::Pi( rho_i, rho_j, drs, ( r_ij, v_ij ), viscousTermsParams );
         *a_i += ( p_term + visco ) * n_j * W * rho_j * ds;
      }
   };

   auto BoundFluid = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j,
         VectorType& r_i, RealType* rho_i, RealType* gamma_i ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         const RealType rho_j = view_rho[ j ];

         const RealType W = SPHKernelFunction::W( drs, h );

         *rho_i += W * m;
         *gamma_i += W * m / rho_j;
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
      RealType gamma_i = 0.f;

      NeighborsLoop::exec( i, r_i, searchInFluid, FluidFluid, v_i, rho_i, p_i, &drho_i, &a_i, &gamma_i );
      NeighborsLoopAnotherSet::exec( i, r_i, searchInBound, FluidBound, v_i, rho_i, p_i, &drho_i, &a_i );

      if( gamma_i > 0.01f ) {
         view_Drho[ i ] = drho_i / gamma_i;
         view_a[ i ] = a_i / gamma_i + gravity;
      }
      else {
         view_Drho[ i ] = 0.f;
         view_a[ i ] = 0.f + gravity;
      }
   };
   SPHParallelFor::exec( 0, numberOfParticles, particleLoop );

   auto particleLoopBoundary = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points_bound[ i ];

      RealType rho_i = 0.f;
      RealType gamma_i = 0.f;

      NeighborsLoopAnotherSet::exec( i, r_i, searchInFluid, BoundFluid, &rho_i, &gamma_i );

      if( gamma_i > 0.01f ){
         view_rho_bound[ i ] = ( rho_i / gamma_i > rho0 ) ? ( rho_i / gamma_i ) : rho0;
		}
      else{
         view_rho_bound[ i ] = rho0;
      }
   };
   SPHParallelFor::exec( 0, numberOfParticles_bound, particleLoopBoundary );
}

template< typename ParticleSystem, typename SPHFluidConfig, typename Variables >
template< typename EquationOfState, typename PhysicalObjectPointer, typename SPHState >
void
WCSPH_BI< ParticleSystem, SPHFluidConfig, Variables >::ComputePressureFromDensity( PhysicalObjectPointer& physicalObject, SPHState& sphState )
{
   auto view_rho = physicalObject->getVariables()->rho.getView();
   auto view_p = physicalObject->getVariables()->p.getView();

   typename EOS::ParamsType eosParams( sphState );

   auto init = [=] __cuda_callable__ ( int i ) mutable
   {
      view_p[ i ] = EquationOfState::DensityToPressure( view_rho[ i ], eosParams );
   };
   Algorithms::parallelFor< DeviceType >( physicalObject->getFirstActiveParticle(), physicalObject->getLastActiveParticle() + 1, init );
}

} // SPH
} // ParticleSystem
} // TNL

