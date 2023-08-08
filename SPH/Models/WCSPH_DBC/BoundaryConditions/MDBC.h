#include "../Interactions.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename ParticleSystem, typename SPHState >
template< typename FluidPointer,
          typename BoudaryPointer,
          typename SPHKernelFunction,
          typename DiffusiveTerm,
          typename ViscousTerm,
          typename EOS,
          typename BCType,
          typename std::enable_if_t< std::is_same_v< BCType, WCSPH_BCTypes::MDBC >, bool > Enabled >
void
WCSPH_DBC< ParticleSystem, SPHState >::updateSolidBoundary( FluidPointer& fluid,
                                                            BoudaryPointer& boundary,
                                                            SPHState& sphState )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   typename ParticleSystem::NeighborsLoopParams searchInFluid( fluid->particles );
   typename ParticleSystem::NeighborsLoopParams searchInBound( boundary->particles );

   /* CONSTANT VARIABLES */
   const RealType searchRadius = fluid->particles->getSearchRadius();
   const RealType h = sphState.h;
   const RealType m = sphState.mass;
   const RealType rho0 = sphState.rho0;

   typename EOS::ParamsType eosParams( sphState );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->particles->getPoints().getView();
   const auto view_rho = fluid->variables->rho.getView();
   const auto view_v = fluid->variables->v.getView();

   const auto view_points_bound = boundary->particles->getPoints().getView();
   auto view_rho_bound = boundary->variables->rho.getView();
   const auto view_v_bound = boundary->variables->v.getView();
   const auto view_ghostNode_bound = boundary->variables->ghostNodes.getView();

   auto BoundFluid = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j,
         VectorType& ghostNode_i, VectorType& v_i, RealType& rho_i, RealType& p_i, Matrix* A_gn, VectorExtendedType* b_gn ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_j - ghostNode_i; //FLUID_POS - GHOSTNODE_POS
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         const VectorType v_j = view_v[ j ];
         const RealType rho_j = view_rho[ j ];

         /* Interaction */
         const VectorType v_ij = v_i - v_j;

         const RealType F = SPHKernelFunction::F( drs, h );
         const RealType W = SPHKernelFunction::W( drs, h );
         const VectorType gradW = r_ij * F;

         const RealType V = m / rho_j;

         Matrix A_local = 0.f;
         VectorExtendedType b_local = 0.f;

         if( SPHConfig::spaceDimension == 2 )
         {
            A_local( 0, 0 ) = W * V;          A_local( 0, 1 ) = r_ij[ 0 ] * W * V;          A_local( 0, 2 ) = r_ij[ 1 ] * W * V;
            A_local( 1, 0 ) = gradW[ 0 ] * V; A_local( 1, 1 ) = r_ij[ 0 ] * gradW[ 0 ] * V; A_local( 1, 2 ) = r_ij[ 1 ] * gradW[ 0 ] * V;
            A_local( 2, 0 ) = gradW[ 1 ] * V; A_local( 2, 1 ) = r_ij[ 0 ] * gradW[ 1 ] * V; A_local( 2, 2 ) = r_ij[ 1 ] * gradW[ 1 ] * V;

            b_local[ 0 ] = W * m; b_local[ 1 ] = gradW[ 0 ] * m; b_local[ 2 ] = gradW[ 1 ] * m;

            *A_gn += A_local;
            *b_gn += b_local;
         }
      }
   };

   auto particleLoopBoundary = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points_bound[ i ];
      const VectorType v_i = view_v_bound[ i ];
      const RealType rho_i = view_rho_bound[ i ];
      const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );
      const VectorType ghostNode_i = view_ghostNode_bound[ i ];

      RealType drho_i = 0.f;
      Matrix A_gn = 0.f;
      VectorExtendedType b_gn = 0.f;

      NeighborsLoop::exec( i, ghostNode_i, searchInFluid, BoundFluid, v_i, rho_i, p_i, &A_gn, &b_gn );

      RealType rho_bound = 0.f;

      if( Matrices::determinant( A_gn ) > 0.001 )
      {
         VectorExtendedType rhoGradRho = Matrices::solve( A_gn, b_gn );
         VectorType r_ign = ghostNode_i - r_i;
         rho_bound = rhoGradRho[ 0 ] + rhoGradRho[ 1 ] * r_ign[ 0 ] + rhoGradRho[ 2 ] * r_ign[ 1 ];
      }
      else if( A_gn( 0, 0 ) > 0.f )
      {
         rho_bound = b_gn[ 0 ] / A_gn( 0, 0 );
      }
      else
      {
         rho_bound = rho0;
      }

      if( rho_bound < rho0 )
         rho_bound = rho0;

      view_rho_bound[ i ] = rho_bound;

   };
   SPHParallelFor::exec( boundary->getFirstActiveParticle(), boundary->getLastActiveParticle() + 1, particleLoopBoundary );

}

} // SPH
} // ParticleSystem
} // TNL

