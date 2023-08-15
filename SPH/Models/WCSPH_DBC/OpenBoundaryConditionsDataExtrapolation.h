#include "Interactions.h"
#include "details.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename ParticleSystem, typename SPHState >
template< typename FluidPointer,
          typename OpenBoudaryPointer,
          typename SPHKernelFunction,
          typename EOS >
void
WCSPH_DBC< ParticleSystem, SPHState >::extrapolateOpenBoundaryData( FluidPointer& fluid,
                                                                    OpenBoudaryPointer& openBoundary,
                                                                    SPHState& sphState )
{
   if constexpr( SPHState::SPHConfig::spaceDimension == 2 )
      extrapolateOpenBoundaryData2D< FluidPointer, OpenBoudaryPointer, SPHKernelFunction, EOS >(
            fluid, openBoundary, sphState );

   if constexpr( SPHState::SPHConfig::spaceDimension == 3 )
      extrapolateOpenBoundaryData3D< FluidPointer, OpenBoudaryPointer, SPHKernelFunction, EOS >(
            fluid, openBoundary, sphState );
}

template< typename ParticleSystem, typename SPHState >
template< typename FluidPointer,
          typename OpenBoudaryPointer,
          typename SPHKernelFunction,
          typename EOS >
void
WCSPH_DBC< ParticleSystem, SPHState >::extrapolateOpenBoundaryData2D( FluidPointer& fluid,
                                                                      OpenBoudaryPointer& openBoundary,
                                                                      SPHState& sphState )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   typename ParticleSystem::NeighborsLoopParams searchInFluid( fluid->particles );
   typename ParticleSystem::NeighborsLoopParams searchInOpenBoundary( openBoundary->particles );

   /* CONSTANT VARIABLES */
   const RealType searchRadius = fluid->particles->getSearchRadius();
   const RealType h = sphState.h;
   const RealType rho0 = sphState.rho0;
   const RealType m = sphState.mass;
   const VectorType gravity = sphState.gravity;
   const RealType extrapolationDetTreshold = openBoundary->parameters.extrapolationDetTreshold;

   typename EOS::ParamsType eosParams( sphState );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->particles->getPoints().getView();
   const auto view_rho = fluid->variables->rho.getView();
   auto view_Drho = fluid->variables->drho.getView();
   const auto view_v = fluid->variables->v.getView();
   auto view_a = fluid->variables->a.getView();

   const auto view_points_openBound = openBoundary->particles->getPoints().getView();
   auto view_rho_openBound = openBoundary->variables->rho.getView();
   auto view_v_openBound = openBoundary->variables->v.getView();

   //temp
   const VectorType bufferPosition = openBoundary->parameters.position;

   auto OpenBoundaryFluid = [=] __cuda_callable__ ( LocalIndexType i,
                                                    LocalIndexType j,
                                                    VectorType& r_i,
                                                    VectorType& v_i,
                                                    RealType& rho_i,
                                                    Matrix* A_gn,
                                                    VectorExtendedType* rho_gradrho_gn,
                                                    VectorExtendedType* vx_gradvx_gn,
                                                    VectorExtendedType* vy_gradvy_gn  ) mutable
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
         const RealType F = SPHKernelFunction::F( drs, h );
         const RealType W = SPHKernelFunction::W( drs, h );
         const VectorType gradW = r_ij * F;

         const RealType V = m / rho_j;

         *A_gn += matrixCorrection2D< Matrix >( W, gradW, r_ij, V );
         *rho_gradrho_gn += getVariableValueAndGradient2D< VectorExtendedType >( W, gradW, rho_j, V );
         *vx_gradvx_gn += getVariableValueAndGradient2D< VectorExtendedType >( W, gradW, v_j[ 0 ], V );
         *vy_gradvy_gn += getVariableValueAndGradient2D< VectorExtendedType >( W, gradW, v_j[ 1 ], V );
      }
   };

   auto particleLoopOpenBoundary = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points_openBound[ i ];
      const VectorType v_i = view_v_openBound[ i ];
      const RealType rho_i = view_rho_openBound[ i ];
      const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );
      const VectorType ghostNode_i = { bufferPosition[ 0 ] - ( r_i[ 0 ] - bufferPosition[ 0 ] ), r_i[ 1 ] }; //FIXME

      Matrix A_gn = 0.f;
      VectorExtendedType rho_gradrho_gn = 0.f;
      VectorExtendedType vx_gradvx_gn = 0.f;
      VectorExtendedType vy_gradvy_gn = 0.f;

      NeighborsLoop::exec(
            i, ghostNode_i, searchInFluid, OpenBoundaryFluid, v_i, rho_i, &A_gn, &rho_gradrho_gn, &vx_gradvx_gn, &vy_gradvy_gn );

      if( Matrices::determinant( A_gn ) > extrapolationDetTreshold )
      {
         VectorType r_ign = ghostNode_i - r_i;

         const VectorExtendedType rho_gradrho = Matrices::solve( A_gn, rho_gradrho_gn );
         const RealType rho_b = rho_gradrho[ 0 ] + rho_gradrho[ 1 ] * r_ign[ 0 ] + rho_gradrho[ 2 ] * r_ign[ 1 ];

         const VectorExtendedType vx_gradvx = Matrices::solve( A_gn, vx_gradvx_gn );
         const VectorExtendedType vy_gradvy = Matrices::solve( A_gn, vy_gradvy_gn );

         const RealType vx_b = vx_gradvx[ 0 ] + vx_gradvx[ 1 ] * r_ign[ 0 ] + vx_gradvx[ 2 ] * r_ign[ 1 ];
         const RealType vy_b = vy_gradvy[ 0 ] + vy_gradvy[ 1 ] * r_ign[ 0 ] + vy_gradvy[ 2 ] * r_ign[ 1 ];

         view_rho_openBound[ i ] = rho_b;
         view_v_openBound[ i ] = { vx_b, vy_b };
      }
      else if( A_gn( 0, 0 ) > 0.f )
      {
         RealType rho_b = rho_gradrho_gn[ 0 ] / A_gn( 0, 0 );
         RealType vx_b = vx_gradvx_gn[ 0 ] / A_gn( 0, 0 );
         RealType vy_b = vy_gradvy_gn[ 0 ] / A_gn( 0, 0 );

         view_rho_openBound[ i ] = rho_b;
         view_v_openBound[ i ] = { vx_b, vy_b };
      }
   };
   SPHParallelFor::exec(
         openBoundary->getFirstActiveParticle(), openBoundary->getLastActiveParticle() + 1, particleLoopOpenBoundary );
}

template< typename ParticleSystem, typename SPHState >
template< typename FluidPointer,
          typename OpenBoudaryPointer,
          typename SPHKernelFunction,
          typename EOS >
void
WCSPH_DBC< ParticleSystem, SPHState >::extrapolateOpenBoundaryData3D( FluidPointer& fluid,
                                                                      OpenBoudaryPointer& openBoundary,
                                                                      SPHState& sphState )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   typename ParticleSystem::NeighborsLoopParams searchInFluid( fluid->particles );
   typename ParticleSystem::NeighborsLoopParams searchInOpenBoundary( openBoundary->particles );

   /* CONSTANT VARIABLES */
   const RealType searchRadius = fluid->particles->getSearchRadius();
   const RealType h = sphState.h;
   const RealType rho0 = sphState.rho0;
   const RealType m = sphState.mass;
   const VectorType gravity = sphState.gravity;
   const RealType extrapolationDetTreshold = openBoundary->parameters.extrapolationDetTreshold;

   typename EOS::ParamsType eosParams( sphState );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->particles->getPoints().getView();
   const auto view_rho = fluid->variables->rho.getView();
   auto view_Drho = fluid->variables->drho.getView();
   const auto view_v = fluid->variables->v.getView();
   auto view_a = fluid->variables->a.getView();

   const auto view_points_openBound = openBoundary->particles->getPoints().getView();
   auto view_rho_openBound = openBoundary->variables->rho.getView();
   auto view_v_openBound = openBoundary->variables->v.getView();

   //temp
   const VectorType bufferPosition = openBoundary->parameters.position;

   auto OpenBoundaryFluid = [=] __cuda_callable__ ( LocalIndexType i,
                                                    LocalIndexType j,
                                                    VectorType& r_i,
                                                    VectorType& v_i,
                                                    RealType& rho_i,
                                                    Matrix* A_gn,
                                                    VectorExtendedType* rho_gradrho_gn,
                                                    VectorExtendedType* vx_gradvx_gn,
                                                    VectorExtendedType* vy_gradvy_gn,
                                                    VectorExtendedType* vz_gradvz_gn  ) mutable
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
         const RealType F = SPHKernelFunction::F( drs, h );
         const RealType W = SPHKernelFunction::W( drs, h );
         const VectorType gradW = r_ij * F;

         const RealType V = m / rho_j;

         *A_gn += matrixCorrection2D< Matrix >( W, gradW, r_ij, V );
         *rho_gradrho_gn += getVariableValueAndGradient2D< VectorExtendedType >( W, gradW, rho_j, V );
         *vx_gradvx_gn += getVariableValueAndGradient2D< VectorExtendedType >( W, gradW, v_j[ 0 ], V );
         *vy_gradvy_gn += getVariableValueAndGradient2D< VectorExtendedType >( W, gradW, v_j[ 1 ], V );
         *vz_gradvz_gn += getVariableValueAndGradient2D< VectorExtendedType >( W, gradW, v_j[ 2 ], V );
      }
   };

   auto particleLoopOpenBoundary = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points_openBound[ i ];
      const VectorType v_i = view_v_openBound[ i ];
      const RealType rho_i = view_rho_openBound[ i ];
      const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );
      const VectorType ghostNode_i = { bufferPosition[ 0 ] - ( r_i[ 0 ] - bufferPosition[ 0 ] ), r_i[ 1 ] }; //FIXME

      Matrix A_gn = 0.f;
      VectorExtendedType rho_gradrho_gn = 0.f;
      VectorExtendedType vx_gradvx_gn = 0.f;
      VectorExtendedType vy_gradvy_gn = 0.f;
      VectorExtendedType vz_gradvz_gn = 0.f;

      NeighborsLoop::exec( i,
                           ghostNode_i,
                           searchInFluid,
                           OpenBoundaryFluid,
                           v_i, rho_i, &A_gn, &rho_gradrho_gn, &vx_gradvx_gn, &vy_gradvy_gn, &vz_gradvz_gn );

      if( Matrices::determinant( A_gn ) > extrapolationDetTreshold )
      {
         VectorType r_ign = ghostNode_i - r_i;

         const VectorExtendedType rho_gradrho = Matrices::solve( A_gn, rho_gradrho_gn );
         const RealType rho_b = rho_gradrho[ 0 ] + rho_gradrho[ 1 ] * r_ign[ 0 ] + rho_gradrho[ 2 ] * r_ign[ 1 ] + rho_gradrho[ 3 ] * r_ign[ 2 ];

         const VectorExtendedType vx_gradvx = Matrices::solve( A_gn, vx_gradvx_gn );
         const VectorExtendedType vy_gradvy = Matrices::solve( A_gn, vy_gradvy_gn );
         const VectorExtendedType vz_gradvz = Matrices::solve( A_gn, vz_gradvz_gn );
         const RealType vx_b = vx_gradvx[ 0 ] + vx_gradvx[ 1 ] * r_ign[ 0 ] + vx_gradvx[ 2 ] * r_ign[ 1 ] + vx_gradvx[ 3 ] * r_ign[ 2 ];
         const RealType vy_b = vy_gradvy[ 0 ] + vy_gradvy[ 1 ] * r_ign[ 0 ] + vy_gradvy[ 2 ] * r_ign[ 1 ] + vy_gradvy[ 3 ] * r_ign[ 2 ];
         const RealType vz_b = vz_gradvz[ 0 ] + vz_gradvz[ 1 ] * r_ign[ 0 ] + vz_gradvz[ 2 ] * r_ign[ 1 ] + vz_gradvz[ 3 ] * r_ign[ 2 ];

         view_rho_openBound[ i ] = rho_b;
         view_v_openBound[ i ] = { vx_b, vy_b, vz_b };
      }
      else if( A_gn( 0, 0 ) > 0.f )
      {
         RealType rho_b = rho_gradrho_gn[ 0 ] / A_gn( 0, 0 );
         RealType vx_b = vx_gradvx_gn[ 0 ] / A_gn( 0, 0 );
         RealType vy_b = vy_gradvy_gn[ 0 ] / A_gn( 0, 0 );
         RealType vz_b = vz_gradvz_gn[ 0 ] / A_gn( 0, 0 );

         view_rho_openBound[ i ] = rho_b;
         view_v_openBound[ i ] = { vx_b, vy_b, vz_b };
      }
   };
   SPHParallelFor::exec(
         openBoundary->getFirstActiveParticle(), openBoundary->getLastActiveParticle() + 1, particleLoopOpenBoundary );
}

} // SPH
} // ParticleSystem
} // TNL

