#include "Interactions.h"
#include "../../customParallelFor.h"


namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename ParticleSystem, typename SPHFluidConfig, typename Variables >
template< typename FluidPointer,
          typename OpenBoudaryPointer,
          typename SPHKernelFunction,
          typename DiffusiveTerm,
          typename ViscousTerm,
          typename EOS,
          typename SPHState >
void
WCSPH_DBC< ParticleSystem, SPHFluidConfig, Variables >::interactionWithOpenBoundary( FluidPointer& fluid,
                                                                                     OpenBoudaryPointer& openBoundary,
                                                                                     SPHState& sphState )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   typename ParticleSystem::NeighborsLoopParams searchInFluid( fluid->particles );
   typename ParticleSystem::NeighborsLoopParams searchInOpenBoundary( openBoundary->particles );

   /* CONSTANT VARIABLES */
   const RealType searchRadius = fluid->particles->getSearchRadius();
   const RealType h = sphState.h;
   const RealType m = sphState.mass;
   const VectorType gravity = sphState.gravity;

   typename DiffusiveTerm::ParamsType diffusiveTermsParams( sphState );
   typename ViscousTerm::ParamsType viscousTermTermsParams( sphState );
   typename EOS::ParamsType eosParams( sphState );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->particles->getPoints().getView();
   const auto view_rho = fluid->variables->rho.getView();
   auto view_Drho = fluid->variables->drho.getView();
   const auto view_v = fluid->variables->v.getView();
   auto view_a = fluid->variables->a.getView();

   const auto view_points_openBound = openBoundary->particles->getPoints().getView();
   const auto view_rho_openBound = openBoundary->variables->rho.getView();
   auto view_Drho_openBound = openBoundary->variables->drho.getView();
   const auto view_v_openBound = openBoundary->variables->v.getView();

   auto FluidOpenBoundary = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j,
         VectorType& r_i, VectorType& v_i, RealType& rho_i, RealType& p_i, RealType* drho_i, VectorType* a_i ) mutable
   {
      const VectorType r_j = view_points_openBound[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if (drs <= searchRadius )
      {
         const VectorType v_j = view_v_openBound[ j ];
         const RealType rho_j = view_rho_openBound[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );

         /* Interaction: */
         const VectorType v_ij = v_i - v_j;

         const RealType F = SPHKernelFunction::F( drs, h );
         const VectorType gradW = r_ij * F;

         const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, drs, diffusiveTermsParams );
         const RealType diffTerm =  psi * ( r_ij, gradW ) * m / rho_j;
         *drho_i += ( v_ij, gradW ) * m - diffTerm;

         const RealType p_term = ( p_i + p_j ) / ( rho_i * rho_j );
         const RealType visco =  ViscousTerm::Pi( rho_i, rho_j, drs, ( r_ij, v_ij ), viscousTermTermsParams );
         *a_i += ( -1.0f ) * ( p_term + visco ) * gradW * m;
      }
   };

   auto OpenBoundaryFluid = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j,
         VectorType& r_i, VectorType& v_i, RealType& rho_i, RealType& p_i, RealType* drho_i, VectorType* a_i ) mutable
   {
      const VectorType r_j = view_points_openBound[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if (drs <= searchRadius )
      {
         const VectorType v_j = view_v_openBound[ j ];
         const RealType rho_j = view_rho_openBound[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );

         /* Interaction: */
         const VectorType v_ij = v_i - v_j;

         const RealType F = SPHKernelFunction::F( drs, h );
         const VectorType gradW = r_ij * F;

         const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, drs, diffusiveTermsParams );
         const RealType diffTerm =  psi * ( r_ij, gradW ) * m / rho_j;
         *drho_i += ( v_ij, gradW ) * m - diffTerm;

         const RealType p_term = ( p_i + p_j ) / ( rho_i * rho_j );
         const RealType visco =  ViscousTerm::Pi( rho_i, rho_j, drs, ( r_ij, v_ij ), viscousTermTermsParams );
         *a_i += ( -1.0f ) * ( p_term + visco ) * gradW * m;
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

      NeighborsLoop::exec( i, r_i, searchInOpenBoundary, FluidOpenBoundary, v_i, rho_i, p_i, &drho_i, &a_i );

      view_Drho[ i ] += drho_i;
      view_a[ i ] += a_i;
   };
   SPHParallelFor::exec( fluid->getFirstActiveParticle(), fluid->getLastActiveParticle() + 1, particleLoop );

   //auto particleLoopBoundary = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   //{
   //   const VectorType r_i = view_points_openBound[ i ];
   //   const VectorType v_i = view_v_openBound[ i ];
   //   const RealType rho_i = view_rho_openBound[ i ];
   //   const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );

   //   RealType drho_i = 0.f;

   //   NeighborsLoop::exec( i, r_i, searchInFluid, BoundFluid, v_i, rho_i, p_i, &drho_i );

   //   view_Drho_openBound[ i ] = drho_i;
   //};
}

template< typename ParticleSystem, typename SPHFluidConfig, typename Variables >
template< typename FluidPointer,
          typename BoundaryPointer,
          typename OpenBoudaryPointer,
          typename SPHKernelFunction,
          typename DiffusiveTerm,
          typename ViscousTerm,
          typename EOS,
          typename SPHState >
void
WCSPH_DBC< ParticleSystem, SPHFluidConfig, Variables >::interactionWithOpenBoundary( FluidPointer& fluid,
                                                                                     BoundaryPointer& boundary,
                                                                                     OpenBoudaryPointer& openBoundary,
                                                                                     SPHState& sphState )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   typename ParticleSystem::NeighborsLoopParams searchInFluid( fluid->particles );
   //typename ParticleSystem::NeighborsLoopParams searchInBound( boundary->particles );
   typename ParticleSystem::NeighborsLoopParams searchInOpenBoundary( openBoundary->particles );

   /* CONSTANT VARIABLES */
   const RealType searchRadius = fluid->particles->getSearchRadius();
   const RealType h = sphState.h;
   const RealType m = sphState.mass;
   const VectorType gravity = sphState.gravity;

   typename DiffusiveTerm::ParamsType diffusiveTermsParams( sphState );
   typename ViscousTerm::ParamsType viscousTermTermsParams( sphState );
   typename EOS::ParamsType eosParams( sphState );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->particles->getPoints().getView();
   const auto view_rho = fluid->variables->rho.getView();
   auto view_Drho = fluid->variables->drho.getView();
   const auto view_v = fluid->variables->v.getView();
   auto view_a = fluid->variables->a.getView();

   const auto view_points_bound = boundary->particles->getPoints().getView();
   const auto view_rho_bound = boundary->variables->rho.getView();
   auto view_Drho_bound = boundary->variables->drho.getView();
   const auto view_v_bound = boundary->variables->v.getView();

   const auto view_points_openBound = openBoundary->particles->getPoints().getView();
   auto view_rho_openBound = openBoundary->variables->rho.getView();
   auto view_v_openBound = openBoundary->variables->v.getView();

   //temp
   const VectorType bufferPosition = openBoundary->parameters.position;

   auto FluidOpenBoundary = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j,
         VectorType& r_i, VectorType& v_i, RealType& rho_i, RealType& p_i, RealType* drho_i, VectorType* a_i ) mutable
   {
      const VectorType r_j = view_points_openBound[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if (drs <= searchRadius )
      {
         const VectorType v_j = view_v_openBound[ j ];
         const RealType rho_j = view_rho_openBound[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );

         /* Interaction: */
         const VectorType v_ij = v_i - v_j;

         const RealType F = SPHKernelFunction::F( drs, h );
         const VectorType gradW = r_ij * F;

         const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, drs, diffusiveTermsParams );
         const RealType diffTerm =  psi * ( r_ij, gradW ) * m / rho_j;
         *drho_i += ( v_ij, gradW ) * m - diffTerm;

         const RealType p_term = ( p_i + p_j ) / ( rho_i * rho_j );
         const RealType visco =  ViscousTerm::Pi( rho_i, rho_j, drs, ( r_ij, v_ij ), viscousTermTermsParams );
         *a_i += ( -1.0f ) * ( p_term + visco ) * gradW * m;
      }
   };

   auto OpenBoundaryFluid = [=] __cuda_callable__ ( LocalIndexType i,
                                                    LocalIndexType j,
                                                    VectorType& r_i,
                                                    VectorType& v_i,
                                                    RealType& rho_i,
                                                    Matrix* A_gn,
                                                    VectorExtendedType* rhogradrho_gn,
                                                    VectorExtendedType* vxgradvx_gn,
                                                    VectorExtendedType* vygradvy_gn  ) mutable
   {
      const VectorType r_j = view_points_openBound[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if (drs <= searchRadius )
      {
         const VectorType v_j = view_v_openBound[ j ];
         const RealType rho_j = view_rho_openBound[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );

         /* Interaction: */
         const RealType F = SPHKernelFunction::F( drs, h );
         const RealType W = SPHKernelFunction::W( drs, h );
         const VectorType gradW = r_ij * F;

         const RealType V = m / rho_j;

         Matrix A_local = 0.f;
         VectorExtendedType rhogradrho_local = 0.f;
         VectorExtendedType vxgradvx_gn_local = 0.f;
         VectorExtendedType vygradvy_gn_local = 0.f;

         if( SPHConfig::spaceDimension == 2 )
         {
            VectorExtendedType vxgradvx_local = 0.f;
            VectorExtendedType vygradvy_local = 0.f;

            A_local( 0, 0 ) = W * V;          A_local( 0, 1 ) = r_ij[ 0 ] * W * V;          A_local( 0, 2 ) = r_ij[ 1 ] * W * V;
            A_local( 1, 0 ) = gradW[ 0 ] * V; A_local( 1, 1 ) = r_ij[ 0 ] * gradW[ 0 ] * V; A_local( 1, 2 ) = r_ij[ 1 ] * gradW[ 0 ] * V;
            A_local( 2, 0 ) = gradW[ 1 ] * V; A_local( 2, 1 ) = r_ij[ 0 ] * gradW[ 1 ] * V; A_local( 2, 2 ) = r_ij[ 1 ] * gradW[ 1 ] * V;

            rhogradrho_local[ 0 ] = W * m; rhogradrho_local[ 1 ] = gradW[ 0 ] * m; rhogradrho_local[ 2 ] = gradW[ 1 ] * m;

            vxgradvx_local[ 0 ] = v_j[ 0 ] * V; vxgradvx_local[ 1 ] = gradW[ 0 ] * v_j[ 0 ] * V; vxgradvx_local[ 2 ] = gradW[ 1 ] * v_j[ 0 ] * V;
            vygradvy_local[ 0 ] = v_j[ 1 ] * V; vygradvy_local[ 1 ] = gradW[ 0 ] * v_j[ 1 ] * V; vygradvy_local[ 2 ] = gradW[ 1 ] * v_j[ 1 ] * V;

            *A_gn += A_local;
            *rhogradrho_gn += rhogradrho_local;
            *vxgradvx_gn += vxgradvx_local;
            *vygradvy_gn += vygradvy_local;
         }
      }
   };

   //EXPERIMETN:
   auto BoundOpenBoundary = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j,
         VectorType& r_i, VectorType& v_i, RealType& rho_i, RealType& p_i, RealType* drho_i ) mutable
   {
      const VectorType r_j = view_points_openBound[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         const VectorType v_j = view_v_openBound[ j ];
         const RealType rho_j = view_rho_openBound[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );

         /* Interaction */
         const VectorType v_ij = v_i - v_j;

         const RealType F = SPHKernelFunction::F( drs, h );
         const VectorType gradW = r_ij * F;

         const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, drs, diffusiveTermsParams );
         const RealType diffTerm =  psi * ( r_ij, gradW ) * m / rho_j;
         *drho_i += ( v_ij, gradW ) * m - diffTerm;
         //printf(" drho in BoundOpenBoundary: %f \n", *drho_i );
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

      NeighborsLoop::exec( i, r_i, searchInOpenBoundary, FluidOpenBoundary, v_i, rho_i, p_i, &drho_i, &a_i );

      view_Drho[ i ] += drho_i;
      view_a[ i ] += a_i;
   };
   SPHParallelFor::exec( fluid->getFirstActiveParticle(), fluid->getLastActiveParticle() + 1, particleLoop );

   //EXPERIMENT;
   auto particleLoopBoundary = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points_bound[ i ];
      const VectorType v_i = view_v_bound[ i ];
      const RealType rho_i = view_rho_bound[ i ];
      const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );

      RealType drho_i = 0.f;

      NeighborsLoop::exec( i, r_i, searchInOpenBoundary, BoundOpenBoundary, v_i, rho_i, p_i, &drho_i );

      view_Drho_bound[ i ] += drho_i;
   };
   SPHParallelFor::exec( boundary->getFirstActiveParticle(), boundary->getLastActiveParticle() + 1, particleLoopBoundary );

   auto particleLoopOpenBoundary = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points_openBound[ i ];
      const VectorType v_i = view_v_openBound[ i ];
      const RealType rho_i = view_rho_openBound[ i ];
      const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );
      const VectorType ghostNode_i = { bufferPosition[ 0 ] - ( r_i[ 0 ] - bufferPosition[ 0 ] ) * 2.f , r_i[ 1 ] }; //FIXME

      Matrix A_gn = 0.f;
      VectorExtendedType rhogradrho_gn = 0.f;
      VectorExtendedType vxgradvx_gn = 0.f;
      VectorExtendedType vygradvy_gn = 0.f;

      NeighborsLoop::exec(
            i, r_i, searchInFluid, OpenBoundaryFluid, v_i, rho_i, &A_gn, &rhogradrho_gn, &vxgradvx_gn, &vygradvy_gn );

      if( Matrices::determinant( A_gn ) > 0.001 )
      {
         VectorType r_ign = ghostNode_i - r_i;

         const VectorExtendedType rhogradrho = Matrices::solve( A_gn, rhogradrho_gn );
         const RealType rho_bound = rhogradrho[ 0 ] + rhogradrho[ 1 ] * r_ign[ 0 ] + rhogradrho[ 2 ] * r_ign[ 1 ];

         const VectorExtendedType vgradvx = Matrices::solve( A_gn, vxgradvx_gn );
         const VectorExtendedType vgradvy = Matrices::solve( A_gn, vygradvy_gn );

         const RealType vx_bound = vgradvx[ 0 ] + vgradvx[ 1 ] * r_ign[ 0 ] + vgradvx[ 2 ] * r_ign[ 1 ];
         const RealType vy_bound = vgradvy[ 0 ] + vgradvy[ 1 ] * r_ign[ 0 ] + vgradvy[ 2 ] * r_ign[ 1 ];

         view_rho_openBound[ i ] = rho_bound;
         const VectorType v_b = { vx_bound, vy_bound };
         view_v_openBound[ i ] = v_b;
      }
      else if( A_gn( 0, 0 ) > 0.f )
      {
         RealType rho_bound = rhogradrho_gn[ 0 ] / A_gn( 0, 0 );
         RealType vx_bound = vxgradvx_gn[ 0 ] / A_gn( 0, 0 );
         RealType vy_bound = vygradvy_gn[ 0 ] / A_gn( 0, 0 );

         view_rho_openBound[ i ] = rho_bound;
         const VectorType v_b = { vx_bound, vy_bound };
         view_v_openBound[ i ] = v_b;
      }
      else
      {

      }


   };

   //TODO: Temp to test.
   if( openBoundary->parameters.identifier == "outlet" )
   {
      std::cout << "Running outlet interpolation." << std::endl;
      SPHParallelFor::exec(
            openBoundary->getFirstActiveParticle(), openBoundary->getLastActiveParticle() + 1, particleLoopOpenBoundary );
      std::cout << "Running outlet interpolation. DONE." << std::endl;

   }

}

template< typename ParticleSystem, typename SPHFluidConfig, typename Variables >
template< typename FluidPointer,
          typename OpenBoudaryPointer,
          typename SPHKernelFunction,
          typename EOS,
          typename SPHState >
void
WCSPH_DBC< ParticleSystem, SPHFluidConfig, Variables >::extrapolateOpenBoundaryData( FluidPointer& fluid,
                                                                                     OpenBoudaryPointer& openBoundary,
                                                                                     SPHState& sphState )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   typename ParticleSystem::NeighborsLoopParams searchInFluid( fluid->particles );
   typename ParticleSystem::NeighborsLoopParams searchInOpenBoundary( openBoundary->particles );

   /* CONSTANT VARIABLES */
   const RealType searchRadius = fluid->particles->getSearchRadius();
   const RealType h = sphState.h;
   const RealType m = sphState.mass;
   const VectorType gravity = sphState.gravity;

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
                                                    VectorExtendedType* rhogradrho_gn,
                                                    VectorExtendedType* vxgradvx_gn,
                                                    VectorExtendedType* vygradvy_gn  ) mutable
   {
      const VectorType r_j = view_points_openBound[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if (drs <= searchRadius )
      {
         const VectorType v_j = view_v_openBound[ j ];
         const RealType rho_j = view_rho_openBound[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );

         /* Interaction: */
         const RealType F = SPHKernelFunction::F( drs, h );
         const RealType W = SPHKernelFunction::W( drs, h );
         const VectorType gradW = r_ij * F;

         const RealType V = m / rho_j;

         Matrix A_local = 0.f;
         VectorExtendedType rhogradrho_local = 0.f;
         VectorExtendedType vxgradvx_gn_local = 0.f;
         VectorExtendedType vygradvy_gn_local = 0.f;

         if( SPHConfig::spaceDimension == 2 )
         {
            VectorExtendedType vxgradvx_local = 0.f;
            VectorExtendedType vygradvy_local = 0.f;

            A_local( 0, 0 ) = W * V;          A_local( 0, 1 ) = r_ij[ 0 ] * W * V;          A_local( 0, 2 ) = r_ij[ 1 ] * W * V;
            A_local( 1, 0 ) = gradW[ 0 ] * V; A_local( 1, 1 ) = r_ij[ 0 ] * gradW[ 0 ] * V; A_local( 1, 2 ) = r_ij[ 1 ] * gradW[ 0 ] * V;
            A_local( 2, 0 ) = gradW[ 1 ] * V; A_local( 2, 1 ) = r_ij[ 0 ] * gradW[ 1 ] * V; A_local( 2, 2 ) = r_ij[ 1 ] * gradW[ 1 ] * V;

            rhogradrho_local[ 0 ] = W * m; rhogradrho_local[ 1 ] = gradW[ 0 ] * m; rhogradrho_local[ 2 ] = gradW[ 1 ] * m;

            vxgradvx_local[ 0 ] = v_j[ 0 ] * V; vxgradvx_local[ 1 ] = gradW[ 0 ] * v_j[ 0 ] * V; vxgradvx_local[ 2 ] = gradW[ 1 ] * v_j[ 0 ] * V;
            vygradvy_local[ 0 ] = v_j[ 1 ] * V; vygradvy_local[ 1 ] = gradW[ 0 ] * v_j[ 1 ] * V; vygradvy_local[ 2 ] = gradW[ 1 ] * v_j[ 1 ] * V;

            *A_gn += A_local;
            *rhogradrho_gn += rhogradrho_local;
            *vxgradvx_gn += vxgradvx_local;
            *vygradvy_gn += vygradvy_local;
         }
      }
   };

   auto particleLoopOpenBoundary = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points_openBound[ i ];
      const VectorType v_i = view_v_openBound[ i ];
      const RealType rho_i = view_rho_openBound[ i ];
      const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );
      const VectorType ghostNode_i = { bufferPosition[ 0 ] - ( r_i[ 0 ] - bufferPosition[ 0 ] ) * 2.f , r_i[ 1 ] }; //FIXME

      Matrix A_gn = 0.f;
      VectorExtendedType rhogradrho_gn = 0.f;
      VectorExtendedType vxgradvx_gn = 0.f;
      VectorExtendedType vygradvy_gn = 0.f;

      NeighborsLoop::exec(
            i, r_i, searchInFluid, OpenBoundaryFluid, v_i, rho_i, &A_gn, &rhogradrho_gn, &vxgradvx_gn, &vygradvy_gn );

      if( Matrices::determinant( A_gn ) > 0.001 )
      {
         VectorType r_ign = ghostNode_i - r_i;

         const VectorExtendedType rhogradrho = Matrices::solve( A_gn, rhogradrho_gn );
         const RealType rho_bound = rhogradrho[ 0 ] + rhogradrho[ 1 ] * r_ign[ 0 ] + rhogradrho[ 2 ] * r_ign[ 1 ];

         const VectorExtendedType vgradvx = Matrices::solve( A_gn, vxgradvx_gn );
         const VectorExtendedType vgradvy = Matrices::solve( A_gn, vygradvy_gn );

         const RealType vx_bound = vgradvx[ 0 ] + vgradvx[ 1 ] * r_ign[ 0 ] + vgradvx[ 2 ] * r_ign[ 1 ];
         const RealType vy_bound = vgradvy[ 0 ] + vgradvy[ 1 ] * r_ign[ 0 ] + vgradvy[ 2 ] * r_ign[ 1 ];

         view_rho_openBound[ i ] = rho_bound;
         const VectorType v_b = { vx_bound, vy_bound };
         view_v_openBound[ i ] = v_b;
      }
      else if( A_gn( 0, 0 ) > 0.f )
      {
         RealType rho_bound = rhogradrho_gn[ 0 ] / A_gn( 0, 0 );
         RealType vx_bound = vxgradvx_gn[ 0 ] / A_gn( 0, 0 );
         RealType vy_bound = vygradvy_gn[ 0 ] / A_gn( 0, 0 );

         view_rho_openBound[ i ] = rho_bound;
         const VectorType v_b = { vx_bound, vy_bound };
         view_v_openBound[ i ] = v_b;
      }
      else
      {

      }

   };

   //TODO: Temp to test.
   if( openBoundary->parameters.identifier == "outlet" )
   {
      std::cout << "Running outlet interpolation." << std::endl;
      SPHParallelFor::exec(
            openBoundary->getFirstActiveParticle(), openBoundary->getLastActiveParticle() + 1, particleLoopOpenBoundary );
      std::cout << "Running outlet interpolation. DONE." << std::endl;

   }

}


} // SPH
} // ParticleSystem
} // TNL

