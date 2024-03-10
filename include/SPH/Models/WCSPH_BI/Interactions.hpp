#include "Interactions.h"

namespace TNL {
namespace SPH {

template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename BoudaryPointer >
void
WCSPH_BI< Particles, ModelConfig >::interaction( FluidPointer& fluid,
                                                 BoudaryPointer& boundary,
                                                 ModelParams& modelParams )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   GlobalIndexType numberOfParticles = fluid->particles->getNumberOfParticles();
   GlobalIndexType numberOfParticles_bound = boundary->particles->getNumberOfParticles();
   const RealType searchRadius = fluid->particles->getSearchRadius();

   typename Particles::NeighborsLoopParams searchInFluid( fluid->particles );
   typename Particles::NeighborsLoopParams searchInBound( boundary->particles );

   /* CONSTANT VARIABLES */
   const RealType h = modelParams.h;
   const RealType m = modelParams.mass;
   const RealType ds = modelParams.boundaryElementSize;
   const RealType rho0 = modelParams.rho0;

   typename DiffusiveTerm::ParamsType diffusiveTermsParams( modelParams );
   typename ViscousTerm::ParamsType viscousTermsParams( modelParams );
   typename EOS::ParamsType eosParams( modelParams );
   typename BoundaryViscousTerm::ParamsType boundaryViscoParams( modelParams );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->particles->getPoints().getView();
   const auto view_rho = fluid->variables->rho.getView();
   auto view_Drho = fluid->variables->drho.getView();
   const auto view_v = fluid->variables->v.getView();
   auto view_a = fluid->variables->a.getView();
   auto view_gamma = fluid->variables->gamma.getView();

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

         const RealType F = KernelFunction::F( drs, h );
         const RealType W = KernelFunction::W( drs, h );
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

         const RealType W = KernelFunction::W( drs, h );

         *drho_i += ( -1.f ) * ( v_ij, n_j ) * W * rho_j * ds;

         const RealType p_term = ( p_i + p_j ) / ( rho_i * rho_j );
         const RealType visco =  ViscousTerm::Pi( rho_i, rho_j, drs, ( r_ij, v_ij ), viscousTermsParams );
         *a_i += ( p_term + visco ) * n_j * W * rho_j * ds + BoundaryViscousTerm::Xi( r_ij, v_ij, n_j, boundaryViscoParams );
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

      TNL::ParticleSystem::NeighborsLoop::exec( i, r_i, searchInFluid, FluidFluid, v_i, rho_i, p_i, &drho_i, &a_i, &gamma_i );
      TNL::ParticleSystem::NeighborsLoopAnotherSet::exec( i, r_i, searchInBound, FluidBound, v_i, rho_i, p_i, &drho_i, &a_i );

      //if( gamma_i > 0.01f ) {
      //   view_Drho[ i ] = drho_i / gamma_i;
      //   view_a[ i ] = a_i / gamma_i + gravity;
      //}
      //else {
      //   view_Drho[ i ] = 0.f;
      //   view_a[ i ] = 0.f + gravity;
      //}

      view_Drho[ i ] = drho_i ;
      view_a[ i ] = a_i;
      view_gamma[ i ] = gamma_i;

   };
   TNL::Algorithms::parallelFor< DeviceType >( fluid->getFirstActiveParticle(), fluid->getLastActiveParticle() + 1, particleLoop );

}

template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename BoudaryPointer >
void
WCSPH_BI< Particles, ModelConfig >::updateSolidBoundary( FluidPointer& fluid,
                                                         BoudaryPointer& boundary,
                                                         ModelParams& modelParams )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   typename Particles::NeighborsLoopParams searchInFluid( fluid->particles );

   /* CONSTANT VARIABLES */
   const RealType searchRadius = fluid->particles->getSearchRadius();
   const RealType h = modelParams.h;
   const RealType m = modelParams.mass;
   const RealType rho0 = modelParams.rho0;

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->particles->getPoints().getView();
   const auto view_rho = fluid->variables->rho.getView();

   const auto view_points_bound = boundary->particles->getPoints().getView();
   auto view_rho_bound = boundary->variables->rho.getView();
   auto view_gamma_bound = boundary->variables->gamma.getView();

   auto BoundFluid = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j,
         VectorType& r_i, RealType* rho_i, RealType* gamma_i ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         const RealType rho_j = view_rho[ j ];

         const RealType W = KernelFunction::W( drs, h );

         *rho_i += W * m;
         *gamma_i += W * m / rho_j;
      }
   };

   auto particleLoopBoundary = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points_bound[ i ];

      RealType rho_i = 0.f;
      RealType gamma_i = 0.f;

      TNL::ParticleSystem::NeighborsLoopAnotherSet::exec( i, r_i, searchInFluid, BoundFluid, &rho_i, &gamma_i );

      view_rho_bound[ i ] = rho_i;
      view_gamma_bound[ i ] = gamma_i;
   };
   TNL::Algorithms::parallelFor< DeviceType >(
         boundary->getFirstActiveParticle(), boundary->getLastActiveParticle() + 1, particleLoopBoundary );
}

template< typename Particles, typename ModelConfig >
template< typename OpenBoundaryPointer, typename BoudaryPointer >
void
WCSPH_BI< Particles, ModelConfig >::updateSolidBoundaryOpenBoundary( BoudaryPointer& boundary,
                                                                     OpenBoundaryPointer& openBoundary,
                                                                     ModelParams& modelParams )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   typename Particles::NeighborsLoopParams searchInOpenBoundary( openBoundary->particles );

   /* CONSTANT VARIABLES */
   const RealType searchRadius = openBoundary->particles->getSearchRadius();
   const RealType h = modelParams.h;
   const RealType m = modelParams.mass;
   const RealType rho0 = modelParams.rho0;

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points_openBound = openBoundary->particles->getPoints().getView();
   const auto view_rho_openBound = openBoundary->variables->rho.getView();

   const auto view_points_bound = boundary->particles->getPoints().getView();
   auto view_rho_bound = boundary->variables->rho.getView();
   auto view_gamma_bound = boundary->variables->gamma.getView();

   auto BoundOpenBoundary = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j,
         VectorType& r_i, RealType* rho_i, RealType* gamma_i ) mutable
   {
      const VectorType r_j = view_points_openBound[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         const RealType rho_j = view_rho_openBound[ j ];

         const RealType W = KernelFunction::W( drs, h );

         *rho_i += W * m;
         *gamma_i += W * m / rho_j;
      }
   };

   auto particleLoopBoundary = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points_bound[ i ];

      RealType rho_i = 0.f;
      RealType gamma_i = 0.f;

      TNL::ParticleSystem::NeighborsLoopAnotherSet::exec( i, r_i, searchInOpenBoundary, BoundOpenBoundary, &rho_i, &gamma_i );

      view_rho_bound[ i ] += rho0;
      view_gamma_bound[ i ] += gamma_i;
   };
   TNL::Algorithms::parallelFor< DeviceType >(
         boundary->getFirstActiveParticle(), boundary->getLastActiveParticle() + 1, particleLoopBoundary );
}

template< typename Particles, typename ModelConfig >
template< typename EquationOfState, typename PhysicalObjectPointer >
void
WCSPH_BI< Particles, ModelConfig >::computePressureFromDensity( PhysicalObjectPointer& physicalObject, ModelParams& modelParams )
{
   auto view_rho = physicalObject->getVariables()->rho.getView();
   auto view_p = physicalObject->getVariables()->p.getView();

   typename EOS::ParamsType eosParams( modelParams );

   auto init = [=] __cuda_callable__ ( int i ) mutable
   {
      view_p[ i ] = EquationOfState::DensityToPressure( view_rho[ i ], eosParams );
   };
   Algorithms::parallelFor< DeviceType >( physicalObject->getFirstActiveParticle(), physicalObject->getLastActiveParticle() + 1, init );
}


template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename OpenBoudaryPointer >
void
WCSPH_BI< Particles, ModelConfig >::interactionWithOpenBoundary( FluidPointer& fluid,
                                                                 OpenBoudaryPointer& openBoundary,
                                                                 ModelParams& modelParams )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   typename Particles::NeighborsLoopParams searchInFluid( fluid->particles );
   typename Particles::NeighborsLoopParams searchInOpenBoundary( openBoundary->particles );

   /* CONSTANT VARIABLES */
   const RealType searchRadius = fluid->particles->getSearchRadius();
   const RealType h = modelParams.h;
   const RealType m = modelParams.mass;

   typename DiffusiveTerm::ParamsType diffusiveTermsParams( modelParams );
   typename ViscousTerm::ParamsType viscousTermsParams( modelParams );
   typename EOS::ParamsType eosParams( modelParams );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->particles->getPoints().getView();
   const auto view_rho = fluid->variables->rho.getView();
   auto view_Drho = fluid->variables->drho.getView();
   const auto view_v = fluid->variables->v.getView();
   auto view_a = fluid->variables->a.getView();
   auto view_gamma = fluid->variables->gamma.getView();

   const auto view_points_openBound = openBoundary->particles->getPoints().getView();
   auto view_rho_openBound = openBoundary->variables->rho.getView();
   auto view_v_openBound = openBoundary->variables->v.getView();

   const auto zoneParticleIndices_view = openBoundary->zone.getParticlesInZone().getConstView();
   const GlobalIndexType numberOfZoneParticles = openBoundary->zone.getNumberOfParticles();

   auto FluidOpenBoundary = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j,
         VectorType& r_i, VectorType& v_i, RealType& rho_i, RealType& p_i, RealType* drho_i, VectorType* a_i, RealType* gamma_i ) mutable
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

         const RealType F = KernelFunction::F( drs, h );
         const RealType W = KernelFunction::W( drs, h );
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

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const GlobalIndexType p = zoneParticleIndices_view[ i ];
      const VectorType r_i = view_points[ p ];
      const VectorType v_i = view_v[ p ];
      const RealType rho_i = view_rho[ p ];
      const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );

      VectorType a_i = 0.f;
      RealType drho_i = 0.f;
      RealType gamma_i = 0.f;

      TNL::ParticleSystem::NeighborsLoopAnotherSet::exec( p, r_i, searchInOpenBoundary, FluidOpenBoundary, v_i, rho_i, p_i, &drho_i, &a_i, &gamma_i );

      view_Drho[ p ] += drho_i;
      view_a[ p ] += a_i;
      view_gamma[ p ] += gamma_i;

   };
   Algorithms::parallelFor< DeviceType >( 0, numberOfZoneParticles, particleLoop );
}

template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename OpenBoudaryPointer >
void
WCSPH_BI< Particles, ModelConfig >::interactionWithBoundaryPatches( FluidPointer& fluid,
                                                                    OpenBoudaryPointer& openBoundary,
                                                                    ModelParams& modelParams )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   typename Particles::NeighborsLoopParams searchInOpenBoundary( openBoundary->particles );

   /* CONSTANT VARIABLES */
   const RealType searchRadius = fluid->particles->getSearchRadius();
   const RealType h = modelParams.h;
   const RealType ds = modelParams.boundaryElementSize;
   const RealType m = modelParams.mass;
   typename DiffusiveTerm::ParamsType diffusiveTermsParams( modelParams );
   typename ViscousTerm::ParamsType viscousTermsParams( modelParams );
   typename EOS::ParamsType eosParams( modelParams );
   typename BoundaryViscousTerm::ParamsType boundaryViscoParams( modelParams );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->particles->getPoints().getView();
   const auto view_rho = fluid->variables->rho.getView();
   auto view_Drho = fluid->variables->drho.getView();
   const auto view_v = fluid->variables->v.getView();
   auto view_a = fluid->variables->a.getView();
   auto view_gamma = fluid->variables->gamma.getView();
   const auto view_points_openBound = openBoundary->particles->getPoints().getView();
   auto view_rho_openBound = openBoundary->variables->rho.getView();
   auto view_v_openBound = openBoundary->variables->v.getView();
   const auto view_n_bound = openBoundary->variables->n.getView();

   const auto zoneParticleIndices_view = openBoundary->zone.getParticlesInZone().getConstView();
   const GlobalIndexType numberOfZoneParticles = openBoundary->zone.getNumberOfParticles();

   auto FluidBound = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j,
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
         const VectorType n_j = view_n_bound[ j ];

         /* Interaction: */
         const VectorType v_ij = v_i - v_j;

         const RealType W = KernelFunction::W( drs, h );

         *drho_i += ( -1.f ) * ( v_ij, n_j ) * W * rho_j * ds;

         const RealType p_term = ( p_i + p_j ) / ( rho_i * rho_j );
         const RealType visco =  ViscousTerm::Pi( rho_i, rho_j, drs, ( r_ij, v_ij ), viscousTermsParams );
         *a_i += ( p_term + visco ) * n_j * W * rho_j * ds + BoundaryViscousTerm::Xi( r_ij, v_ij, n_j, boundaryViscoParams );
      }
   };

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const GlobalIndexType p = zoneParticleIndices_view[ i ];
      const VectorType r_i = view_points[ p ];
      const VectorType v_i = view_v[ p ];
      const RealType rho_i = view_rho[ p ];
      const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );

      VectorType a_i = 0.f;
      RealType drho_i = 0.f;

      TNL::ParticleSystem::NeighborsLoopAnotherSet::exec(
            p, r_i, searchInOpenBoundary, FluidBound, v_i, rho_i, p_i, &drho_i, &a_i );

      view_Drho[ p ] += drho_i;
      view_a[ p ] += a_i;

   };
   Algorithms::parallelFor< DeviceType >( 0, numberOfZoneParticles, particleLoop );
}

template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename BoundaryPointer>
void
WCSPH_BI< Particles, ModelConfig >::finalizeInteraction( FluidPointer& fluid,
                                                         BoundaryPointer& boundary,
                                                         ModelParams& modelParams )
{
   //finalize fluid-boundary interactions
   const VectorType gravity = modelParams.gravity;

   auto view_Drho = fluid->variables->drho.getView();
   auto view_a = fluid->variables->a.getView();
   auto view_gamma = fluid->variables->gamma.getView();

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const RealType gamma_i = view_gamma[ i ];

      if( gamma_i > 0.01f ) {
         view_Drho[ i ] = view_Drho[ i ] / gamma_i;
         view_a[ i ] = view_a[ i ] / gamma_i + gravity;
      }
      else {
         view_Drho[ i ] = 0.f;
         view_a[ i ] = 0.f + gravity;
      }
   };
   TNL::Algorithms::parallelFor< DeviceType >(
         fluid->getFirstActiveParticle(), fluid->getLastActiveParticle() + 1, particleLoop );

   //finalize boundary-fluid interactions
   const RealType rho0 = modelParams.rho0;

   auto view_rho_bound = boundary->variables->rho.getView();
   const auto view_gamma_bound = boundary->variables->gamma.getConstView();

   auto particleLoopBoundary = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const RealType gamma_i = view_gamma_bound[ i ];
      const RealType rho_i = view_rho_bound[ i ];

      if( gamma_i > 0.01f ) {
         view_rho_bound[ i ] = ( rho_i / gamma_i > rho0 ) ? ( rho_i / gamma_i ) : rho0;
         //printf( "Gamma_i: %f, Rho_i: %f, product: %f", gamma_i, rho_i, ( rho_i / gamma_i ) );
      }
      else {
         view_rho_bound[ i ] = rho0;
      }
   };
   TNL::Algorithms::parallelFor< DeviceType >(
         boundary->getFirstActiveParticle(), boundary->getLastActiveParticle() + 1, particleLoopBoundary );
}

template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename BoudaryPointer, typename OpenBoundaryPointer >
void
WCSPH_BI< Particles, ModelConfig >::interactWithPeriodicBoundary( FluidPointer& fluid,
                                                                  BoudaryPointer& boundary,
                                                                  OpenBoundaryPointer& openBoundary,
                                                                  ModelParams& modelParams,
                                                                  const VectorType shift )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   GlobalIndexType numberOfParticles = fluid->particles->getNumberOfParticles();
   GlobalIndexType numberOfParticles_bound = boundary->particles->getNumberOfParticles();
   const RealType searchRadius = fluid->particles->getSearchRadius();

   typename Particles::NeighborsLoopParams searchInFluid( fluid->particles );
   typename Particles::NeighborsLoopParams searchInBound( boundary->particles );

   /* CONSTANT VARIABLES */
   const RealType h = modelParams.h;
   const RealType m = modelParams.mass;
   const RealType ds = modelParams.boundaryElementSize;
   const RealType rho0 = modelParams.rho0;

   typename DiffusiveTerm::ParamsType diffusiveTermsParams( modelParams );
   typename ViscousTerm::ParamsType viscousTermsParams( modelParams );
   typename EOS::ParamsType eosParams( modelParams );
   typename BoundaryViscousTerm::ParamsType boundaryViscoParams( modelParams );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->particles->getPoints().getView();
   const auto view_rho = fluid->variables->rho.getView();
   auto view_Drho = fluid->variables->drho.getView();
   const auto view_v = fluid->variables->v.getView();
   auto view_a = fluid->variables->a.getView();
   auto view_gamma = fluid->variables->gamma.getView();

   const auto view_points_bound = boundary->particles->getPoints().getView();
   auto view_rho_bound = boundary->variables->rho.getView();
   const auto view_v_bound = boundary->variables->v.getView();
   const auto view_n_bound = boundary->variables->n.getView();

   const auto zoneParticleIndices_view = openBoundary->zone.getParticlesInZone().getConstView();
   const GlobalIndexType numberOfZoneParticles = openBoundary->zone.getNumberOfParticles();

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

         const RealType F = KernelFunction::F( drs, h );
         const RealType W = KernelFunction::W( drs, h );
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

         const RealType W = KernelFunction::W( drs, h );

         *drho_i += ( -1.f ) * ( v_ij, n_j ) * W * rho_j * ds;

         const RealType p_term = ( p_i + p_j ) / ( rho_i * rho_j );
         const RealType visco =  ViscousTerm::Pi( rho_i, rho_j, drs, ( r_ij, v_ij ), viscousTermsParams );
         *a_i += ( p_term + visco ) * n_j * W * rho_j * ds + BoundaryViscousTerm::Xi( r_ij, v_ij, n_j, boundaryViscoParams );
      }
   };

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const GlobalIndexType p = zoneParticleIndices_view[ i ];
      const VectorType r_i = view_points[ p ] + shift;
      const VectorType v_i = view_v[ p ];
      const RealType rho_i = view_rho[ p ];
      const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );
      VectorType a_i = 0.f;
      RealType drho_i = 0.f;
      RealType gamma_i = 0.f;

      TNL::ParticleSystem::NeighborsLoop::exec( p, r_i, searchInFluid, FluidFluid, v_i, rho_i, p_i, &drho_i, &a_i, &gamma_i );
      TNL::ParticleSystem::NeighborsLoopAnotherSet::exec( p, r_i, searchInBound, FluidBound, v_i, rho_i, p_i, &drho_i, &a_i );

      view_Drho[ p ] += drho_i ;
      view_a[ p ] += a_i;
      view_gamma[ p ] += gamma_i;
   };
   Algorithms::parallelFor< DeviceType >( 0, numberOfZoneParticles, particleLoop );

}

template< typename Particles, typename ModelConfig >
template< typename OpenBoundaryPointer, typename BoudaryPointer, typename FluidPointer >
void
WCSPH_BI< Particles, ModelConfig >::updateSolidBoundaryPeriodicBoundary( FluidPointer& fluid,
                                                                         BoudaryPointer& boundary,
                                                                         OpenBoundaryPointer& openBoundary,
                                                                         ModelParams& modelParams,
                                                                         const VectorType shift )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   typename Particles::NeighborsLoopParams searchInOpenBoundary( openBoundary->particles );

   /* CONSTANT VARIABLES */
   const RealType searchRadius = fluid->particles->getSearchRadius();
   const RealType h = modelParams.h;
   const RealType m = modelParams.mass;
   const RealType rho0 = modelParams.rho0;

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->particles->getPoints().getView();
   const auto view_rho = fluid->variables->rho.getView();

   const auto view_points_bound = boundary->particles->getPoints().getView();
   auto view_rho_bound = boundary->variables->rho.getView();
   auto view_gamma_bound = boundary->variables->gamma.getView();

   const auto zoneParticleIndices_view = openBoundary->zone.getParticlesInZone().getConstView();
   const GlobalIndexType numberOfZoneParticles = openBoundary->zone.getNumberOfParticles();

   auto BoundFluid = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j,
         VectorType& r_i, RealType* rho_i, RealType* gamma_i ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         const RealType rho_j = view_rho[ j ];

         const RealType W = KernelFunction::W( drs, h );

         *rho_i += W * m;
         *gamma_i += W * m / rho_j;
      }
   };

   auto particleLoopBoundary = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const GlobalIndexType p = zoneParticleIndices_view[ i ];
      const VectorType r_i = view_points_bound[ p ] + shift;

      RealType rho_i = 0.f;
      RealType gamma_i = 0.f;

      TNL::ParticleSystem::NeighborsLoopAnotherSet::exec( p, r_i, searchInOpenBoundary, BoundFluid, &rho_i, &gamma_i );

      view_rho_bound[ p ] += rho_i;
      view_gamma_bound[ p ] += gamma_i;
   };
   Algorithms::parallelFor< DeviceType >( 0, numberOfZoneParticles, particleLoopBoundary );
}

} // SPH
} // TNL
