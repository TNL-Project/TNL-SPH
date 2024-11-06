#include "BoundaryConditionsTypes.h"
#include "Interactions.h"
#include <type_traits>

namespace TNL {
namespace SPH {


template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename BoudaryPointer >
void
WCSPH_MFD< Particles, ModelConfig >::computeMomentumMatrices( FluidPointer& fluid, BoudaryPointer& boundary )
{

   auto FluidFluid = [ = ] __cuda_callable__( LocalIndexType i,
                                              LocalIndexType j,
                                              VectorType & r_i,
                                              MomentumMatrixType *M_i ) mutable
   {
      const BaseVectorType W_ji = ABFs::eval( r_ji );
      const BaseVectorType X_ji = TaylorMonomials::eval( r_ji );
      *M_i += tensorProduct( X_ji, W_ji );
   };

   auto particleLoopConsistent = [ = ] __cuda_callable__( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points[ i ];
      MomentumMatrixType M_i = 0.f;

      Particles::NeighborsLoop::exec( i, r_i, searchInFluid, FluidFluid, v_i, rho_i, p_i, &drho_i, &a_i, &gamma_i );

      momentumMatrix_view[ i ] = M_i;
   };
}

template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename BoudaryPointer >
void
WCSPH_MFD< Particles, ModelConfig >::computeABFsCoefs( const BaseVectorType& W_ji
                                                       const MomentumMatrixType& M_i,
                                                       const BaseVectorType& Cd )
{
   return ( W_ji )
}


template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename BoudaryPointer >
void
WCSPH_BI< Particles, ModelConfig >::interaction( FluidPointer& fluid, BoudaryPointer& boundary, ModelParams& modelParams )
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
   const auto view_elementSize_bound = boundary->variables->elementSize.getConstView();

   auto FluidFluid = [ = ] __cuda_callable__( LocalIndexType i,
                                              LocalIndexType j,
                                              VectorType & r_i,
                                              VectorType & v_i,
                                              RealType & rho_i,
                                              RealType & p_i,
                                              RealType * drho_i,
                                              VectorType * a_i,
                                              RealType * gamma_i ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius ) {
         const VectorType v_j = view_v[ j ];
         const RealType rho_j = view_rho[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );

         const VectorType v_ij = v_i - v_j;

         const BaseVectorType W_ji = ABFs::eval( r_ji );
         const VectorType w_ji_grad = { ( W_ji, psi_i_d1x ), ( W_ji, psi_i_d1y ) };
         const RealType w_ji_lap = W_ji * psi_i_d2x + W_ji * psi_i_2y;


         const RealType V_j = m / rho_j;

         //const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, r_ij, drs, diffusiveTermsParams );
         //const RealType diffTerm = psi * ( r_ij, gradW ) * V_j;
         *drho_i += rho_i * ( v_j - v_i, w_ji_grad ); // * V_j - diffTerm;

         const VectorType grad_p = ( p_j - p_i ) * w_ji_grad;
         const VectorType visco_term = ( mu / rho_i ) * ( v_j - v_i ) * w_ji_lap;

         *a_i += ( -1.0f / rho_i ) * grad_p + visco_term;

         *gamma_i += W * m / rho_j;
      }
   };

   auto FluidBoundConsistent = [ = ] __cuda_callable__( LocalIndexType i,
                                                        LocalIndexType j,
                                                        VectorType & r_i,
                                                        VectorType & v_i,
                                                        RealType & rho_i,
                                                        RealType & p_i,
                                                        RealType * drho_i,
                                                        VectorType * a_i ) mutable
   {
      const VectorType r_j = view_points_bound[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius ) {
         const VectorType v_j = view_v_bound[ j ];
         const RealType rho_j = view_rho_bound[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );
         const VectorType n_j = view_n_bound[ j ];
         const RealType ds_j = view_elementSize_bound[ j ];

         const VectorType v_ij = v_i - v_j;

         const BaseVectorType W_ji =
         const VectorType

         *drho_i += ( -1.f ) * ( v_ij, n_j ) * W * rho_j * ds_j;

         const VectorType grad_p = ( p_j - p_i ) * w_ji_grad;

         const VectorType visco_term = ViscousTerm::BI_Pi( drs, r_ij, v_ij, rho_i, rho_j, W, n_j, ds_j, viscousTermsParams );
         const VectorType bvt = BoundaryViscousTerm::Xi( r_ij, v_ij, n_j, boundaryViscoParams );
         //FIXME: The signs are fucked, because I used inner normals.
         //       Correct is of course: -1/rho * grad + visco
         *a_i += ( 1.f / rho_i ) * grad_p - visco_term  + bvt;
      }
   };

   auto particleLoopConsistent = [ = ] __cuda_callable__( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points[ i ];
      const VectorType v_i = view_v[ i ];
      const RealType rho_i = view_rho[ i ];
      const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );
      const MomentumMatrixType M_i = momentumMatrix_view[ i ];
      VectorType dvdt_i = 0.f;
      RealType drhodt_i = 0.f;

      // get gradient ABFs coefficients
      const BaseVectorType psi_i_d1x = MFD::psi_dx( M_i, 1 );
      const BaseVectorType psi_i_d1y = MFD::psi_dy( M_i, 1 );
      // get laplace ABFs coefficients
      const BaseVectorType psi_i_d2x = MFD::psi_dx( M_i, 2 );
      const BaseVectorType psi_i_d2y = MFD::psi_dy( M_i, 2 );

      Particles::NeighborsLoop::exec( i, r_i, searchInFluid, FluidFluid, v_i, rho_i, p_i, &drho_i, &a_i, &gamma_i );
      Particles::NeighborsLoopAnotherSet::exec( i, r_i, searchInBound, FluidBoundConsistent, v_i, rho_i, p_i, &drho_i, &a_i );

      view_drhodt[ i ] = drhodt_i;
      view_dvdt[ i ] = dvdt_i;
   };
   fluid->particles->forAll( particleLoopConsistent );
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

   auto BoundFluid = [ = ] __cuda_callable__( LocalIndexType i,
                                              LocalIndexType j,
                                              VectorType& r_i,
                                              RealType*
                                              rho_i,
                                              RealType* gamma_i ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius ) {
         const RealType rho_j = view_rho[ j ];

         const RealType W = KernelFunction::W( drs, h );

         *rho_i += W * m;
         *gamma_i += W * m / rho_j;
      }
   };

   auto BoundFluidConservative = [ = ] __cuda_callable__( LocalIndexType i,
                                                          LocalIndexType j,
                                                          VectorType& r_i,
                                                          RealType*
                                                          rho_i,
                                                          RealType* gamma_i ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius ) {
         const RealType rho_j = view_rho[ j ];

         const RealType W = KernelFunction::W( drs, h );

         *rho_i += 2.f * W * m;
         *gamma_i += W * m / rho_j;
      }
   };

   auto particleLoopBoundary = [ = ] __cuda_callable__( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points_bound[ i ];
      RealType rho_i = 0.f;
      RealType gamma_i = 0.f;

      //Particles::NeighborsLoopAnotherSet::exec( i, r_i, searchInFluid, BoundFluid, &rho_i, &gamma_i );
      Particles::NeighborsLoopAnotherSet::exec( i, r_i, searchInFluid, BoundFluidConservative, &rho_i, &gamma_i );

      view_rho_bound[ i ] = rho_i;
      view_gamma_bound[ i ] = gamma_i;
   };
   boundary->particles->forAll( particleLoopBoundary );

   if constexpr( Model::ModelConfigType::SPHConfig::numberOfPeriodicBuffers > 0 ) {
      for( long unsigned int i = 0; i < std::size( boundary->periodicPatches ); i++ ) {
         const auto zoneParticleIndices_view = boundary->periodicPatches[ i ]->particleZone.getParticlesInZone().getConstView();
         const GlobalIndexType numberOfZoneParticles = boundary->periodicPatches[ i ]->particleZone.getNumberOfParticles();
         const VectorType shift = boundary->periodicPatches[ i ]->config.shift;

         auto periodicParticleLoopBoundary = [ = ] __cuda_callable__( LocalIndexType i ) mutable
         {
            const GlobalIndexType p = zoneParticleIndices_view[ i ];
            const VectorType r_i = view_points_bound[ p ] + shift;
            RealType rho_i = 0.f;
            RealType gamma_i = 0.f;

            Particles::NeighborsLoopAnotherSet::exec( p, r_i, searchInFluid, BoundFluid, &rho_i, &gamma_i );

            view_rho_bound[ p ] += rho_i;
            view_gamma_bound[ p ] += gamma_i;
         };
         Algorithms::parallelFor< DeviceType >( 0, numberOfZoneParticles, periodicParticleLoopBoundary );
      }
   }
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

   auto BoundOpenBoundary =
      [ = ] __cuda_callable__(
         LocalIndexType i, LocalIndexType j, VectorType & r_i, RealType * rho_i, RealType * gamma_i ) mutable
   {
      const VectorType r_j = view_points_openBound[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius ) {
         const RealType rho_j = view_rho_openBound[ j ];

         const RealType W = KernelFunction::W( drs, h );

         *rho_i += W * m;
         *gamma_i += W * m / rho_j;
      }
   };

   auto particleLoopBoundary = [ = ] __cuda_callable__( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points_bound[ i ];

      RealType rho_i = 0.f;
      RealType gamma_i = 0.f;

      Particles::NeighborsLoopAnotherSet::exec( i, r_i, searchInOpenBoundary, BoundOpenBoundary, &rho_i, &gamma_i );

      view_rho_bound[ i ] += rho_i;
      view_gamma_bound[ i ] += gamma_i;
   };
   boundary->particles->forAll( particleLoopBoundary );
}

template< typename Particles, typename ModelConfig >
template< typename EquationOfState, typename PhysicalObjectPointer >
void
WCSPH_BI< Particles, ModelConfig >::computePressureFromDensity( PhysicalObjectPointer& physicalObject,
                                                                ModelParams& modelParams )
{
   auto view_rho = physicalObject->getVariables()->rho.getView();
   auto view_p = physicalObject->getVariables()->p.getView();

   typename EOS::ParamsType eosParams( modelParams );

   auto evalPressure = [ = ] __cuda_callable__( int i ) mutable
   {
      view_p[ i ] = EquationOfState::DensityToPressure( view_rho[ i ], eosParams );
   };
   physicalObject->particles->forAll( evalPressure ); //TODO: forloop?
}

template< typename Particles, typename ModelConfig >
template< typename FluidPointer >
void
WCSPH_BI< Particles, ModelConfig >::filterDensity( FluidPointer& fluid, ModelParams& modelParams )
{
   //TODO: This requires ParticleType template due to neighbor loop. I don't like this.
   DensityFilter::template filterDensityOverlaps< ParticlesType >( fluid, modelParams );
   //DensityFilter::template filterDensity< ParticlesType >( fluid, modelParams );
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

   auto FluidOpenBoundary = [ = ] __cuda_callable__( LocalIndexType i,
                                                     LocalIndexType j,
                                                     VectorType & r_i,
                                                     VectorType & v_i,
                                                     RealType & rho_i,
                                                     RealType & p_i,
                                                     RealType * drho_i,
                                                     VectorType * a_i,
                                                     RealType * gamma_i ) mutable
   {
      const VectorType r_j = view_points_openBound[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius ) {
         const VectorType v_j = view_v_openBound[ j ];
         const RealType rho_j = view_rho_openBound[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );

         const VectorType v_ij = v_i - v_j;

         const RealType F = KernelFunction::F( drs, h );
         const RealType W = KernelFunction::W( drs, h );
         const VectorType gradW = r_ij * F;
         const RealType V_j = m / rho_j;

         const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, r_ij, drs, diffusiveTermsParams );
         const RealType diffTerm = psi * ( r_ij, gradW ) * V_j;
         *drho_i += rho_i * ( v_ij, gradW ) * V_j - diffTerm;

         const VectorType grad_p = ( p_i + p_j ) * gradW * V_j;
         const VectorType visco_term = ViscousTerm::Pi( drs, r_ij, v_ij, rho_i, rho_j, gradW, V_j, viscousTermsParams );
         *a_i += ( -1.0f / rho_i ) * grad_p + visco_term;

         *gamma_i += W * m / rho_j;
      }
   };

   auto particleLoop = [ = ] __cuda_callable__( LocalIndexType i ) mutable
   {
      const GlobalIndexType p = zoneParticleIndices_view[ i ];
      const VectorType r_i = view_points[ p ];
      const VectorType v_i = view_v[ p ];
      const RealType rho_i = view_rho[ p ];
      const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );

      VectorType a_i = 0.f;
      RealType drho_i = 0.f;
      RealType gamma_i = 0.f;

      Particles::NeighborsLoopAnotherSet::exec(
         p, r_i, searchInOpenBoundary, FluidOpenBoundary, v_i, rho_i, p_i, &drho_i, &a_i, &gamma_i );

      view_Drho[ p ] += drho_i;
      view_a[ p ] += a_i;
      view_gamma[ p ] += gamma_i;
   };
   Algorithms::parallelFor< DeviceType >( 0, numberOfZoneParticles, particleLoop );
}

//FIXME: WTF is this function
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

   auto FluidBound = [ = ] __cuda_callable__( LocalIndexType i,
                                              LocalIndexType j,
                                              VectorType & r_i,
                                              VectorType & v_i,
                                              RealType & rho_i,
                                              RealType & p_i,
                                              RealType * drho_i,
                                              VectorType * a_i ) mutable
   {
      const VectorType r_j = view_points_openBound[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius ) {
         const VectorType v_j = view_v_openBound[ j ];
         const RealType rho_j = view_rho_openBound[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );
         const VectorType n_j = view_n_bound[ j ];

         /* Interaction: */
         const VectorType v_ij = v_i - v_j;

         const RealType W = KernelFunction::W( drs, h );

         *drho_i += ( -1.f ) * ( v_ij, n_j ) * W * rho_j * ds;

         const RealType p_term = ( p_i + p_j ) / ( rho_i * rho_j );
         const RealType visco = ViscousTerm::Pi( rho_i, rho_j, drs, ( r_ij, v_ij ), viscousTermsParams );
         *a_i += ( p_term + visco ) * n_j * W * rho_j * ds + BoundaryViscousTerm::Xi( r_ij, v_ij, n_j, boundaryViscoParams );
      }
   };

   auto particleLoop = [ = ] __cuda_callable__( LocalIndexType i ) mutable
   {
      const GlobalIndexType p = zoneParticleIndices_view[ i ];
      const VectorType r_i = view_points[ p ];
      const VectorType v_i = view_v[ p ];
      const RealType rho_i = view_rho[ p ];
      const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );

      VectorType a_i = 0.f;
      RealType drho_i = 0.f;

      Particles::NeighborsLoopAnotherSet::exec( p, r_i, searchInOpenBoundary, FluidBound, v_i, rho_i, p_i, &drho_i, &a_i );

      view_Drho[ p ] += drho_i;
      view_a[ p ] += a_i;
   };
   Algorithms::parallelFor< DeviceType >( 0, numberOfZoneParticles, particleLoop );
}

template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename BoundaryPointer >
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

   auto finalizeInteractionConsistent = [ = ] __cuda_callable__( LocalIndexType i ) mutable
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

   auto finalizeInteractionConservative = [ = ] __cuda_callable__( LocalIndexType i ) mutable
   {
      view_a[ i ] += gravity;
   };

   if constexpr( std::is_same_v< typename ModelConfig::BCType, WCSPH_BCTypes::BIConsistent_numeric> )
      fluid->particles->forAll( finalizeInteractionConsistent );
   else if constexpr( std::is_same_v< typename ModelConfig::BCType, WCSPH_BCTypes::BIConservative_numeric> )
      fluid->particles->forAll( finalizeInteractionConservative );
}

template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename BoundaryPointer >
void
WCSPH_BI< Particles, ModelConfig >::finalizeBoundaryInteraction( FluidPointer& fluid,
                                                                 BoundaryPointer& boundary,
                                                                 ModelParams& modelParams )
{
   //finalize boundary-fluid interactions
   const RealType rho0 = modelParams.rho0;

   auto view_rho_bound = boundary->variables->rho.getView();
   const auto view_gamma_bound = boundary->variables->gamma.getConstView();

   auto particleLoopBoundary = [ = ] __cuda_callable__( LocalIndexType i ) mutable
   {
      //const RealType gamma_i = view_gamma_bound[ i ];
      const RealType gamma_i = 1;
      const RealType rho_i = view_rho_bound[ i ];

      if( gamma_i > 0.01f ) {
         view_rho_bound[ i ] = ( rho_i / gamma_i > rho0 ) ? ( rho_i / gamma_i ) : rho0;
      }
      else {
         view_rho_bound[ i ] = rho0;
      }
   };
   boundary->particles->forAll( particleLoopBoundary );
}

}  //namespace SPH
}  //namespace TNL
