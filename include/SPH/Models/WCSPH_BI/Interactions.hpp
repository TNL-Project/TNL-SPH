#include "BoundaryConditionsTypes.h"
#include "Interactions.h"
#include <type_traits>

namespace TNL {
namespace SPH {

template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename BoudaryPointer >
void
WCSPH_BI< Particles, ModelConfig >::interaction( FluidPointer& fluid, BoudaryPointer& boundary, ModelParams& modelParams )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   GlobalIndexType numberOfParticles = fluid->getParticles()->getNumberOfParticles();
   GlobalIndexType numberOfParticles_bound = boundary->getParticles()->getNumberOfParticles();
   const RealType searchRadius = fluid->getParticles()->getSearchRadius();

   typename Particles::NeighborsLoopParams searchInFluid( fluid->getParticles() );
   typename Particles::NeighborsLoopParams searchInBound( boundary->getParticles() );

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
   const auto view_points = fluid->getParticles()->getPoints().getView();
   const auto view_rho = fluid->getVariables()->rho.getView();
   auto view_Drho = fluid->getVariables()->drho.getView();
   const auto view_v = fluid->getVariables()->v.getView();
   auto view_a = fluid->getVariables()->a.getView();
   auto view_gamma = fluid->getVariables()->gamma.getView();

   const auto view_points_bound = boundary->getParticles()->getPoints().getView();
   auto view_rho_bound = boundary->getVariables()->rho.getView();
   const auto view_v_bound = boundary->getVariables()->v.getView();
   const auto view_n_bound = boundary->getVariables()->n.getView();
   const auto view_elementSize_bound = boundary->getVariables()->elementSize.getConstView();

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

         // there is reason why is the gradient and volume merged and written this way
         const VectorType gradWV_j = r_ij * KernelFunction::F( drs, h ) * m / rho_j;

         const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, r_ij, drs, diffusiveTermsParams );
         const RealType diffTerm = psi * ( r_ij, gradWV_j );
         *drho_i += rho_i * ( v_ij, gradWV_j ) - diffTerm;

         const VectorType grad_p = ( p_i + p_j ) * gradWV_j;
         const VectorType visco_term = ViscousTerm::Pi( drs, r_ij, v_ij, rho_i, rho_j, gradWV_j, viscousTermsParams );
         *a_i += ( -1.0f / rho_i ) * grad_p + visco_term;

         *gamma_i +=  KernelFunction::W( drs, h ) * m / rho_j;
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
         const RealType S_j = view_elementSize_bound[ j ];
         const VectorType v_ij = v_i - v_j;

         // there is reason why is the kernel and area merged and written this way
         const RealType WS_j = S_j * KernelFunction::W( drs, h ) ;

         *drho_i += ( -1.f ) * rho_j * ( v_ij, n_j ) * WS_j;

         const VectorType grad_p = ( p_i + p_j ) * n_j * WS_j;
         const VectorType visco_term = ViscousTerm::BI_Pi( drs, r_ij, v_ij, rho_i, rho_j, n_j, WS_j, viscousTermsParams );
         const VectorType bvt = BoundaryViscousTerm::Xi( r_ij, v_ij, n_j, boundaryViscoParams );
         //FIXME: The signs are fucked, because I used inner normals.
         //       Correct is of course: -1/rho * grad + visco
         *a_i += ( 1.f / rho_i ) * grad_p - visco_term  + bvt;
      }
   };

   auto FluidBoundConservative = [ = ] __cuda_callable__( LocalIndexType i,
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
         const VectorType n_j =  view_n_bound[ j ];
         const RealType S_j = view_elementSize_bound[ j ];
         const VectorType v_ij = v_i - v_j;

         // there is reason why is the kernel and area merged and written this way
         const RealType WS_j = S_j * KernelFunction::W( drs, h ) ;

         //FIXME: The signs are fucked, because I used inner normals.
         //       Correct is of course: -1/rho * grad + visco
         *drho_i += ( -2.f ) * rho_i * ( v_i - v_j, n_j ) * WS_j;

         const VectorType grad_p = 2.f * p_i * n_j * WS_j;

         const VectorType visco_term = ViscousTerm::BI_Pi( drs, r_ij, v_ij, rho_i, rho_j, n_j, WS_j, viscousTermsParams );
         const VectorType bvt = BoundaryViscousTerm::Xi( r_ij, v_ij, n_j, boundaryViscoParams );

         //FIXME: The signs are fucked, because I used inner normals.
         //       Correct is of course: -1/rho * grad + visco
         *a_i += ( +1.f / rho_i ) * grad_p - visco_term  + bvt;
      }
   };

   auto particleLoopConsistent = [ = ] __cuda_callable__( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points[ i ];
      const VectorType v_i = view_v[ i ];
      const RealType rho_i = view_rho[ i ];
      const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );
      VectorType a_i = 0.f;
      RealType drho_i = 0.f;
      RealType gamma_i = 0.f;

      Particles::NeighborsLoop::exec( i, r_i, searchInFluid, FluidFluid, v_i, rho_i, p_i, &drho_i, &a_i, &gamma_i );
      Particles::NeighborsLoopAnotherSet::exec( i, r_i, searchInBound, FluidBoundConsistent, v_i, rho_i, p_i, &drho_i, &a_i );

      view_Drho[ i ] = drho_i;
      view_a[ i ] = a_i;
      view_gamma[ i ] = gamma_i;
   };

   auto particleLoopConservative = [ = ] __cuda_callable__( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points[ i ];
      const VectorType v_i = view_v[ i ];
      const RealType rho_i = view_rho[ i ];
      const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );
      VectorType a_i = 0.f;
      RealType drho_i = 0.f;
      RealType gamma_i = 0.f;

      Particles::NeighborsLoop::exec( i, r_i, searchInFluid, FluidFluid, v_i, rho_i, p_i, &drho_i, &a_i, &gamma_i );
      Particles::NeighborsLoopAnotherSet::exec( i, r_i, searchInBound, FluidBoundConservative, v_i, rho_i, p_i, &drho_i, &a_i );

      view_Drho[ i ] = drho_i;
      view_a[ i ] = a_i;
      view_gamma[ i ] = gamma_i;
   };

   if constexpr( std::is_same_v< typename ModelConfig::BCType, WCSPH_BCTypes::BIConsistent_numeric > ||
                 std::is_same_v< typename ModelConfig::BCType, WCSPH_BCTypes::BIConsistent_numeric_interpolated >)
      fluid->getParticles()->forAll( particleLoopConsistent );
   else if constexpr( std::is_same_v< typename ModelConfig::BCType, WCSPH_BCTypes::BIConservative_numeric> )
      fluid->getParticles()->forAll( particleLoopConservative );

   if( fluid->periodicPatches.size() > 0 ) {
      for( long unsigned int i = 0; i < std::size( fluid->periodicPatches ); i++ ) {
         const auto zoneParticleIndices_view = fluid->periodicPatches[ i ]->particleZone.getParticlesInZone().getConstView();
         const GlobalIndexType numberOfZoneParticles = fluid->periodicPatches[ i ]->particleZone.getNumberOfParticles();
         const VectorType shift = fluid->periodicPatches[ i ]->config.shift;

         auto periodicParticleLoop = [ = ] __cuda_callable__( LocalIndexType i ) mutable
         {
            const GlobalIndexType p = zoneParticleIndices_view[ i ];
            const VectorType r_i = view_points[ p ] + shift;
            const VectorType v_i = view_v[ p ];
            const RealType rho_i = view_rho[ p ];
            const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );
            VectorType a_i = 0.f;
            RealType drho_i = 0.f;
            RealType gamma_i = 0.f;

            Particles::NeighborsLoop::exec( p, r_i, searchInFluid, FluidFluid, v_i, rho_i, p_i, &drho_i, &a_i, &gamma_i );
            Particles::NeighborsLoopAnotherSet::exec( p, r_i, searchInBound, FluidBoundConsistent, v_i, rho_i, p_i, &drho_i, &a_i );

            view_Drho[ p ] += drho_i;
            view_a[ p ] += a_i;
            view_gamma[ p ] += gamma_i;
         };
         Algorithms::parallelFor< DeviceType >( 0, numberOfZoneParticles, periodicParticleLoop );
      }
   }
}

template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename BoudaryPointer >
requires std::is_same_v< typename ModelConfig::BCType, WCSPH_BCTypes::BIConsistent_numeric > ||
         std::is_same_v< typename ModelConfig::BCType, WCSPH_BCTypes::BIConservative_numeric >
void
WCSPH_BI< Particles, ModelConfig >::updateSolidBoundary( FluidPointer& fluid,
                                                         BoudaryPointer& boundary,
                                                         ModelParams& modelParams )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   typename Particles::NeighborsLoopParams searchInFluid( fluid->getParticles() );

   /* CONSTANT VARIABLES */
   const RealType searchRadius = fluid->getParticles()->getSearchRadius();
   const RealType h = modelParams.h;
   const RealType m = modelParams.mass;
   const RealType rho0 = modelParams.rho0;

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->getParticles()->getPoints().getView();
   const auto view_rho = fluid->getVariables()->rho.getView();

   const auto view_points_bound = boundary->getParticles()->getPoints().getView();
   auto view_rho_bound = boundary->getVariables()->rho.getView();
   auto view_gamma_bound = boundary->getVariables()->gamma.getView();

   auto BoundFluidConsistent = [ = ] __cuda_callable__( LocalIndexType i,
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

   auto particleLoopBoundaryConsistent = [ = ] __cuda_callable__( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points_bound[ i ];
      RealType rho_i = 0.f;
      RealType gamma_i = 0.f;

      Particles::NeighborsLoopAnotherSet::exec( i, r_i, searchInFluid, BoundFluidConsistent, &rho_i, &gamma_i );

      view_rho_bound[ i ] = rho_i;
      view_gamma_bound[ i ] = gamma_i;
   };

   auto particleLoopBoundaryConservative = [ = ] __cuda_callable__( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points_bound[ i ];
      RealType rho_i = 0.f;
      RealType gamma_i = 0.f;

      Particles::NeighborsLoopAnotherSet::exec( i, r_i, searchInFluid, BoundFluidConservative, &rho_i, &gamma_i );

      view_rho_bound[ i ] = rho_i;
      view_gamma_bound[ i ] = gamma_i;
   };

   if constexpr( std::is_same_v< typename ModelConfig::BCType, WCSPH_BCTypes::BIConsistent_numeric> )
      boundary->getParticles()->forAll( particleLoopBoundaryConsistent );
   else if constexpr( std::is_same_v< typename ModelConfig::BCType, WCSPH_BCTypes::BIConservative_numeric> )
      boundary->getParticles()->forAll( particleLoopBoundaryConservative );

   if( boundary->periodicPatches.size() > 0 ) {
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

            Particles::NeighborsLoopAnotherSet::exec( p, r_i, searchInFluid, BoundFluidConsistent, &rho_i, &gamma_i );

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
   typename Particles::NeighborsLoopParams searchInOpenBoundary( openBoundary->getParticles() );

   /* CONSTANT VARIABLES */
   const RealType searchRadius = openBoundary->getParticles()->getSearchRadius();
   const RealType h = modelParams.h;
   const RealType m = modelParams.mass;
   const RealType rho0 = modelParams.rho0;

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points_openBound = openBoundary->getParticles()->getPoints().getView();
   const auto view_rho_openBound = openBoundary->getVariables()->rho.getView();

   const auto view_points_bound = boundary->getParticles()->getPoints().getView();
   auto view_rho_bound = boundary->getVariables()->rho.getView();
   auto view_gamma_bound = boundary->getVariables()->gamma.getView();

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
   boundary->getParticles()->forAll( particleLoopBoundary );
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
   physicalObject->getParticles()->forAll( evalPressure ); //TODO: forloop?
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
   typename Particles::NeighborsLoopParams searchInFluid( fluid->getParticles() );
   typename Particles::NeighborsLoopParams searchInOpenBoundary( openBoundary->getParticles() );

   /* CONSTANT VARIABLES */
   const RealType searchRadius = fluid->getParticles()->getSearchRadius();
   const RealType h = modelParams.h;
   const RealType m = modelParams.mass;

   typename DiffusiveTerm::ParamsType diffusiveTermsParams( modelParams );
   typename ViscousTerm::ParamsType viscousTermsParams( modelParams );
   typename EOS::ParamsType eosParams( modelParams );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->getParticles()->getPoints().getView();
   const auto view_rho = fluid->getVariables()->rho.getView();
   auto view_Drho = fluid->getVariables()->drho.getView();
   const auto view_v = fluid->getVariables()->v.getView();
   auto view_a = fluid->getVariables()->a.getView();
   auto view_gamma = fluid->getVariables()->gamma.getView();

   const auto view_points_openBound = openBoundary->getParticles()->getPoints().getView();
   auto view_rho_openBound = openBoundary->getVariables()->rho.getView();
   auto view_v_openBound = openBoundary->getVariables()->v.getView();

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

         // there is reason why is the gradient and volume merged and written this way
         const VectorType gradWV_j = r_ij * KernelFunction::F( drs, h ) * m / rho_j;

         const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, r_ij, drs, diffusiveTermsParams );
         const RealType diffTerm = psi * ( r_ij, gradWV_j );
         *drho_i += rho_i * ( v_ij, gradWV_j ) - diffTerm;

         const VectorType grad_p = ( p_i + p_j ) * gradWV_j;
         const VectorType visco_term = ViscousTerm::Pi( drs, r_ij, v_ij, rho_i, rho_j, gradWV_j, viscousTermsParams );
         *a_i += ( -1.0f / rho_i ) * grad_p + visco_term;

         *gamma_i +=  KernelFunction::W( drs, h ) * m / rho_j;
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
   typename Particles::NeighborsLoopParams searchInOpenBoundary( openBoundary->getParticles() );

   /* CONSTANT VARIABLES */
   const RealType searchRadius = fluid->getParticles()->getSearchRadius();
   const RealType h = modelParams.h;
   const RealType ds = modelParams.boundaryElementSize;
   const RealType m = modelParams.mass;
   typename DiffusiveTerm::ParamsType diffusiveTermsParams( modelParams );
   typename ViscousTerm::ParamsType viscousTermsParams( modelParams );
   typename EOS::ParamsType eosParams( modelParams );
   typename BoundaryViscousTerm::ParamsType boundaryViscoParams( modelParams );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->getParticles()->getPoints().getView();
   const auto view_rho = fluid->getVariables()->rho.getView();
   auto view_Drho = fluid->getVariables()->drho.getView();
   const auto view_v = fluid->getVariables()->v.getView();
   auto view_a = fluid->getVariables()->a.getView();
   auto view_gamma = fluid->getVariables()->gamma.getView();
   const auto view_points_openBound = openBoundary->getParticles()->getPoints().getView();
   auto view_rho_openBound = openBoundary->getVariables()->rho.getView();
   auto view_v_openBound = openBoundary->getVariables()->v.getView();
   const auto view_n_bound = openBoundary->getVariables()->n.getView();

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
         const VectorType v_ij = v_i - v_j;

         // there is reason why is the kernel and area merged and written this way
         const RealType WS_j = ds * KernelFunction::W( drs, h ) ;

         *drho_i += ( -1.f ) * rho_j * ( v_ij, n_j ) * WS_j;

         const RealType p_term = ( p_i + p_j ) / ( rho_i * rho_j );
         const RealType visco = ViscousTerm::Pi( rho_i, rho_j, drs, ( r_ij, v_ij ), viscousTermsParams );
         *a_i += ( p_term + visco ) * n_j * rho_j * WS_j + BoundaryViscousTerm::Xi( r_ij, v_ij, n_j, boundaryViscoParams );

         /*
         const VectorType grad_p = ( p_i + p_j ) * n_j * WS_j;
         const VectorType visco_term = ViscousTerm::BI_Pi( drs, r_ij, v_ij, rho_i, rho_j, n_j, WS_j, viscousTermsParams );
         const VectorType bvt = BoundaryViscousTerm::Xi( r_ij, v_ij, n_j, boundaryViscoParams );
         //FIXME: The signs are fucked, because I used inner normals.
         //       Correct is of course: -1/rho * grad + visco
         *a_i += ( 1.f / rho_i ) * grad_p - visco_term  + bvt;
         */
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

   auto view_Drho = fluid->getVariables()->drho.getView();
   auto view_a = fluid->getVariables()->a.getView();
   auto view_gamma = fluid->getVariables()->gamma.getView();

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
      fluid->getParticles()->forAll( finalizeInteractionConsistent );
   else if constexpr( std::is_same_v< typename ModelConfig::BCType, WCSPH_BCTypes::BIConservative_numeric> )
      fluid->getParticles()->forAll( finalizeInteractionConservative );
}

template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename BoundaryPointer >
requires std::is_same_v< typename ModelConfig::BCType, WCSPH_BCTypes::BIConsistent_numeric >
void
WCSPH_BI< Particles, ModelConfig >::finalizeBoundaryInteraction( FluidPointer& fluid,
                                                                 BoundaryPointer& boundary,
                                                                 ModelParams& modelParams )
{
   //finalize boundary-fluid interactions
   const RealType rho0 = modelParams.rho0;

   auto view_rho_bound = boundary->getVariables()->rho.getView();
   const auto view_gamma_bound = boundary->getVariables()->gamma.getConstView();

   auto particleLoopBoundary = [ = ] __cuda_callable__( LocalIndexType i ) mutable
   {
      const RealType gamma_i = view_gamma_bound[ i ];
      const RealType rho_i = view_rho_bound[ i ];

      if( gamma_i > 0.01f ) {
         view_rho_bound[ i ] = ( rho_i / gamma_i > rho0 ) ? ( rho_i / gamma_i ) : rho0;
      }
      else {
         view_rho_bound[ i ] = ( rho_i > rho0 ) ? rho_i : rho0;
      }

   };
   boundary->getParticles()->forAll( particleLoopBoundary );
}

template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename BoundaryPointer >
requires std::is_same_v< typename ModelConfig::BCType, WCSPH_BCTypes::BIConservative_numeric >
void
WCSPH_BI< Particles, ModelConfig >::finalizeBoundaryInteraction( FluidPointer& fluid,
                                                                 BoundaryPointer& boundary,
                                                                 ModelParams& modelParams )
{}

}  //namespace SPH
}  //namespace TNL
