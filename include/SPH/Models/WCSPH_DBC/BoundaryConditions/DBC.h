#include "../Interactions.h"

namespace TNL {
namespace SPH {

template< typename Particles, typename ModelConfig >
template< typename FluidPointer,
          typename BoudaryPointer,
          typename BCType,
          std::enable_if_t< std::is_same_v< BCType, WCSPH_BCTypes::DBC >, bool > Enabled >
void
WCSPH_DBC< Particles, ModelConfig >::updateSolidBoundary( FluidPointer& fluid,
                                                          BoudaryPointer& boundary,
                                                          ModelParams& modelParams )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   //typename Particles::NeighborsLoopParams searchInFluid( fluid->particles );
   auto searchInFluid = boundary->getParticles()->getSearchToken( fluid->getParticles() );

   /* CONSTANT VARIABLES */
   const RealType searchRadius = fluid->particles->getSearchRadius();
   const RealType h = modelParams.h;
   const RealType m = modelParams.mass;
   typename DiffusiveTerm::ParamsType diffusiveTermsParams( modelParams );
   typename EOS::ParamsType eosParams( modelParams );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->particles->getPoints().getView();
   const auto view_rho = fluid->variables->rho.getView();
   const auto view_v = fluid->variables->v.getView();

   const auto view_points_bound = boundary->particles->getPoints().getView();
   const auto view_rho_bound = boundary->variables->rho.getView();
   auto view_Drho_bound = boundary->variables->drho.getView();
   const auto view_v_bound = boundary->variables->v.getView();

   auto BoundFluid = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j,
         VectorType& r_i, VectorType& v_i, RealType& rho_i, RealType& p_i, RealType* drho_i ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         const VectorType v_j = view_v[ j ];
         const RealType rho_j = view_rho[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );

         /* Interaction */
         const VectorType v_ij = v_i - v_j;

         const RealType F = KernelFunction::F( drs, h );
         const VectorType gradW = r_ij * F;

         const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, drs, diffusiveTermsParams );
         const RealType diffTerm =  psi * ( r_ij, gradW ) * m / rho_j;
         *drho_i += ( v_ij, gradW ) * m - diffTerm;
      }
   };

   auto particleLoopBoundary = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points_bound[ i ];
      const VectorType v_i = view_v_bound[ i ];
      const RealType rho_i = view_rho_bound[ i ];
      const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );

      RealType drho_i = 0.f;
      Particles::NeighborsLoopAnotherSet::exec( i, r_i, searchInFluid, BoundFluid, v_i, rho_i, p_i, &drho_i );
      view_Drho_bound[ i ] = drho_i;
   };
   boundary->particles->forAll( particleLoopBoundary );

   if constexpr( Model::ModelConfigType::SPHConfig::numberOfPeriodicBuffers > 0 ){
      for( long unsigned int i = 0; i < std::size( boundary->periodicPatches ); i++ ){

         const auto zoneParticleIndices_view = boundary->periodicPatches[ i ]->particleZone.getParticlesInZone().getConstView();
         const GlobalIndexType numberOfZoneParticles = boundary->periodicPatches[ i ]->particleZone.getNumberOfParticles();
         const VectorType shift = boundary->periodicPatches[ i ]->config.shift;

         auto periodicParticleLoopBoundary = [=] __cuda_callable__ ( LocalIndexType i ) mutable
         {
            const GlobalIndexType p = zoneParticleIndices_view[ i ];
            const VectorType r_i = view_points_bound[ p ] + shift;
            const VectorType v_i = view_v_bound[ p ];
            const RealType rho_i = view_rho_bound[ p ];
            const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );

            RealType drho_i = 0.f;
            Particles::NeighborsLoop::exec( p, r_i, searchInFluid, BoundFluid, v_i, rho_i, p_i, &drho_i );
            view_Drho_bound[ p ] += drho_i;
         };

         Algorithms::parallelFor< DeviceType >( 0, numberOfZoneParticles, periodicParticleLoopBoundary );
      }
   }
}

template< typename Particles, typename ModelConfig >
template< typename OpenBoundaryPointer,
          typename BoudaryPointer,
          typename BCType,
          std::enable_if_t< std::is_same_v< BCType, WCSPH_BCTypes::DBC >, bool > Enabled >
void
WCSPH_DBC< Particles, ModelConfig >::updateSolidBoundaryOpenBoundary( BoudaryPointer& boundary,
                                                                      OpenBoundaryPointer& openBoundary,
                                                                      ModelParams& modelParams )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   typename Particles::NeighborsLoopParams searchInOpenBoundary( openBoundary->particles );

   /* CONSTANT VARIABLES */
   const RealType searchRadius = openBoundary->particles->getSearchRadius();
   const RealType h = modelParams.h;
   const RealType m = modelParams.mass;
   typename DiffusiveTerm::ParamsType diffusiveTermsParams( modelParams );
   typename EOS::ParamsType eosParams( modelParams );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points_openBound = openBoundary->particles->getPoints().getView();
   const auto view_rho_openBound = openBoundary->variables->rho.getView();
   const auto view_v_openBound = openBoundary->variables->v.getView();

   const auto view_points_bound = boundary->particles->getPoints().getView();
   const auto view_rho_bound = boundary->variables->rho.getView();
   auto view_Drho_bound = boundary->variables->drho.getView();
   const auto view_v_bound = boundary->variables->v.getView();

   auto BoundFluid = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j,
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

         const RealType F = KernelFunction::F( drs, h );
         const VectorType gradW = r_ij * F;

         const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, drs, diffusiveTermsParams );
         const RealType diffTerm =  psi * ( r_ij, gradW ) * m / rho_j;
         *drho_i += ( v_ij, gradW ) * m - diffTerm;
      }
   };

   auto particleLoopBoundary = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points_bound[ i ];
      const VectorType v_i = view_v_bound[ i ];
      const RealType rho_i = view_rho_bound[ i ];
      const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );
      RealType drho_i = 0.f;

      Particles::NeighborsLoopAnotherSet::exec( i, r_i, searchInOpenBoundary, BoundFluid, v_i, rho_i, p_i, &drho_i );
      view_Drho_bound[ i ] += drho_i;
   };
   boundary->particles->forAll( particleLoopBoundary );

}

template< typename Particles, typename ModelConfig >
template< typename FluidPointer,
          typename BoundaryPointer,
          typename BCType,
          typename std::enable_if_t< std::is_same_v< BCType, WCSPH_BCTypes::DBC >, bool > Enabled >

void
WCSPH_DBC< Particles, ModelConfig >::finalizeBoundaryInteraction( FluidPointer& fluid,
                                                                  BoundaryPointer& boundary,
                                                                  ModelParams& modelParams )
{

}

} // SPH
} // TNL

