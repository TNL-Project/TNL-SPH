#include "Interactions.h"
#include <execution>
#include "details.h"

namespace TNL {
namespace SPH {

template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename BoudaryPointer >
void
WCSPH_DBC< Particles, ModelConfig >::interaction( FluidPointer& fluid,
                                                  BoudaryPointer& boundary,
                                                  ModelParams& modelParams )
{
   // searchable objects
   typename Particles::NeighborsLoopParams searchInFluid( fluid->particles );
   typename Particles::NeighborsLoopParams searchInBound( boundary->particles );

   // load constant variables
   const RealType searchRadius = fluid->particles->getSearchRadius();
   const RealType h = modelParams.h;
   const RealType m = modelParams.mass;
   const VectorType gravity = modelParams.gravity;
   typename DiffusiveTerm::ParamsType diffusiveTermsParams( modelParams );
   typename ViscousTerm::ParamsType viscousTermTermsParams( modelParams );
   typename EOS::ParamsType eosParams( modelParams );

   // load variables
   const auto view_points = fluid->particles->getPoints().getConstView();
   const auto view_rho = fluid->variables->rho.getConstView();
   auto view_Drho = fluid->variables->drho.getView();
   const auto view_v = fluid->variables->v.getConstView();
   auto view_a = fluid->variables->a.getView();

   const auto view_points_bound = boundary->particles->getPoints().getConstView();
   const auto view_rho_bound = boundary->variables->rho.getConstView();
   const auto view_v_bound = boundary->variables->v.getConstView();

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

         /* Interaction: */
         const VectorType v_ij = v_i - v_j;

         const RealType F = KernelFunction::F( drs, h );
         const VectorType gradW = r_ij * F;

         const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, drs, diffusiveTermsParams );
         const RealType diffTerm =  psi * ( r_ij, gradW ) * m / rho_j;
         *drho_i += ( v_ij, gradW ) * m - diffTerm;

         const RealType p_term = ( p_i + p_j ) / ( rho_i * rho_j );
         const RealType visco =  ViscousTerm::Pi( rho_i, rho_j, drs, ( r_ij, v_ij ), viscousTermTermsParams );
         *a_i += ( -1.0f ) * ( p_term + visco ) * gradW * m;
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

         /* Interaction: */
         const VectorType v_ij = v_i - v_j;

         const RealType F = KernelFunction::F( drs, h );
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

      TNL::ParticleSystem::NeighborsLoop::exec( i, r_i, searchInFluid, FluidFluid, v_i, rho_i, p_i, &drho_i, &a_i );
      TNL::ParticleSystem::NeighborsLoopAnotherSet::exec( i, r_i, searchInBound, FluidBound, v_i, rho_i, p_i, &drho_i, &a_i );

      view_Drho[ i ] = drho_i;
      a_i += gravity;
      view_a[ i ] = a_i;
   };
   fluid->particles->forAll( particleLoop );

   if constexpr( Model::ModelConfigType::SPHConfig::numberOfPeriodicBuffers > 0 ){
      for( long unsigned int i = 0; i < std::size( fluid->periodicPatches ); i++ ){

         const auto zoneParticleIndices_view = fluid->periodicPatches[ i ]->particleZone.getParticlesInZone().getConstView();
         const GlobalIndexType numberOfZoneParticles = fluid->periodicPatches[ i ]->particleZone.getNumberOfParticles();
         const VectorType shift = fluid->periodicPatches[ i ]->config.shift;

         auto periodicParticleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
         {
            const GlobalIndexType p = zoneParticleIndices_view[ i ];
            const VectorType r_i = view_points[ p ] + shift;
            const VectorType v_i = view_v[ p ];
            const RealType rho_i = view_rho[ p ];
            const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );
            VectorType a_i = 0.f;
            RealType drho_i = 0.f;

            TNL::ParticleSystem::NeighborsLoop::exec(
                  p, r_i, searchInFluid, FluidFluid, v_i, rho_i, p_i, &drho_i, &a_i );
            TNL::ParticleSystem::NeighborsLoopAnotherSet::exec(
                  p, r_i, searchInBound, FluidBound, v_i, rho_i, p_i, &drho_i, &a_i );

            view_Drho[ p ] += drho_i;
            view_a[ p ] += a_i;
         };
         Algorithms::parallelFor< DeviceType >( 0, numberOfZoneParticles, periodicParticleLoop );
      }
   }
}

template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename OpenBoudaryPointer >
void
WCSPH_DBC< Particles, ModelConfig >::interactionWithOpenBoundary( FluidPointer& fluid,
                                                                  OpenBoudaryPointer& openBoundary,
                                                                  ModelParams& modelParams )
{
   // searchable object
   typename Particles::NeighborsLoopParams searchInOpenBoundary( openBoundary->particles );

   // load constant variables
   const RealType searchRadius = fluid->particles->getSearchRadius();
   const RealType h = modelParams.h;
   const RealType m = modelParams.mass;
   typename DiffusiveTerm::ParamsType diffusiveTermsParams( modelParams );
   typename ViscousTerm::ParamsType viscousTermTermsParams( modelParams );
   typename EOS::ParamsType eosParams( modelParams );

   // load variables
   const auto view_points = fluid->particles->getPoints().getView();
   const auto view_rho = fluid->variables->rho.getView();
   auto view_Drho = fluid->variables->drho.getView();
   const auto view_v = fluid->variables->v.getView();
   auto view_a = fluid->variables->a.getView();

   const auto view_points_openBound = openBoundary->particles->getPoints().getView();
   auto view_rho_openBound = openBoundary->variables->rho.getView();
   auto view_v_openBound = openBoundary->variables->v.getView();

   const auto zoneParticleIndices_view = openBoundary->zone.getParticlesInZone().getConstView();
   const GlobalIndexType numberOfZoneParticles = openBoundary->zone.getNumberOfParticles();

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

         const RealType F = KernelFunction::F( drs, h );
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
      const GlobalIndexType p = zoneParticleIndices_view[ i ];
      const VectorType r_i = view_points[ p ];
      const VectorType v_i = view_v[ p ];
      const RealType rho_i = view_rho[ p ];
      const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );
      VectorType a_i = 0.f;
      RealType drho_i = 0.f;

      TNL::ParticleSystem::NeighborsLoopAnotherSet::exec(
            p, r_i, searchInOpenBoundary, FluidOpenBoundary, v_i, rho_i, p_i, &drho_i, &a_i );

      view_Drho[ p ] += drho_i;
      view_a[ p ] += a_i;
   };
   Algorithms::parallelFor< DeviceType >( 0, numberOfZoneParticles, particleLoop );
}

template< typename Particles, typename ModelConfig >
template< typename EquationOfState, typename PhysicalObjectPointer >
void
WCSPH_DBC< Particles, ModelConfig >::computePressureFromDensity( PhysicalObjectPointer& physicalObject,
                                                                 ModelParams& modelParams )
{
   auto view_rho = physicalObject->getVariables()->rho.getView();
   auto view_p = physicalObject->getVariables()->p.getView();
   typename EquationOfState::ParamsType eosParams( modelParams );

   auto evalPressure = [=] __cuda_callable__ ( int i ) mutable
   {
      view_p[ i ] = EquationOfState::DensityToPressure( view_rho[ i ], eosParams );
   };
   physicalObject->particles->forAll( evalPressure ); //TODO: forloop?
}

template< typename Particles, typename ModelConfig >
template< typename FluidPointer,
          typename BoundaryPointer,
          typename BCType,
          typename std::enable_if_t< std::is_same_v< BCType, WCSPH_BCTypes::DBC >, bool > Enabled >

void
WCSPH_DBC< Particles, ModelConfig >::finalizeInteraction( FluidPointer& fluid,
                                                          BoundaryPointer& boundary,
                                                          ModelParams& modelParams )
{

}

template< typename Particles, typename ModelConfig >
template< typename FluidPointer,
          typename BoundaryPointer,
          typename BCType,
          typename std::enable_if_t< std::is_same_v< BCType, WCSPH_BCTypes::MDBC >, bool > Enabled >

void
WCSPH_DBC< Particles, ModelConfig >::finalizeInteraction( FluidPointer& fluid,
                                                          BoundaryPointer& boundary,
                                                          ModelParams& modelParams )
{
   const RealType rho0 = modelParams.rho0;

   auto view_rho_bound = boundary->variables->rho.getView();
   const auto view_points_bound = boundary->particles->getPoints().getConstView();
   const auto view_ghostNode_bound = boundary->variables->ghostNodes.getConstView();
   const auto view_rhoGradRhoGhostNode_bound = boundary->variables->rhoGradRho_gn.getConstView();
   const auto view_correctionMatrices_bound = boundary->variables->cMatrix_gn.getConstView();

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points_bound[ i ];
      const VectorType ghostNode_i = view_ghostNode_bound[ i ];
      const Matrix cMatrix_gn = view_correctionMatrices_bound[ i ];
      const VectorExtendedType rhoGradRho_gn = view_rhoGradRhoGhostNode_bound[ i ];
      RealType rho_bound = 0.f;

      if( Matrices::determinant( cMatrix_gn ) > 0.001 ) {
         VectorExtendedType cRhoGradRho = Matrices::solve( cMatrix_gn, rhoGradRho_gn );
         VectorType r_ign = ghostNode_i - r_i;
         rho_bound = cRhoGradRho[ 0 ] + cRhoGradRho[ 1 ] * r_ign[ 0 ] + cRhoGradRho[ 2 ] * r_ign[ 1 ];
      }
      else if( cMatrix_gn( 0, 0 ) > 0.f ) {
         rho_bound = rhoGradRho_gn[ 0 ] / cMatrix_gn( 0, 0 );
      }
      else {
         rho_bound = rho0;
      }

      if( rho_bound < rho0 )
         rho_bound = rho0;

      view_rho_bound[ i ] = rho_bound;
   };
   boundary->particles->forAll( particleLoop ); //TODO: forloop?
}

} // SPH
} // TNL
