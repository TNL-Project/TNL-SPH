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
   auto searchInFluid = fluid->getParticles()->getSearchToken( fluid->particles );
   auto searchInBound = fluid->getParticles()->getSearchToken( boundary->particles );

   Particles* _fluidParticles = &fluid->particles.template getData< DeviceType >();
   Particles* _boundParticles = &fluid->particles.template getData< DeviceType >();
   Variables* _fluidVars = &fluid->variables.template modifyData< DeviceType >();
   Variables* _boundVars = &boundary->variables.template modifyData< DeviceType >();
   ModelParams* _modelParams = &modelParams;

   auto fluidFluid = [ _fluidParticles, _fluidVars ] __cuda_callable__
      ( LocalIndexType i, LocalIndexType j, VectorType& r_i, InteractionState & interastionState ) mutable
   {
      SPHModel::fluidFluidInteraction( *_particles, *_variables, *_interactionState );
   };

   auto FluidBound = [ _boundParticles, _boundVars ] __cuda_callable__
      ( LocalIndexType i, LocalIndexType j, VectorType& r_i, InteractionState & interactionState ) mutable
   {
      SPHModel::fluidBoundaryInteraction( *_boundParticles, *_boundVars, *_interactionState );
   };

   auto particleLoop = [ _particles, _vars, _modelParams ] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      SPHModel::InteractionState interactionState( *_particles, *_variables, i, modelParams );
      Particles::NeighborsLoop::exec( i, interactionState.r_i, searchInFluid, fluidFluid, interactionState );
      Particles::NeighborsLoopAnotherSet::exec( i, interactionState.r_i, searchInBound, fluidBound, interactionState );
      interactionState.writeFluid( *_particles, *_variables, i );
   };
   fluid->particles->forAll( particleLoop );
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

      Particles::NeighborsLoopAnotherSet::exec(
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
{}

template< typename Particles, typename ModelConfig >
template< typename FluidPointer,
          typename BoundaryPointer,
          typename BCType,
          typename std::enable_if_t< std::is_same_v< BCType, WCSPH_BCTypes::MDBC >, bool > Enabled >

void
WCSPH_DBC< Particles, ModelConfig >::finalizeInteraction( FluidPointer& fluid,
                                                          BoundaryPointer& boundary,
                                                          ModelParams& modelParams )
{}

} // SPH
} // TNL
