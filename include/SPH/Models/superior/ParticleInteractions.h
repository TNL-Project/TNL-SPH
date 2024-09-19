#include "../../SPHTraits.h"

namespace TNL {
namespace SPH {

template< typename ModelConfig >
class DBCModelInteractionState
{
public:
   using SPHTraitsType = SPHFluidTraits< typename ModelConfig::SPHBaseConfig >;
   using IndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   // exchangable model terms
   using KernelFunction = typename ModelConfig::KernelFunction;
   using DiffusiveTerm = typename ModelConfig::DiffusiveTerm;
   using ViscousTerm = typename ModelConfig::ViscousTerm;
   using EOS = typename ModelConfig::EOS;

   template< typename Particles, typename Variabels, typename ModelParams >
   DBCModelInteractionState( const Particles & particles,
                             const Variabels & vars,
                             const ModelParams & modelParams,
                             const IndexType i )
   {
      searchRadius = particles.getSearchRadius();
      m = modelParams.m;
      h = modelParams.h;

      diffusiveTermParams.init( modelParams );
      viscousTermTermParams.init( modelParams );
      eosParams.init( modelParams );

      r_i = particles.getPoint( i );
      v_i = vars.v[ i ];
      rho_i = vars.rho[ i ];
      p_i = EOS::DensityToPressure( rho_i, eosParams );
   }

   template< typename Particles, typename Variabels >
   void
   write( const Particles & particles, Variabels & vars, const IndexType & i )
   {
      vars.drhodt_i[ i ] = drhodt_i;
      dvdt_i += gravity;
      vars.dvdt[ i ] = dvdt_i;
   }

   // constants
   RealType searchRadius;
   RealType m;
   RealType h;
   VectorType gravity;

   // exchangable model terms constants
   typename DiffusiveTerm::ParamsType diffusiveTermParams;
   typename ViscousTerm::ParamsType viscousTermTermParams;
   typename EOS::ParamsType eosParams;

   // particle i loaded variables
   VectorType r_i;
   VectorType v_i;
   RealType rho_i;
   RealType p_i;

   // variables summed over neighbors j
   VectorType dvdt_i = 0.f;
   RealType drhodt_i = 0.f;
};

template< typename ModelConfig >
class DBCModelInteractions
{
public:
   using SPHTraitsType = SPHFluidTraits< typename ModelConfig::SPHBaseConfig >;
   using IndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   using InteractionState =  DBCModelInteractionState< ModelConfig >;

   // exchangable model terms
   using KernelFunction = typename ModelConfig::KernelFunction;
   using DiffusiveTerm = typename ModelConfig::DiffusiveTerm;
   using ViscousTerm = typename ModelConfig::ViscousTerm;
   using EOS = typename ModelConfig::EOS;

   template< typename Particles, typename Variables >
   __cuda_callable__
   static void
   fluidFluidInteraction( const Particles & particles,
                          const Variables & vars,
                          const IndexType & j,
                          InteractionState & state )
   {
      const VectorType r_j = particles.getPoint( j );
      const VectorType r_ij = state.r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if (drs <= state.searchRadius )
      {
         const VectorType v_j = vars.v[ j ];
         const RealType rho_j = vars.rho[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, state.eosParams );

         const VectorType v_ij = state.v_i - v_j;

         const RealType F = KernelFunction::F( drs, state.h );
         const VectorType gradW = r_ij * F;

         const RealType psi = DiffusiveTerm::Psi( state.rho_i, rho_j, drs, state.diffusiveTermsParams );
         const RealType diffTerm =  psi * ( r_ij, gradW ) * state.m / rho_j;
         state.drhodt_i += ( v_ij, gradW ) * state.m - diffTerm;

         const RealType p_term = ( state.p_i + p_j ) / ( state.rho_i * rho_j );
         const RealType visco =  ViscousTerm::Pi( state.rho_i, rho_j, drs, ( r_ij, v_ij ), state.viscousTermTermsParams );
         state.dvdt_i += ( -1.0f ) * ( p_term + visco ) * gradW * state.m;
      }
   }

   template< typename Particles, typename Variables >
   __cuda_callable__
   static void
   fluidBoundaryInteraction( const Particles & boundParticles,
                             const Variables & boundVars,
                             const IndexType & j,
                             InteractionState & state )
   {
      const VectorType r_j = boundParticles.getPoint( j );
      const VectorType r_ij = state.r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if (drs <= state.searchRadius )
      {
         const VectorType v_j = boundVars.v[ j ];
         const RealType rho_j = boundVars.rho[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, state.eosParams );

         const VectorType v_ij = state.v_i - v_j;

         const RealType F = KernelFunction::F( drs, state.h );
         const VectorType gradW = r_ij * F;

         const RealType psi = DiffusiveTerm::Psi( state.rho_i, rho_j, drs, state.diffusiveTermsParams );
         const RealType diffTerm =  psi * ( r_ij, gradW ) * state.m / rho_j;
         state.drhodt_i += ( v_ij, gradW ) * state.m - diffTerm;

         const RealType p_term = ( state.p_i + p_j ) / ( state.rho_i * rho_j );
         const RealType visco =  ViscousTerm::Pi( state.rho_i, rho_j, drs, ( r_ij, v_ij ), state.viscousTermTermsParams );
         state.dvdt_i += ( -1.0f ) * ( p_term + visco ) * gradW * state.m;
      }
    }
};

}
}

