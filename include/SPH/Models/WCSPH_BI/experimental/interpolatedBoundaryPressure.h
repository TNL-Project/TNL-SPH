#include "../Interactions.h"

namespace TNL {
namespace SPH {

template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename BoudaryPointer >
requires std::is_same_v< typename ModelConfig::BCType, WCSPH_BCTypes::BIConsistent_numeric_interpolated >
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
   const RealType dynamicViscosity = modelParams.dynamicViscosity;
   const VectorType gravity = modelParams.gravity;
   const RealType preventZero = modelParams.eps * h * h;
   const RealType mgvt_const = ( ParticlesType::spaceDimension == 2 ) ? 8.f : 10.f;
   typename ModelParams::EOS::ParamsType eosParams( modelParams );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->getParticles()->getPoints().getView();
   const auto view_rho = fluid->getVariables()->rho.getView();
   const auto view_v = boundary->getVariables()->rho.getView();

   const auto view_points_bound = boundary->getParticles()->getPoints().getView();
   auto view_rho_bound = boundary->getVariables()->rho.getView();
   auto view_v_bound = boundary->getVariables()->v.getView();
   auto view_dvdt_bound = boundary->getVariables()->a.getView();
   auto view_gamma_bound = boundary->getVariables()->gamma.getView();

   auto getVelocityLaplacian = [ = ] __cuda_callable__( LocalIndexType i,
                                                        LocalIndexType j,
                                                        VectorType& r_i,
                                                        VectorType& v_i,
                                                        VectorType* lap_v_i,
                                                        RealType* gamma_i ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius ) {
         const RealType rho_j = view_rho[ j ];
         const VectorType v_j = view_v[ j ];
         const VectorType v_ij = v_i - v_j;
         const RealType WV_j = KernelFunction::W( drs, h ) * m / rho_j;

         // MGVT:
         //const VectorType gradWV_j = r_ij * KernelFunction::F( drs, h ) * m / rho_j;
         //*lap_v_i -= mgvt_const * ( r_ij, v_ij ) / ( drs * drs + preventZero ) * gradWV_j;

         // Morris - free slip:
         const RealType FV_j = KernelFunction::F( drs, h ) * m / rho_j;
         *lap_v_i -=  2.f * FV_j * v_ij;

         *gamma_i += WV_j;
      }
   };

   auto interpolate = [ = ] __cuda_callable__( LocalIndexType i,
                                               LocalIndexType j,
                                               VectorType& r_i,
                                               VectorType& grad_p_rho_i,
                                               RealType* p_i ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius ) {
         const RealType rho_j = view_rho[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );
         const RealType WV_j = KernelFunction::W( drs, h ) * m / rho_j;;

         *p_i += ( p_j + ( grad_p_rho_i, r_ij ) ) * WV_j;
      }

   };

   auto particleLoopBoundaryConservative = [ = ] __cuda_callable__( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points_bound[ i ];
      const VectorType v_i = view_v_bound[ i ];
      const VectorType dvdt_i = view_dvdt_bound[ i ];
      RealType rho_i = view_rho_bound[ i ];
      RealType gamma_i = 0.f;
      VectorType lap_v_i = 0.f;
      RealType p_i = 0.f;

      Particles::NeighborsLoopAnotherSet::exec( i, r_i, searchInFluid, getVelocityLaplacian, v_i, &lap_v_i, &gamma_i );
      view_gamma_bound[ i ] = gamma_i;

      // interpolation shepard
      if( gamma_i < 1e-3f )
         gamma_i = 1.f;

      const VectorType grad_p_i = dynamicViscosity / rho0 * lap_v_i / gamma_i + gravity - dvdt_i;
      const VectorType grad_p_rho_i = rho_i * grad_p_i;
      Particles::NeighborsLoopAnotherSet::exec( i, r_i, searchInFluid, interpolate, grad_p_rho_i, &p_i );

      p_i /= gamma_i;
      rho_i = EOS::pressureToDensity( p_i, eosParams );
      view_rho_bound[ i ] = rho_i;
   };
   boundary->getParticles()->forAll( particleLoopBoundaryConservative );
}

template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename BoundaryPointer >
requires std::is_same_v< typename ModelConfig::BCType, WCSPH_BCTypes::BIConsistent_numeric_interpolated >
void
WCSPH_BI< Particles, ModelConfig >::finalizeBoundaryInteraction( FluidPointer& fluid,
                                                                 BoundaryPointer& boundary,
                                                                 ModelParams& modelParams )
{}

}  //namespace SPH
}  //namespace TNL

