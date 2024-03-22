#include "../Interactions.h"
#include "../details.h"

namespace TNL {
namespace SPH {

template< typename Particles, typename ModelConfig >
template< typename FluidPointer,
          typename BoudaryPointer,
          typename BCType,
          typename std::enable_if_t< std::is_same_v< BCType, WCSPH_BCTypes::MDBC >, bool > Enabled >
void
WCSPH_DBC< Particles, ModelConfig >::updateSolidBoundary( FluidPointer& fluid,
                                                          BoudaryPointer& boundary,
                                                          ModelParams& modelParams )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   typename Particles::NeighborsLoopParams searchInFluid( fluid->particles );
   typename Particles::NeighborsLoopParams searchInBound( boundary->particles );

   /* CONSTANT VARIABLES */
   const RealType searchRadius = fluid->particles->getSearchRadius();
   const RealType h = modelParams.h;
   const RealType m = modelParams.mass;
   const RealType rho0 = modelParams.rho0;

   typename EOS::ParamsType eosParams( modelParams );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->particles->getPoints().getView();
   const auto view_rho = fluid->variables->rho.getView();

   const auto view_points_bound = boundary->particles->getPoints().getView();
   auto view_rho_bound = boundary->variables->rho.getView();
   const auto view_v_bound = boundary->variables->v.getView();
   const auto view_ghostNode_bound = boundary->variables->ghostNodes.getView();
   auto view_rhoGradRhoGhostNode_bound = boundary->variables->rhoGradRho_gn.getView();
   auto view_correctionMatrices_bound = boundary->variables->cMatrix_gn.getView();

   auto BoundFluid = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j,
         VectorType& ghostNode_i, VectorType& v_i, RealType& rho_i, RealType& p_i, Matrix* A_gn, VectorExtendedType* b_gn ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_j - ghostNode_i; //FLUID_POS - GHOSTNODE_POS
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         const RealType rho_j = view_rho[ j ];

         /* Interaction */
         const RealType F = KernelFunction::F( drs, h );
         const RealType W = KernelFunction::W( drs, h );
         const VectorType gradW = r_ij * F;

         const RealType V = m / rho_j;

         Matrix A_local = 0.f;
         VectorExtendedType b_local = 0.f;

         //TODO: We could bamybe handle this automatically by overloading the matrixCorrection and gVVG function
         if( SPHConfig::spaceDimension == 2 )
         {
            *A_gn += matrixCorrection2D< Matrix >( W, gradW, r_ij, V );
            *b_gn += getVariableValueAndGradient2D< VectorExtendedType >( W, gradW, rho_j, V );
         }
      }
   };

   auto particleLoopBoundary = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points_bound[ i ];
      const VectorType v_i = view_v_bound[ i ];
      const RealType rho_i = view_rho_bound[ i ];
      const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );
      const VectorType ghostNode_i = view_ghostNode_bound[ i ];

      RealType drho_i = 0.f;
      Matrix A_gn = 0.f;
      VectorExtendedType b_gn = 0.f;

      TNL::ParticleSystem::NeighborsLoop::exec( i, ghostNode_i, searchInFluid, BoundFluid, v_i, rho_i, p_i, &A_gn, &b_gn );

      view_rhoGradRhoGhostNode_bound[ i ] = b_gn;
      view_correctionMatrices_bound[ i ] = A_gn;
   };
   TNL::Algorithms::parallelFor< DeviceType >(
         boundary->getFirstActiveParticle(), boundary->getLastActiveParticle() + 1, particleLoopBoundary );

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
            const VectorType ghostNode_i = view_ghostNode_bound[ p ] + shift;

            RealType drho_i = 0.f;
            Matrix A_gn = 0.f;
            VectorExtendedType b_gn = 0.f;

            TNL::ParticleSystem::NeighborsLoop::exec( p, ghostNode_i, searchInFluid, BoundFluid, v_i, rho_i, p_i, &A_gn, &b_gn );

            view_rhoGradRhoGhostNode_bound[ p ] += b_gn;
            view_correctionMatrices_bound[ p ] += A_gn;
         };
         Algorithms::parallelFor< DeviceType >( 0, numberOfZoneParticles, periodicParticleLoopBoundary );
      }
   }
}

template< typename Particles, typename ModelConfig >
template< typename OpenBoundaryPointer,
          typename BoudaryPointer,
          typename BCType,
          typename std::enable_if_t< std::is_same_v< BCType, WCSPH_BCTypes::MDBC >, bool > Enabled >
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
   const RealType rho0 = modelParams.rho0;

   typename EOS::ParamsType eosParams( modelParams );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points_openBound = openBoundary->particles->getPoints().getView();
   const auto view_rho_openBound = openBoundary->variables->rho.getView();

   const auto view_points_bound = boundary->particles->getPoints().getView();
   auto view_rho_bound = boundary->variables->rho.getView();
   const auto view_v_bound = boundary->variables->v.getView();
   const auto view_ghostNode_bound = boundary->variables->ghostNodes.getView();
   auto view_rhoGradRhoGhostNode_bound = boundary->variables->rhoGradRho_gn.getView();
   auto view_correctionMatrices_bound = boundary->variables->cMatrix_gn.getView();

   auto BoundOpenBoundary = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j,
         VectorType& ghostNode_i, VectorType& v_i, RealType& rho_i, RealType& p_i, Matrix* A_gn, VectorExtendedType* b_gn ) mutable
   {
      const VectorType r_j = view_points_openBound[ j ];
      const VectorType r_ij = r_j - ghostNode_i; //FLUID_POS - GHOSTNODE_POS
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         const RealType rho_j = view_rho_openBound[ j ];

         /* Interaction */
         const RealType F = KernelFunction::F( drs, h );
         const RealType W = KernelFunction::W( drs, h );
         const VectorType gradW = r_ij * F;

         const RealType V = m / rho_j;

         Matrix A_local = 0.f;
         VectorExtendedType b_local = 0.f;

         //TODO: We could bamybe handle this automatically by overloading the matrixCorrection and gVVG function
         if( SPHConfig::spaceDimension == 2 )
         {
            *A_gn += matrixCorrection2D< Matrix >( W, gradW, r_ij, V );
            *b_gn += getVariableValueAndGradient2D< VectorExtendedType >( W, gradW, rho_j, V );
         }
      }
   };

   auto particleLoopBoundary = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points_bound[ i ];
      const VectorType v_i = view_v_bound[ i ];
      const RealType rho_i = view_rho_bound[ i ];
      const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );
      const VectorType ghostNode_i = view_ghostNode_bound[ i ];

      RealType drho_i = 0.f;
      Matrix A_gn = 0.f;
      VectorExtendedType b_gn = 0.f;

      TNL::ParticleSystem::NeighborsLoop::exec(
            i, ghostNode_i, searchInOpenBoundary, BoundOpenBoundary, v_i, rho_i, p_i, &A_gn, &b_gn );

      view_rhoGradRhoGhostNode_bound[ i ] += b_gn;
      view_correctionMatrices_bound[ i ] += A_gn;
   };
   TNL::Algorithms::parallelFor< DeviceType >(
         boundary->getFirstActiveParticle(), boundary->getLastActiveParticle() + 1, particleLoopBoundary );
}

} // SPH
} // TNL

