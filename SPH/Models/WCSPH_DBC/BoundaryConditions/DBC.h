#pragma once

#include "../Variables.h"
#include "../../../SPHTraits.h"
#include "../../../../Particles/neighborSearchLoop.h"
#include "../../../customParallelFor.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {
namespace BoundaryConditions {

template< typename SPHConfig >
class DBCVariables : public SPHFluidVariables< SPHConfig >
{
   public:
   using BaseType = SPHFluidVariables< SPHConfig >;
   using GlobalIndexType = typename BaseType::GlobalIndexType;

   DBCVariables( GlobalIndexType size )
   : SPHFluidVariables< SPHConfig >( size ) {};
};

template< typename SPHConfig >
class DBC
{
   public:
   using Variables = DBCVariables< SPHConfig >;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using DeviceType = typename SPHConfig::DeviceType;

   using LocalIndexType = typename SPHTraitsType::LocalIndexType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using ScalarType = typename SPHTraitsType::ScalarType;
   using VectorType = typename SPHTraitsType::VectorType;
   using IndexVectorType = typename SPHTraitsType::IndexVectorType;

   template< typename ParticleSystem,
             typename FluidPointer,
             typename BoudaryPointer,
             typename SPHKernelFunction,
             typename DiffusiveTerm,
             typename EOS,
             typename SPHState >
   static void
   update( FluidPointer& fluid, BoudaryPointer& boundary, SPHState& sphState )
   {
      /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
      typename ParticleSystem::NeighborsLoopParams searchInFluid( fluid->particles );
      typename ParticleSystem::NeighborsLoopParams searchInBound( boundary->particles );

      /* CONSTANT VARIABLES */
      const RealType searchRadius = fluid->particles->getSearchRadius();
      const RealType h = sphState.h;
      const RealType m = sphState.mass;

      typename DiffusiveTerm::ParamsType diffusiveTermsParams( sphState );
      typename EOS::ParamsType eosParams( sphState );

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

            const RealType F = SPHKernelFunction::F( drs, h );
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

         NeighborsLoop::exec( i, r_i, searchInFluid, BoundFluid, v_i, rho_i, p_i, &drho_i );

         view_Drho_bound[ i ] = drho_i;
      };
      SPHParallelFor::exec( boundary->getFirstActiveParticle(), boundary->getLastActiveParticle() + 1, particleLoopBoundary );
   }

};

} // BoundaryConditions
} // SPH
} // ParticleSystem
} // TNL

