#pragma once
#include <TNL/Algorithms/parallelFor.h>

#include "../SPHTraits.h"

namespace TNL {
namespace SPH {

template< typename ParticleSystem, typename SPHConfig>
class ElasticBounce
{
public:

   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using DeviceType = typename SPHConfig::DeviceType;

   using LocalIndexType = typename SPHTraitsType::LocalIndexType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

template< typename FluidPointer, typename BoudaryPointer, typename SPHState >
static void
boundaryCorrection( FluidPointer& fluid,
                    BoudaryPointer& boundary,
                    SPHState& sphState,
                    const RealType& dt )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   GlobalIndexType numberOfParticles = fluid->particles->getNumberOfParticles();
   const RealType searchRadius = fluid->particles->getSearchRadius();

   typename ParticleSystem::NeighborsLoopParams searchInBound( boundary->particles );

   /* CONSTANT VARIABLES */
   const RealType dp = sphState.dp;
   const RealType r_box = sphState.r_boxFactor * dp;
   const RealType minimalDistanceFactor = sphState.minimalDistanceFactor;
   const RealType elasticFactor = sphState.elasticFactor;

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->particles->getPoints().getView();
   auto view_v = fluid->variables->v.getView();
   auto view_a = fluid->variables->a.getView();

   const auto view_points_bound = boundary->particles->getPoints().getView();
   const auto view_v_bound = boundary->variables->v.getView();
   const auto view_n_bound = boundary->variables->n.getView();

   auto elasticBounce = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j,
         VectorType& r_i, VectorType* ve_i, VectorType* ae_i ) mutable
   {
      const VectorType r_j = view_points_bound[ j ];
      const VectorType r_ji = r_j - r_i;
      const RealType drs = l2Norm( r_ji );
      if (drs <= searchRadius )
      {
         const VectorType v_j = view_v_bound[ j ];
         const VectorType n_j = ( -1.f ) * view_n_bound[ j ];

         const RealType r0 = ( r_ji, n_j );
         // Particle is behind the wall
         if( r0 < 0.f )
            return;

         const VectorType r_ji_box = r_ji - r0 * n_j;
         // Particle is too far from boundary element
         if( ( r_ji_box, r_ji_box ) >= r_box * r_box )
            return;

         const RealType v_n = ( *ve_i - v_j, n_j );
         const RealType dvdt_n = ( *ae_i, n_j ); //TODO: include dvdt_j
         const RealType r_n = dt * v_n + 0.5f * dt * dt *  dvdt_n ;

         // Particle is already running away from boundary
         if( r_n < 0.f )
            return;

         // Reflect the particle
         if( r0 - r_n <= minimalDistanceFactor * dp ){
            *ve_i = ( *ve_i ) - ( 1.f + elasticFactor ) * v_n * n_j;
            *ae_i = ( *ae_i ) - ( 1.f + elasticFactor ) * dvdt_n * n_j;
         }
      }
   };

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points[ i ];
      const VectorType v_i = view_v[ i ];
      const VectorType dvdt_i = view_a[ i ];

      VectorType ve_i = v_i;
      VectorType ae_i = dvdt_i;

      ParticleSystem::NeighborsLoopAnotherSet::exec( i, r_i, searchInBound, elasticBounce, &ve_i, &ae_i );

      view_v[ i ] = ve_i;
      view_a[ i ] = ae_i;

   };
   Algorithms::parallelFor< DeviceType >( 0, numberOfParticles, particleLoop );

}
};


} // SPH
} // TNL

