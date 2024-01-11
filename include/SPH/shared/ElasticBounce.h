#include <TNL/Algorithms/parallelFor.h>

#include "../SPHTraits.h"
#include "../../Particles/neighborSearchLoop.h"

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
                    SPHState& sphState )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   GlobalIndexType numberOfParticles = fluid->particles->getNumberOfParticles();
   const RealType searchRadius = fluid->particles->getSearchRadius();

   typename ParticleSystem::NeighborsLoopParams searchInBound( boundary->particles );

   /* CONSTANT VARIABLES */
   const RealType r_box = sphState.r_box;
   const RealType minimalDistanceFactor = sphState.minimalDistanceFactor;
   const RealType elasticFactor = sphState.elasticFactor;
   const RealType dt = sphState.dtInit; //TODO: Use dynamic timestep
   const RealType dp = sphState.dp;

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
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if (drs <= searchRadius )
      {
         const VectorType v_j = view_v_bound[ j ];
         const VectorType n_j = view_n_bound[ j ];

         const RealType r0 = ( r_ij, n_j );

         if( r0 < 0.f ) //Particle is behind the wall
            return;

         if( ( r_ij, r_ij ) >= r_box * r_box ) //Particle is to far from boundary element
            return;

         const VectorType v_ij = *ve_i - v_j;
         const RealType v_n = ( -1.f ) * ( v_ij, n_j );
         const RealType dvdt_n = ( -1.f ) * ( *ae_i, n_j ); //at this point, we assume boundary with no dvdt_j

         const RealType r_n = dt * v_n + 0.5f * dt * dt *  dvdt_n ;

         if( r_n < 0.f ) //Particle is already running away from boundary
            return;

         if( r0 - r_n <= minimalDistanceFactor * dp )
         {
            *ve_i = ( *ve_i ) - ( 1.f + elasticFactor ) * ( -1.f ) * v_n * n_j; // (-1.f) due to inner normal
            *ae_i = ( *ae_i ) - ( 1.f + elasticFactor ) * ( -1.f ) * dvdt_n * n_j;
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

