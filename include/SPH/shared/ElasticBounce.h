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
      GlobalIndexType numberOfParticles = fluid->getParticles()->getNumberOfParticles();
      const RealType searchRadius = fluid->getParticles()->getSearchRadius();

      typename ParticleSystem::NeighborsLoopParams searchInBound( boundary->getParticles() );

      /* CONSTANT VARIABLES */
      const RealType dp = sphState.dp;
      const RealType r_box = sphState.r_boxFactor * dp;
      const RealType minimalDistanceFactor = sphState.minimalDistanceFactor;
      const RealType elasticFactor = sphState.elasticFactor;

      /* VARIABLES AND FIELD ARRAYS */
      const auto view_points = fluid->getParticles()->getPoints().getView();
      auto view_v = fluid->getVariables()->v.getView();
      auto view_a = fluid->getVariables()->a.getView();

      const auto view_points_bound = boundary->getParticles()->getPoints().getView();
      const auto view_v_bound = boundary->getVariables()->v.getView();
      const auto view_n_bound = boundary->getVariables()->n.getView();

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

   template< typename FluidPointer, typename BoudaryPointer, typename SPHState >
   static void
   boundaryCorrectionPST( FluidPointer& fluid,
                          BoudaryPointer& boundary,
                          SPHState& sphState,
                          const RealType& dt )
   {
      // do nothing
   }

};

template< typename ParticleSystem, typename SPHConfig>
class ElasticBounceLight
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
      GlobalIndexType numberOfParticles = fluid->getParticles()->getNumberOfParticles();
      const RealType searchRadius = fluid->getParticles()->getSearchRadius();

      typename ParticleSystem::NeighborsLoopParams searchInBound( boundary->getParticles() );

      /* CONSTANT VARIABLES */
      const RealType dp = sphState.dp;
      const RealType r_boxFactor = sphState.r_boxFactor;
      const RealType minimalDistanceFactor = sphState.minimalDistanceFactor;
      const RealType elasticFactor = sphState.elasticFactor;

      /* VARIABLES AND FIELD ARRAYS */
      const auto view_points = fluid->getParticles()->getPoints().getView();
      auto view_v = fluid->getVariables()->v.getView();
      auto view_a = fluid->getVariables()->a.getView();

      const auto view_points_bound = boundary->getParticles()->getPoints().getView();
      const auto view_v_bound = boundary->getVariables()->v.getView();
      const auto view_n_bound = boundary->getVariables()->n.getView();
      const auto view_elementSize_bound = boundary->getVariables()->elementSize.getConstView();

      auto elasticBounce = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j,
            VectorType& r_i,  VectorType& v_i, VectorType* ve_i, VectorType* dvdte_i ) mutable
      {
         const VectorType r_j = view_points_bound[ j ];
         const VectorType r_ji = r_j - r_i;
         const RealType drs = l2Norm( r_ji );
         if (drs <= searchRadius )
         {
            const VectorType v_j = view_v_bound[ j ];
            const VectorType n_j = ( -1.f ) * view_n_bound[ j ];
            const RealType s_j = view_elementSize_bound[ j ];

            const RealType rn = ( r_ji, n_j );
            // Particle is behind the wall
            if( rn < 0.f )
               return;

            RealType r_box;
            if( SPHConfig::spaceDimension == 2 )
               r_box = r_boxFactor * s_j;
            else if( SPHConfig::spaceDimension == 3 )
               r_box = r_boxFactor * sqrt( s_j );
            const VectorType r_ji_box = r_ji - rn * n_j;
            // Particle is too far from boundary element
            if( ( r_ji_box, r_ji_box ) >= r_box * r_box )
               return;

            const RealType drn = dt * ( *ve_i, n_j );
            // Particle is already running away from boundary
            if( drn < 0.f )
               return;

            // Reflect the particle if is inside or is entering effective area
            if( rn - drn <= minimalDistanceFactor * dp ){
               const VectorType v_updated = v_i + dt * ( *dvdte_i );
               const VectorType v_updated_r = v_updated - 2.f * ( v_updated, n_j ) * n_j;
               *dvdte_i = ( v_updated_r - v_i ) / dt;
               *ve_i = v_i + 0.5f * dt * ( *dvdte_i );
            }
         }
      };

      auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
      {
         const VectorType r_i = view_points[ i ];
         const VectorType v_i = view_v[ i ];
         const VectorType dvdt_i = view_a[ i ];

         VectorType ve_i = v_i + 0.5f * dt * dvdt_i;
         VectorType dvdte_i = dvdt_i;

         ParticleSystem::NeighborsLoopAnotherSet::exec( i, r_i, searchInBound, elasticBounce, v_i, &ve_i, &dvdte_i );

         view_a[ i ] = dvdte_i;
      };
      Algorithms::parallelFor< DeviceType >( 0, numberOfParticles, particleLoop );

   }

   template< typename FluidPointer, typename BoudaryPointer, typename SPHState >
   static void
   boundaryCorrectionPST( FluidPointer& fluid,
                          BoudaryPointer& boundary,
                          SPHState& sphState,
                          const RealType& dt )
   {
      GlobalIndexType numberOfParticles = fluid->getParticles()->getNumberOfParticles();
      const RealType searchRadius = fluid->getParticles()->getSearchRadius();
      typename ParticleSystem::NeighborsLoopParams searchInBound( boundary->getParticles() );

      const RealType m = sphState.mass;
      const RealType r_boxFactor = sphState.r_boxFactor;

      auto view_points = fluid->getParticles()->getPoints().getView();
      const auto view_rho = fluid->getVariables()->rho.getConstView();
      const auto view_points_bound = boundary->getParticles()->getPoints().getView();
      const auto view_n_bound = boundary->getVariables()->n.getConstView();
      const auto view_elementSize_bound = boundary->getVariables()->elementSize.getConstView();

      auto pstLoop = [=] __cuda_callable__ ( LocalIndexType i,
                                             LocalIndexType j,
                                             VectorType& r_i,
                                             RealType& Ri,
                                             VectorType* r_i_mod ) mutable
      {
         const VectorType r_j = view_points_bound[ j ];
         //const VectorType r_ji = r_j - ( *r_i_mod );
         const VectorType r_ji = r_j - view_points[ i ];
         const RealType drs = l2Norm( r_ji );
         if (drs <= searchRadius )
         {
            const VectorType n_j = ( -1.f ) * view_n_bound[ j ];
            const RealType r_ji_n = ( r_ji, n_j );
            const RealType s_j = view_elementSize_bound[ j ];

            // Particle is too far from boundary or it is righ placed
            if( std::fabs( r_ji_n ) > Ri )
               return;

            RealType Rj;
            if( SPHConfig::spaceDimension == 2 )
               Rj = r_boxFactor * s_j;
            else if( SPHConfig::spaceDimension == 3 )
               Rj = r_boxFactor * sqrt( s_j );
            const VectorType r_ji_t = r_ji - r_ji_n * n_j;

            // Particle is too far from boundary element
            if( ( r_ji_t, r_ji_t ) >= Rj * Rj )
               return;

            //*r_i_mod += ( r_ji_n - Ri ) * n_j;
            view_points[ i ] += ( r_ji_n - Ri ) * n_j;
         }
      };

      auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
      {
         const VectorType r_i = view_points[ i ];
         const RealType rho_i = view_rho[ i ];
         const RealType Ri = 0.5f * pow( m / rho_i, 1.f / SPHConfig::spaceDimension );
         VectorType r_i_mod = r_i;

         ParticleSystem::NeighborsLoopAnotherSet::exec( i, r_i, searchInBound, pstLoop, Ri, &r_i_mod );

         //view_points[ i ] = r_i_mod;
      };
      Algorithms::parallelFor< DeviceType >( 0, numberOfParticles, particleLoop );
   }
};


} // SPH
} // TNL

