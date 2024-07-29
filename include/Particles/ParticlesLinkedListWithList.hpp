#include "ParticlesLinkedList.h"
#include "details/details.h"

namespace TNL {
namespace ParticleSystem {

template < typename ParticleConfig, typename DeviceType >
void
ParticlesLinkedList< ParticleConfig, DeviceType >::buildParticleList()
{
   auto neighborList_view = this->neighborList.getView();
   auto compareDistance = [=] __cuda_callable__ ( LocalIndexType i,
                                                  LocalIndexType j,
                                                  VectorType& r_i
                                                  GlobalIndexType& globIdxStart,
                                                  GlobalIndexType* neihgborsCount ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         globalIdx = globIdxStart + *neihgborsCount;
         neighborList_view[ globalIdx + 1 ] = j;
         *neihgborsCount++;
      }
   };

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points_bound[ i ];
      const GlobalIndexType globIdxStart = i * neighborsCountLimit;
      GlobalIndexType neihgborsCount = 0;

      TNL::ParticleSystem::NeighborsLoopAnotherSet::exec(
            i, r_i, searchInFluid, compareDistance, globIdxStart, neihgborsCount );

      neighborList_view[ globalIdx ] = neihgborsCount;
   };
   this->particles->forAll( particleLoop );
}

template< typename NeighborsLoopParams, typename Function, typename... FunctionArgs >
__cuda_callable__
static void
exec( typename NeighborsLoopParams::GlobalIndexType i,
      typename NeighborsLoopParams::PointType r_i,
      const NeighborsLoopParams& params,
      Function f, FunctionArgs... args )
{
   const GlobalIndexType globalIdx = i * params.neighborsCountLimit;
   const GlobalIndexType neighborsCount = params.neighborListView[ globalIdx ];
   for( int p = 0; p < neighborsCount; p++ ){
      int j = params.neighborListView[ globalIdx + p ];
      if( i == j ){ continue; }
      f( i, j, r_i, args... );
   }
}

} //namespace TNL
} //namespace Particles
