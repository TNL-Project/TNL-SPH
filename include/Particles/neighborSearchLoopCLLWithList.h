#include <climits>
namespace TNL {
namespace ParticleSystem {

struct NeighborsLoopCellLinkedListWithList
{
   template< typename NeighborsLoopParams, typename Function, typename... FunctionArgs >
   __cuda_callable__
   static void
   exec( typename NeighborsLoopParams::GlobalIndexType i,
         typename NeighborsLoopParams::PointType r_i,
         const NeighborsLoopParams& params,
         Function f, FunctionArgs... args )
   {
      //const int globalIdx = i * params.neighborsCountLimit;
      //const int neighborsCount = params.neighborListStorageView[ globalIdx ];

      const int neighborsCount = params.neighborListStorageView[ i ];
      const int numberOfParticles = params.numberOfParticles;

      for( int p = 0; p < neighborsCount; p++ ){
         //int j = params.neighborListStorageView[ globalIdx + p + 1 ];
         int j = params.neighborListStorageView[ numberOfParticles * ( p + 1 ) + i ];
         //if( j == INT_MAX ){ break; }
         //if( i == 0 )
         //   printf(" << i: %d, j: %d p: %d, idx: %d  >> ", i, j, p, numberOfParticles * ( p + 1 ) + i );
         if( i == j ){ continue; }
         f( i, j, r_i, args... );
      }
   }
};

struct NeighborsLoopCellLinkedListWithListAnitherSet
{
   template< typename NeighborsLoopParams, typename Function, typename... FunctionArgs >
   __cuda_callable__
   static void
   exec( typename NeighborsLoopParams::GlobalIndexType i,
         typename NeighborsLoopParams::PointType r_i,
         const NeighborsLoopParams& params,
         Function f, FunctionArgs... args )
   {
      const int globalIdx = i * params.neighborsCountLimit;
      const int particlesToSearchLabel = params.particleSetLabel;
      const int numberOfParticles = params.numberOfParticles;

      // in case we search in different set, offset to store neighbors
      int anotherSetOffset = 0;
      for( int k = 0; k < particlesToSearchLabel; k++ )
         //anotherSetOffset += params.neighborListStorageView[ globalIdx + anotherSetOffset ] + 1;
         anotherSetOffset += params.neighborListStorageView[ anotherSetOffset * numberOfParticles + i ] + 1;

      //const int neighborsCount = params.neighborListStorageView[ globalIdx + anotherSetOffset ];
      const int neighborsCount = params.neighborListStorageView[ numberOfParticles * anotherSetOffset + i ];

      for( int p = 0; p < neighborsCount; p++ ){
         //int j = params.neighborListStorageView[ globalIdx + anotherSetOffset + p + 1 ];
         int j = params.neighborListStorageView[ numberOfParticles * ( p + anotherSetOffset + 1 ) + i ];
         //if( i == 0 )
         //   printf(" << LAB: %d i: %d, j: %d p: %d, idx: %d nbc: %d asos: %d >> ", particlesToSearchLabel, i, j, p, numberOfParticles * ( p + anotherSetOffset + 1 ) + i, neighborsCount, anotherSetOffset );

         //if( i == j ){ continue; }
         f( i, j, r_i, args... );
      }
   }
};

} //namespace TNL
} //namespace Particles

