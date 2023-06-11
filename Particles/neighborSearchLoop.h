#pragma once

namespace TNL {
namespace ParticleSystem {
namespace SPH { //this is just temp

struct NeighborsLoop
{
   template< typename NeighborsLoopParams, typename Function, typename... FunctionArgs >
   __cuda_callable__
   static void
   exec( //const typename NeighborsLoopParams::GlobalIndexType& i,
         typename NeighborsLoopParams::GlobalIndexType i,
         //const typename NeighborsLoopParams::PointType& r_i,
         typename NeighborsLoopParams::PointType r_i,
         const NeighborsLoopParams& params,
         Function f, FunctionArgs... args )
   {
      //const typename NeighborsLoopParams gridIndex = NeighborsLoopParams::CellIndexer ...
      const typename NeighborsLoopParams::IndexVectorType gridIndex = TNL::floor( ( r_i - params.gridOrigin ) / params.searchRadius );

      for( int cj = gridIndex[ 1 ] -1; cj <= gridIndex[ 1 ] + 1; cj++ ){
         for( int ci = gridIndex[ 0 ] - 1; ci <= gridIndex[ 0 ] + 1; ci++ ){
            const unsigned int neighborCell = NeighborsLoopParams::CellIndexer::EvaluateCellIndex( ci, cj, params.gridSize );
            const typename NeighborsLoopParams::PairIndexType firstLastParticle = params.view_firstLastCellParticle[ neighborCell ];
            int j = firstLastParticle[ 0 ];
            int j_end = firstLastParticle[ 1 ];
            if( j_end >= params.numberOfParticles )
             	j_end = -1;
            while( ( j <= j_end ) ){
               if( i == j ){ j++; continue; }
               f( i, j, r_i, args... ); //added r_i
               //f( params.i, j, args... );
               j++;
            } //while over particle in cell
         } //for cells in x direction
      } //for cells in y direction
   }
};

} // SPH
} // ParticleSystem
} // TNL

