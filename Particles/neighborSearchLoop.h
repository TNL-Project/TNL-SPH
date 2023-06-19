#pragma once

namespace TNL {
namespace ParticleSystem {
namespace SPH { //this is just temp

struct NeighborsLoop2D
{
   template< typename NeighborsLoopParams, typename Function, typename... FunctionArgs >
   __cuda_callable__
   static void
   exec( typename NeighborsLoopParams::GlobalIndexType i,
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

struct NeighborsLoop2DAnotherSet
{
   template< typename NeighborsLoopParams, typename Function, typename... FunctionArgs >
   __cuda_callable__
   static void
   exec( typename NeighborsLoopParams::GlobalIndexType i,
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
               f( i, j, r_i, args... ); //added r_i
               //f( params.i, j, args... );
               j++;
            } //while over particle in cell
         } //for cells in x direction
      } //for cells in y direction
   }
};

struct NeighborsLoop3D
{
   template< typename NeighborsLoopParams, typename Function, typename... FunctionArgs >
   __cuda_callable__
   static void
   exec( typename NeighborsLoopParams::GlobalIndexType i,
         typename NeighborsLoopParams::PointType r_i,
         const NeighborsLoopParams& params,
         Function f, FunctionArgs... args )
   {
      //const typename NeighborsLoopParams gridIndex = NeighborsLoopParams::CellIndexer ...
      const typename NeighborsLoopParams::IndexVectorType gridIndex = TNL::floor( ( r_i - params.gridOrigin ) / params.searchRadius );

      for( int ck = gridIndex[ 2 ] -1; ck <= gridIndex[ 2 ] + 1; ck++ ){
         for( int cj = gridIndex[ 1 ] -1; cj <= gridIndex[ 1 ] + 1; cj++ ){
            for( int ci = gridIndex[ 0 ] - 1; ci <= gridIndex[ 0 ] + 1; ci++ ){
               const unsigned int neighborCell = NeighborsLoopParams::CellIndexer::EvaluateCellIndex( ci, cj, ck, params.gridSize );
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
      } //for cells in z direction
   }
};

struct NeighborsLoop3DAnotherSet
{
   template< typename NeighborsLoopParams, typename Function, typename... FunctionArgs >
   __cuda_callable__
   static void
   exec( typename NeighborsLoopParams::GlobalIndexType i,
         typename NeighborsLoopParams::PointType r_i,
         const NeighborsLoopParams& params,
         Function f, FunctionArgs... args )
   {
      //const typename NeighborsLoopParams gridIndex = NeighborsLoopParams::CellIndexer ...
      const typename NeighborsLoopParams::IndexVectorType gridIndex = TNL::floor( ( r_i - params.gridOrigin ) / params.searchRadius );

      for( int ck = gridIndex[ 2 ] -1; ck <= gridIndex[ 2 ] + 1; ck++ ){
         for( int cj = gridIndex[ 1 ] -1; cj <= gridIndex[ 1 ] + 1; cj++ ){
            for( int ci = gridIndex[ 0 ] - 1; ci <= gridIndex[ 0 ] + 1; ci++ ){
               const unsigned int neighborCell = NeighborsLoopParams::CellIndexer::EvaluateCellIndex( ci, cj, ck, params.gridSize );
               const typename NeighborsLoopParams::PairIndexType firstLastParticle = params.view_firstLastCellParticle[ neighborCell ];
               int j = firstLastParticle[ 0 ];
               int j_end = firstLastParticle[ 1 ];
               if( j_end >= params.numberOfParticles )
                	j_end = -1;
               while( ( j <= j_end ) ){
                  f( i, j, r_i, args... ); //added r_i
                  //f( params.i, j, args... );
                  j++;
               } //while over particle in cell
            } //for cells in x direction
         } //for cells in y direction
      } //for cells in z direction
   }
};

struct NeighborsLoop
{
   template< typename NeighborsLoopParams, typename Function, typename... FunctionArgs >
   __cuda_callable__
   static void
   exec( typename NeighborsLoopParams::GlobalIndexType i,
         typename NeighborsLoopParams::PointType r_i,
         const NeighborsLoopParams& params,
         Function f, FunctionArgs... args )
   {
      if constexpr( NeighborsLoopParams::PointType::getSize() == 1 ) {
         static_assert( NeighborsLoopParams::PointType::getSize() == 1, "loop over neighbors is not implemented for 1 dimension" );
      }
      else if constexpr( NeighborsLoopParams::PointType::getSize() == 2 ) {
         NeighborsLoop2D::exec( i, r_i, params, f, args... );
      }
      else if constexpr( NeighborsLoopParams::PointType::getSize() == 3 ) {
         NeighborsLoop3D::exec( i, r_i, params, f, args... );
      }
      else {
         static_assert( NeighborsLoopParams::PointType::getSize() <= 3, "loop over neighbors is not implemented yet for 4 or more dimensions" );
      }
   }
};

struct NeighborsLoopAnotherSet
{
   template< typename NeighborsLoopParams, typename Function, typename... FunctionArgs >
   __cuda_callable__
   static void
   exec( typename NeighborsLoopParams::GlobalIndexType i,
         typename NeighborsLoopParams::PointType r_i,
         const NeighborsLoopParams& params,
         Function f, FunctionArgs... args )
   {
      if constexpr( NeighborsLoopParams::PointType::getSize() == 1 ) {
         static_assert( NeighborsLoopParams::PointType::getSize() == 1, "loop over neighbors is not implemented for 1 dimension" );
      }
      else if constexpr( NeighborsLoopParams::PointType::getSize() == 2 ) {
         NeighborsLoop2DAnotherSet::exec( i, r_i, params, f, args... );
      }
      else if constexpr( NeighborsLoopParams::PointType::getSize() == 3 ) {
         NeighborsLoop3DAnotherSet::exec( i, r_i, params, f, args... );
      }
      else {
         static_assert( NeighborsLoopParams::PointType::getSize() <= 3, "loop over neighbors is not implemented yet for 4 or more dimensions" );
      }
   }
};

} // SPH
} // ParticleSystem
} // TNL

