#include "neighbourSearch.h"

namespace TNL {
namespace ParticleSystem {

template< typename ParticleConfig, typename ParticleSystem >
const typename NeighborSearch< ParticleConfig, ParticleSystem >::PairIndexArrayType&
NeighborSearch< ParticleConfig, ParticleSystem >::getCellFirstLastParticleList() const
{
   return firstLastCellParticle;
}

template< typename ParticleConfig, typename ParticleSystem >
typename NeighborSearch< ParticleConfig, ParticleSystem >::PairIndexArrayType&
NeighborSearch< ParticleConfig, ParticleSystem >::getCellFirstLastParticleList()
{
   return firstLastCellParticle;
}

template< typename ParticleConfig, typename ParticleSystem >
void
NeighborSearch< ParticleConfig, ParticleSystem >::resetListWithIndices
()
{
   auto view_firstLastCellParticle = this->firstLastCellParticle.getView();
   auto init = [=] __cuda_callable__ ( int i ) mutable
   {
      view_firstLastCellParticle[ i ] = INT_MAX ;
   };
   Algorithms::ParallelFor< DeviceType >::exec( 0, this->firstLastCellParticle.getSize(), init );
}

template< typename ParticleConfig, typename ParticleSystem >
void
NeighborSearch< ParticleConfig, ParticleSystem >::particlesToCells
()
{
   GlobalIndexType numberOfParticles = particles->getNumberOfParticles();
   if( numberOfParticles == 0 ) //temp
      return;
   auto view_firstLastCellParticle = this->firstLastCellParticle.getView();
   const auto view_particleCellIndex = this->particles->getParticleCellIndices().getView();

   //resolve first particle
   view_firstLastCellParticle.setElement( view_particleCellIndex.getElement( 0 ), { 0, ( view_particleCellIndex.getElement( 0 ) != view_particleCellIndex.getElement( 0+1 ) ) ? 0 : INT_MAX } ) ;

   auto init = [=] __cuda_callable__ ( int i ) mutable
   {
      if( view_particleCellIndex[ i ] != view_particleCellIndex[ i-1 ] )
         view_firstLastCellParticle[  view_particleCellIndex[ i ] ][ 0 ] = i ;
      if( view_particleCellIndex[ i ] != view_particleCellIndex[ i+1 ] )
         view_firstLastCellParticle[  view_particleCellIndex[ i ] ][ 1 ] =  i ;
   };
   Algorithms::ParallelFor< DeviceType >::exec( 1, numberOfParticles - 1, init );

   //resolve last partile
   view_firstLastCellParticle.setElement( view_particleCellIndex.getElement( numberOfParticles - 1 ), { ( view_particleCellIndex.getElement( numberOfParticles -1 ) != view_particleCellIndex.getElement( numberOfParticles-2 ) ) ? numberOfParticles-1 : INT_MAX, numberOfParticles - 1 } );

}

template< typename ParticleConfig, typename ParticleSystem >
void
NeighborSearch< ParticleConfig, ParticleSystem >::searchForNeighbors()
{
   NeighborSearch< ParticleConfig, ParticleSystem >::particlesToCells();
   NeighborSearch< ParticleConfig, ParticleSystem >::searchForNeighborsWithForEach();
}

template< typename ParticleConfig, typename ParticleSystem >
template< typename Function, typename... FunctionArgs >
__cuda_callable__
void
NeighborSearch< ParticleConfig, ParticleSystem >::loopOverNeighbors( const GlobalIndexType i, const GlobalIndexType& numberOfParticles, const GlobalIndexType& gridX, const GlobalIndexType& gridY, const PairIndexArrayView& view_firstLastCellParticle, const CellIndexArrayView& view_particleCellIndex, Function f, FunctionArgs... args )
{
   for( int cj = gridY -1; cj <= gridY + 1; cj++ ){
      for( int ci = gridX - 1; ci <= gridX + 1; ci++ ){
         const unsigned int neighborCell = ParticleSystem::CellIndexer::EvaluateCellIndex( ci, cj );
         const PairIndexType firstLastParticle= view_firstLastCellParticle[ neighborCell ];
         int j = firstLastParticle[ 0 ];
         int j_end = firstLastParticle[ 1 ];
         if( j_end >= numberOfParticles )
          	j_end = -1;
         while( ( j <= j_end ) ){
            if( i == j ){ j++; continue; }
            f( i, j, args... );
            j++;
         } //while over particle in cell
      } //for cells in x direction
   } //for cells in y direction
}

template< typename ParticleConfig, typename ParticleSystem >
template< typename Function, typename... FunctionArgs >
__cuda_callable__
void
NeighborSearch< ParticleConfig, ParticleSystem >::loopOverNeighbors( const GlobalIndexType i, const GlobalIndexType& numberOfParticles, const GlobalIndexType& gridX, const GlobalIndexType& gridY, const GlobalIndexType& gridZ, const PairIndexArrayView& view_firstLastCellParticle, const CellIndexArrayView& view_particleCellIndex, Function f, FunctionArgs... args )
{
   for( int cj = gridZ -1; cj <= gridZ + 1; cj++ ){
      for( int cj = gridY -1; cj <= gridY + 1; cj++ ){
         for( int ci = gridX - 1; ci <= gridX + 1; ci++ ){
            const unsigned int neighborCell = ParticleSystem::CellIndexer::EvaluateCellIndex( ci, cj );
            const PairIndexType firstLastParticle= view_firstLastCellParticle[ neighborCell ];
            int j = firstLastParticle[ 0 ];
            int j_end = firstLastParticle[ 1 ];
            if( j_end >= numberOfParticles )
             	j_end = -1;
            while( ( j <= j_end ) ){
               if( i == j ){ j++; continue; }
               f( i, j, args... );
               j++;
            } //while over particle in cell
         } //for cells in z direction
      } //for cells in y direction
   } //for cells in x direction
}

template< typename ParticleConfig, typename ParticleSystem >
template< typename Function, typename... FunctionArgs >
__cuda_callable__
void
NeighborSearch< ParticleConfig, ParticleSystem >::loopOverNeighbors( const GlobalIndexType i, const GlobalIndexType& numberOfParticles, const PairIndexArrayView& view_firstLastCellParticle, const CellIndexArrayView& view_particleCellIndex, Function f, FunctionArgs... args )
{
   static constexpr GlobalIndexType numberOfCellsInX = ParticleSystem::Config::gridXsize; //FIXIT
   const unsigned int activeCell = view_particleCellIndex[ i ];

   for( int ci = -1; ci <= 1; ci++ ){
      for( int cj = -1; cj <= 1; cj++ ){
         const unsigned int neighborCell = activeCell + cj * numberOfCellsInX + ci;
         //const unsigned int neighborCell = ParticleSystem::CellIndexer::EvaluateCellIndex( ... );
         const PairIndexType firstLastParticle= view_firstLastCellParticle[ neighborCell ];
         int j = firstLastParticle[ 0 ];
         int j_end = firstLastParticle[ 1 ];
         if( j_end >= numberOfParticles )
          	j_end = -1;
         while( ( j <= j_end ) ){
            if( i == j ){ j++; continue; }
            f( i, j, args... );
            j++;
         } //while over particle in cell
      } //for cells in y direction
   } //for cells in x direction
}

template< typename ParticleConfig, typename ParticleSystem >
void
NeighborSearch< ParticleConfig, ParticleSystem >::searchForNeighborsWithForEach()
{
   //needed for neighbor list
   const auto view_points = this->particles->getPoints().getView();
   const RealType searchRadius = this->particles->getSearchRadius();
   auto view_neighborsCount = this->particles->getNeighborsCountList().getView();
   auto view_neighbors = this->particles->getNeighborsList().getView();

   // needed for negihbor loop
   const auto view_firstLastCellParticle = getCellFirstLastParticleList().getView();
   const auto view_particleCellIndex = particles->getParticleCellIndices().getView();
   GlobalIndexType numberOfParticles = particles->getNumberOfParticles();

   using PointType = typename ParticleSystem::PointType;
   static constexpr RealType gridXbegin = ParticleSystem::Config::gridXbegin; //FIXIT
   static constexpr RealType gridYbegin = ParticleSystem::Config::gridYbegin; //FIXIT
   static constexpr GlobalIndexType numberOfCellsInX = ParticleSystem::Config::gridXsize; //FIXIT


   auto compareTwoParticles = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j  ) mutable
   {
      if( ( l2Norm( view_points[ i ] - view_points[ j ] ) < searchRadius ) && ( i != j ) )
      {
         view_neighbors[ ( ParticleConfig::maxOfNeigborsPerParticle )*i + view_neighborsCount[ i ] ] = j;
         view_neighborsCount[ i ]++;
      }
   };

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const PointType r_i = view_points[ i ];
      const int gridIndexI = TNL::floor( ( r_i[ 0 ] - gridXbegin ) / searchRadius );
      const int gridIndexJ = TNL::floor( ( r_i[ 1 ] - gridYbegin ) / searchRadius );
      const unsigned int neighborCell =  gridIndexJ * numberOfCellsInX + gridIndexI;
      this->loopOverNeighbors( i, numberOfParticles, gridIndexI, gridIndexJ, view_firstLastCellParticle, view_particleCellIndex, compareTwoParticles );
   };
   Algorithms::ParallelFor< DeviceType >::exec( 0, particles->getNumberOfParticles(), particleLoop );
}

} // ParticleSystem
} // TNL

