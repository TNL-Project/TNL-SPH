#include "neighbourSearch.h"

namespace TNL {
namespace ParticleSystem {

template< typename ParticleConfig, typename ParticleSystem >
const typename ParticleSystem::CellIndexArrayType&
NeighborSearch< ParticleConfig, ParticleSystem >::getCellFirstParticleList() const
{
   return firstCellParticle;
}

template< typename ParticleConfig, typename ParticleSystem >
typename ParticleSystem::CellIndexArrayType&
NeighborSearch< ParticleConfig, ParticleSystem >::getCellFirstParticleList()
{
   return firstCellParticle;
}

template< typename ParticleConfig, typename ParticleSystem >
const typename ParticleSystem::CellIndexArrayType&
NeighborSearch< ParticleConfig, ParticleSystem >::getCellLastParticleList() const
{
   return lastCellParticle;
}

template< typename ParticleConfig, typename ParticleSystem >
typename ParticleSystem::CellIndexArrayType&
NeighborSearch< ParticleConfig, ParticleSystem >::getCellLastParticleList()
{
   return lastCellParticle;
}

template< typename ParticleConfig, typename ParticleSystem >
__cuda_callable__
void
NeighborSearch< ParticleConfig, ParticleSystem >::resetListWithIndices
()
{
   auto view_firstCellParticle = this->firstCellParticle.getView();
   auto init = [=] __cuda_callable__ ( int i ) mutable
   {
      view_firstCellParticle[ i ] = INT_MAX ;
   };
   Algorithms::ParallelFor< DeviceType >::exec( 0, this->firstCellParticle.getSize(), init );
}

template< typename ParticleConfig, typename ParticleSystem >
__cuda_callable__
void
NeighborSearch< ParticleConfig, ParticleSystem >::particlesToCells
()
{
   GlobalIndexType numberOfParticles = particles->getNumberOfParticles();
   auto view_firstCellParticle = this->firstCellParticle.getView();
   auto view_lastCellParticle = this->lastCellParticle.getView();
   const auto view_particleCellIndex = this->particles->getParticleCellIndices().getView();

   //resolve first particle
   view_firstCellParticle.setElement( view_particleCellIndex.getElement( 0 ) ,  0 );
   view_lastCellParticle.setElement( view_particleCellIndex.getElement( 0 ), ( view_particleCellIndex.getElement( 0 ) != view_particleCellIndex.getElement( 0+1 ) ) ? 0 : INT_MAX );

   auto init = [=] __cuda_callable__ ( int i ) mutable
   {
      if( view_particleCellIndex[ i ] != view_particleCellIndex[ i-1 ] )
         view_firstCellParticle[  view_particleCellIndex[ i ] ] = i ;
      if( view_particleCellIndex[ i ] != view_particleCellIndex[ i+1 ] )
         view_lastCellParticle[  view_particleCellIndex[ i ] ] =  i ;
   };
   Algorithms::ParallelFor< DeviceType >::exec( 1, numberOfParticles - 1, init );

   //resolve last partile
   view_firstCellParticle.setElement( view_particleCellIndex.getElement( numberOfParticles - 1 ), ( view_particleCellIndex.getElement( numberOfParticles -1 ) != view_particleCellIndex.getElement( numberOfParticles-2 ) ) ? numberOfParticles-1 : INT_MAX );
   view_lastCellParticle.setElement(  view_particleCellIndex.getElement( numberOfParticles-1 ), numberOfParticles - 1 );
}

template< typename ParticleConfig, typename ParticleSystem >
void
NeighborSearch< ParticleConfig, ParticleSystem >::searchForNeighborsExampleLoop()
{
   const CellIndexArrayView view_firstCellParticle = this->firstCellParticle.getView();
   const CellIndexArrayView view_particleCellIndex = this->particles->getParticleCellIndices().getView();
   const auto view_points = this->particles->getPoints().getView();

   auto view_neighborsCount = this->particles->getNeighborsCountList().getView();
   auto view_neighbors = this->particles->getNeighborsList().getView();
   const RealType searchRadius = this->particles->getSearchRadius();

   GlobalIndexType numberOfParticles  = particles->getNumberOfParticles();
   static constexpr GlobalIndexType _numberOfCells = ParticleConfig::gridXsize; //FIX THIS

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i  ) mutable
   {
      const unsigned int activeCell = view_particleCellIndex[ i ];

      for( int ci = -1; ci <= 1; ci++ ){
         for( int cj = -1; cj <= 1; cj++ ){
            const unsigned int neighborCell = activeCell + cj * _numberOfCells + ci;
            int j = view_firstCellParticle[ neighborCell ];
            while( ( j < numberOfParticles ) && ( j >= 0 ) && ( view_particleCellIndex[ j ] == neighborCell ) ){

               /* START OF LOOP OVER NEIGHBROS */
               if( ( l2Norm( view_points[ i ] - view_points[ j ] ) < searchRadius ) && ( i != j ) )
               {
                  view_neighbors[ ( ParticleConfig::maxOfNeigborsPerParticle )*i + view_neighborsCount[ i ] ] = j;
                  view_neighborsCount[ i ]++;
               }
               /* END OF LOOP OVER NEIGHBROS */

               j++;
            }
         }
      }
   };
   Algorithms::ParallelFor< DeviceType >::exec( 0, numberOfParticles, particleLoop );
}

template< typename ParticleConfig, typename ParticleSystem >
void
NeighborSearch< ParticleConfig, ParticleSystem >::searchForNeighbors()
{
   NeighborSearch< ParticleConfig, ParticleSystem >::particlesToCells();

   //NeighborSearch< ParticleConfig, ParticleSystem >::searchForNeighborsExampleLoop();
   //NeighborSearch< ParticleConfig, ParticleSystem >::searchForNeighborsWithForAll();
   NeighborSearch< ParticleConfig, ParticleSystem >::searchForNeighborsWithForEach();
}

template< typename ParticleConfig, typename ParticleSystem >
template< typename Function, typename... FunctionArgs >
__cuda_callable__
void
NeighborSearch< ParticleConfig, ParticleSystem >::loopOverParticlesAndNeighbors( Function f, FunctionArgs... args )
{
   static constexpr GlobalIndexType numberOfCellsInX = ParticleSystem::Config::gridXsize; //FIXIT
   const auto view_firstCellParticle = getCellFirstParticleList().getView();
   const auto view_particleCellIndex = particles->getParticleCellIndices().getView();
   GlobalIndexType numberOfParticles = particles->getNumberOfParticles();

   //C++20 solution: auto particleLoop = [=, ... args = std::forward< FunctionArgs >( args )] __cuda_callable__ ( LocalIndexType i  ) mutable
   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i, FunctionArgs... args ) mutable
   {
      const unsigned int activeCell = view_particleCellIndex[ i ];

      for( int ci = -1; ci <= 1; ci++ ){
         for( int cj = -1; cj <= 1; cj++ ){
            const unsigned int neighborCell = activeCell + cj * numberOfCellsInX + ci;
            int j = view_firstCellParticle[ neighborCell ];

            while( ( j < numberOfParticles ) && ( view_particleCellIndex[ j ] == neighborCell ) ){
               if( i == j ){ j++; continue; }
               f( i, j, args... );
               j++;
            } //while over particle in cell
         } //for cells in y direction
      } //for cells in x direction
   };
   Algorithms::ParallelFor< DeviceType >::exec( 0, particles->getNumberOfParticles(), particleLoop, args... );
}

template< typename ParticleConfig, typename ParticleSystem >
void
NeighborSearch< ParticleConfig, ParticleSystem >::searchForNeighborsWithForAll()
{
   const auto view_points = this->particles->getPoints().getView();
   const RealType searchRadius = this->particles->getSearchRadius();
   auto view_neighborsCount = this->particles->getNeighborsCountList().getView();
   auto view_neighbors = this->particles->getNeighborsList().getView();

   auto compareTwoParticles = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j  ) mutable
   {
      if( ( l2Norm( view_points[ i ] - view_points[ j ] ) < searchRadius ) && ( i != j ) )
      {
         view_neighbors[ ( ParticleConfig::maxOfNeigborsPerParticle )*i + view_neighborsCount[ i ] ] = j;
         view_neighborsCount[ i ]++;
      }
   };
   this->loopOverParticlesAndNeighbors( compareTwoParticles );
}

template< typename ParticleConfig, typename ParticleSystem >
template< typename Function, typename... FunctionArgs >
__cuda_callable__
void
NeighborSearch< ParticleConfig, ParticleSystem >::loopOverNeighbors( const GlobalIndexType i, const GlobalIndexType& numberOfParticles, const CellIndexArrayView& view_firstCellParticle, const CellIndexArrayView& view_lastCellParticle, const CellIndexArrayView& view_particleCellIndex, Function f, FunctionArgs... args )
{
   static constexpr GlobalIndexType numberOfCellsInX = ParticleSystem::Config::gridXsize; //FIXIT
   const unsigned int activeCell = view_particleCellIndex[ i ];

   for( int ci = -1; ci <= 1; ci++ ){
      for( int cj = -1; cj <= 1; cj++ ){
         const unsigned int neighborCell = activeCell + cj * numberOfCellsInX + ci;
         int j = view_firstCellParticle[ neighborCell ];
         int j_end = view_lastCellParticle[ neighborCell ];
         if( j_end >= numberOfParticles )
          	j_end = -1;
         while( ( j <= j_end ) ){ //test
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
   const auto view_firstCellParticle = getCellFirstParticleList().getView();
   const auto view_lastCellParticle = getCellFirstParticleList().getView();
   const auto view_particleCellIndex = particles->getParticleCellIndices().getView();
   GlobalIndexType numberOfParticles = particles->getNumberOfParticles();

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
      this->loopOverNeighbors( i, numberOfParticles, view_firstCellParticle, view_lastCellParticle, view_particleCellIndex, compareTwoParticles );
   };
   Algorithms::ParallelFor< DeviceType >::exec( 0, particles->getNumberOfParticles(), particleLoop );
}

} // ParticleSystem
} // TNL

