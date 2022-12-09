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
   GlobalIndexType nPtcs = particles->getNumberOfParticles();

   auto view_firstCellParticle = this->firstCellParticle.getView();
   auto view_lastCellParticle = this->lastCellParticle.getView();
   //auto view_particleCellIndex = this->particles->getParticleCellIndices().getData(); //works for loop
   auto view_particleCellIndex = this->particles->getParticleCellIndices().getView();

	 //RESET LIST
	 //view_particleCellIndex = INT_MAX;

   //: //resolve first particle by hand
   //view_firstCellParticle[  view_particleCellIndex[ 0 ] ] = 0;
   view_firstCellParticle.setElement( view_particleCellIndex.getElement( 0 ) ,  0 );
   //view_lastCellParticle[  view_particleCellIndex[ 0 ] ] = ( view_particleCellIndex[ 0 ] != view_particleCellIndex[ 0+1 ] ) ? 0 : -1;
	 view_lastCellParticle.setElement( view_particleCellIndex.getElement( 0 ), ( view_particleCellIndex.getElement( 0 ) != view_particleCellIndex.getElement( 0+1 ) ) ? 0 : INT_MAX );

   auto init = [=] __cuda_callable__ ( int i ) mutable
   {
   		if( view_particleCellIndex[ i ] != view_particleCellIndex[ i-1 ] )
         view_firstCellParticle[  view_particleCellIndex[ i ] ] = i ;
     	if( view_particleCellIndex[ i ] != view_particleCellIndex[ i+1 ] )
        view_lastCellParticle[  view_particleCellIndex[ i ] ] =  i ;
   };
   Algorithms::ParallelFor< DeviceType >::exec( 1, this->particles->getNumberOfParticles() - 1, init );

   //: //resolve last particle by hand
   //view_firstCellParticle[  view_particleCellIndex[ nPtcs-1 ] ] = ( view_particleCellIndex[ nPtcs -1 ] != view_particleCellIndex[ nPtcs-2 ] ) ? nPtcs-1 : -1;
   view_firstCellParticle.setElement( view_particleCellIndex.getElement( nPtcs - 1 ), ( view_particleCellIndex.getElement( nPtcs -1 ) != view_particleCellIndex.getElement( nPtcs-2 ) ) ? nPtcs-1 : INT_MAX );
   //view_lastCellParticle[  view_particleCellIndex[ this->particles->getNumberOfParticles()-1 ] ] = this->particles->getNumberOfParticles()-1;
   view_lastCellParticle.setElement(  view_particleCellIndex.getElement( nPtcs-1 ), nPtcs - 1 );

}


template< typename ParticleConfig, typename ParticleSystem >
__cuda_callable__
void
NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells
//( LocalIndexType centralCell, LocalIndexType neighborCell )
//( LocalIndexType centralCell, LocalIndexType neighborCell, CellIndexArrayView view_firstCellParticle , CellIndexArrayView view_lastCellParticle, PointTypeArrayView view_particles )
(LocalIndexType centralCell, LocalIndexType neighborCell, CellIndexArrayView view_firstCellParticle , CellIndexArrayView view_lastCellParticle, PointTypeArrayView view_particles, NeighborsCountView view_neighborsCount, NeighborsView view_neighbors,  RealType searchRadius)
{
	 /**
		* Check this.
		* Maybe NeighborSearch::,...->getCellFirstParticleList().getView()
		* could do the work.
		*/
   //auto view_firstCellParticle = this->getCellFirstParticleList().getView();
   //auto view_lastCellParticle = this->getCellLastParticleList().getView();
   //auto view_particles = this->particles->getPoints().getView();
	 //RealType searchRadius = this->particles->getSearchRadius();

   if(view_firstCellParticle[ neighborCell ] != -1)
   {
	 	for( LocalIndexType i = ( LocalIndexType ) view_firstCellParticle[ neighborCell ];
	 			 i < view_lastCellParticle[ neighborCell ] + 1;
	 			 i++ )
	 	for( LocalIndexType j = ( LocalIndexType ) view_firstCellParticle[ centralCell ];
	 			 j < view_lastCellParticle[ centralCell ] + 1;
	 			 j++ )
    //auto f = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j ) mutable
    {
     	 //printf("Particle i: %f, %f particle j: %f, %f, l2Norm: %f\n", view_particles[ i ][0], view_particles[ i ][1], view_particles[ j ][0], view_particles[ j ][1], l2Norm(view_particles[ i ] - view_particles[ j ]));
       if( ( l2Norm( view_particles[ i ] - view_particles[ j ] ) < searchRadius ) && ( i != j ) )
			 {
        	//this->particles->setNeighbor(static_cast< LocalIndexType>( j ), static_cast< LocalIndexType> ( i ) ); //i is nbs, j is central!
				 	/* Temp. This is done directly in particle method, but needs to be fixed. */
  				view_neighbors[ ( ParticleConfig::maxOfNeigborsPerParticle )*j + view_neighborsCount[ j ] ] = i; //i is nbs, j is central!
  				view_neighborsCount[ j ]++;
			 }
    }//;
	 	//printf(" view_firstCellParticle[ neighborCell ]: %d, view_firstCellParticle[ centralCell ]: %d, view_lastCellParticle[ neighborCell ] + 1: %d, LocalIndexType ) view_lastCellParticle[ centralCell ] + 1: %d\n", ( LocalIndexType ) view_firstCellParticle[ neighborCell ] , ( LocalIndexType ) view_firstCellParticle[ centralCell ], ( LocalIndexType ) view_lastCellParticle[ neighborCell ] + 1, ( LocalIndexType ) view_lastCellParticle[ centralCell ] + 1 );

	 	/* Too many loops. */
    //Algorithms::ParallelFor2D< DeviceType >::exec(
    //    ( LocalIndexType ) view_firstCellParticle[ neighborCell ],
    //    ( LocalIndexType ) view_firstCellParticle[ centralCell ],
    //    ( LocalIndexType ) view_lastCellParticle[ neighborCell ] + 1,
    //    ( LocalIndexType ) view_lastCellParticle[ centralCell ] + 1,
    //    f );
   }
}

template< typename ParticleConfig, typename ParticleSystem >
//__cuda_callable__
void
NeighborSearch< ParticleConfig, ParticleSystem >::runCycleOverGrid()
{
   //const auto view_firstCellParticle = this->firstCellParticle.getView();
   //const auto view_lastCellParticle = this->lastCellParticle.getView();
   const CellIndexArrayView view_firstCellParticle = this->firstCellParticle.getView();
   const CellIndexArrayView view_lastCellParticle = this->lastCellParticle.getView();
	 const PointTypeArrayView view_points = this->particles->getPoints().getView();

   auto view_particleCellIndex = this->particles->getParticleCellIndices().getView();

	 /* TMEP */
	 NeighborsCountView view_neighborsCount = this->particles->getNeighborsCountList().getView();
	 NeighborsView view_neighbors = this->particles->getNeighborsList().getView();
	 RealType searchRadius = this->particles->getSearchRadius();

	 /*
   auto initF = [=] __cuda_callable__ ( myCell centralCell  ) mutable
	 {
		 const GlobalIndexType i = centralCell.getIndex();
		 if( view_firstCellParticle[ i ] != -1 )
		 {
		  	//#include "somethingNotNice.h"
				//if( centralCell.isBoundary() == false )

				const LocalIndexType mp = centralCell.template getNeighbourEntity< 2, -1, 1 >().getIndex();
				//NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mp);
				this->getNeighborsFromTwoCells(i, mp, view_firstCellParticle, view_lastCellParticle, view_points, view_neighborsCount, view_neighbors, searchRadius );
				//this->getNeighborsFromTwoCells(i, mp );

				const LocalIndexType zp = centralCell.template getNeighbourEntity< 2, 0, 1 >().getIndex();
				//NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zp);
				this->getNeighborsFromTwoCells(i, zp, view_firstCellParticle, view_lastCellParticle, view_points, view_neighborsCount, view_neighbors, searchRadius );
				//this->getNeighborsFromTwoCells(i, zp );

				const LocalIndexType pp = centralCell.template getNeighbourEntity< 2, 1, 1 >().getIndex();
				//NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pp);
				this->getNeighborsFromTwoCells(i, pp, view_firstCellParticle, view_lastCellParticle, view_points, view_neighborsCount, view_neighbors, searchRadius );
				//this->getNeighborsFromTwoCells(i, pp );

				const LocalIndexType mz = centralCell.template getNeighbourEntity< 2, -1, 0 >().getIndex();
				//NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mz);
				this->getNeighborsFromTwoCells(i, mz, view_firstCellParticle, view_lastCellParticle, view_points, view_neighborsCount, view_neighbors, searchRadius );
				//this->getNeighborsFromTwoCells(i, mz );

				const LocalIndexType zz = centralCell.template getNeighbourEntity< 2, 0, 0 >().getIndex();
				//NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zz);
				this->getNeighborsFromTwoCells(i, zz, view_firstCellParticle, view_lastCellParticle, view_points, view_neighborsCount, view_neighbors, searchRadius );
				//this->getNeighborsFromTwoCells(i, zz );

				const LocalIndexType pz = centralCell.template getNeighbourEntity< 2, 1, 0 >().getIndex();
				//NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pz);
				this->getNeighborsFromTwoCells(i, pz, view_firstCellParticle, view_lastCellParticle, view_points, view_neighborsCount, view_neighbors, searchRadius );
				//this->getNeighborsFromTwoCells(i, pz );

				const LocalIndexType mm = centralCell.template getNeighbourEntity< 2, -1, -1 >().getIndex();
				//NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mm);
				this->getNeighborsFromTwoCells(i, mm, view_firstCellParticle, view_lastCellParticle, view_points, view_neighborsCount, view_neighbors, searchRadius );
				//this->getNeighborsFromTwoCells(i, mm );

				const LocalIndexType zm = centralCell.template getNeighbourEntity< 2, 0, -1 >().getIndex();
				//NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zm);
				this->getNeighborsFromTwoCells(i, zm, view_firstCellParticle, view_lastCellParticle, view_points, view_neighborsCount, view_neighbors, searchRadius );
				//this->getNeighborsFromTwoCells(i, zm );

				const LocalIndexType pm = centralCell.template getNeighbourEntity< 2, 1, -1 >().getIndex();
				//NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pm);
				this->getNeighborsFromTwoCells(i, pm, view_firstCellParticle, view_lastCellParticle, view_points, view_neighborsCount, view_neighbors, searchRadius);
				//this->getNeighborsFromTwoCells(i, pm );
				//printf("[%d, %d, %d, %d, %d, %d, %d, %d, %d]\n", mp, zp, pp, mz, zz, pz, mm, zm, pm);
		 }
	 };
	 particles->grid->template forInteriorEntities< 2 >( initF );
	 */

	 GlobalIndexType N  = particles->getNumberOfParticles();
   static constexpr GlobalIndexType _numberOfCells = ParticleConfig::gridXsize;

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i  ) mutable
	 {

		 const unsigned int activeCell = view_particleCellIndex[ i ];

		 for( int ci = -1; ci <= 1; ci++ ){
			 for( int cj = -1; cj <= 1; cj++ ){

				 const unsigned int neighborCell = activeCell + cj * _numberOfCells + ci;
				 int j = view_firstCellParticle[ neighborCell ]; //USE INT MAX

					//printf("Particle i = %d in cell %d and particle j = %d in cell %d.\n", i, activeCell, j, neighborCell);

				 //while( ( j < N ) && ( j >= 0 ) && ( view_particleCellIndex[ j ] == neighborCell ) ){
				 while( ( j < N ) && ( j >= 0 ) && ( view_particleCellIndex[ j ] == neighborCell ) ){

     	 		//printf("Particle i: %f, %f particle j: %f, %f, l2Norm: %f\n", view_points[ i ][0], view_points[ i ][1], view_points[ j ][0], view_points[ j ][1], l2Norm(view_points[ i ] - view_points[ j ]));

					 	// --- DO STUFF ---
       			if( ( l2Norm( view_points[ i ] - view_points[ j ] ) < searchRadius ) && ( i != j ) )
			 			{
       			 	//this->particles->setNeighbor(static_cast< LocalIndexType>( j ), static_cast< LocalIndexType> ( i ) ); //i is nbs, j is central!
							/*
  		 			 	view_neighbors[ ( ParticleConfig::maxOfNeigborsPerParticle )*j + view_neighborsCount[ j ] ] = i; //i is nbs, j is central!
  		 			 	view_neighborsCount[ j ]++;
							*/
  		 			 	view_neighbors[ ( ParticleConfig::maxOfNeigborsPerParticle )*i + view_neighborsCount[ i ] ] = j; //i is nbs, j is central! - nope
  		 			 	view_neighborsCount[ i ]++;
			 			}
						// --- PERFORM INTERACTIONS ---

					 j++;
				 }

			 }
		 }


	 };
	 Algorithms::ParallelFor< DeviceType >::exec( 0, particles->getNumberOfParticles(), particleLoop );

}

template< typename ParticleConfig, typename ParticleSystem >
//__cuda_callable__
void
NeighborSearch< ParticleConfig, ParticleSystem >::searchForNeighbors()
{
   NeighborSearch< ParticleConfig, ParticleSystem >::particlesToCells();
	 cudaDeviceSynchronize();
   NeighborSearch< ParticleConfig, ParticleSystem >::runCycleOverGrid();
}

template< typename ParticleConfig, typename ParticleSystem >
template< typename Function, typename... FunctionArgs >
__cuda_callable__
void
NeighborSearch< ParticleConfig, ParticleSystem >::neighborParticleLoop( GlobalIndexType i, Function f, FunctionArgs... args )
{
   static constexpr GlobalIndexType numberOfCellsInX = ParticleSystem::Config::gridXsize; //FIXIT
	 const auto view_firstCellParticle = this->neighborSearch->getCellFirstParticleList().getView();
	 const auto view_particleCellIndex = particles->getParticleCellIndices().getView();
	 GlobalIndexType numberOfParticles = particles->getNumberOfParticles();

	 const unsigned int activeCell = view_particleCellIndex[ i ];
	 for( int ci = -1; ci <= 1; ci++ ){
	    for( int cj = -1; cj <= 1; cj++ ){

		 		 const unsigned int neighborCell = activeCell + cj * numberOfCellsInX + ci;
		 		 int j = view_firstCellParticle[ neighborCell ];

		 		 while( ( j < numberOfParticles ) && ( j >= 0 ) && ( view_particleCellIndex[ j ] == neighborCell ) ){
					  if( i == j ){ j++; continue; }
				 		f( args... );
		 		 		j++;
		 		 } //while over particle in cell
			} //for cells in y direction
	 } //for cells in x direction
}

template< typename ParticleConfig, typename ParticleSystem >
void
NeighborSearch< ParticleConfig, ParticleSystem >::searchForNeighborsClean()
{
	 const auto view_points = this->particles->getPoints().getView();
	 const RealType searchRadius = this->particles->getSearchRadius();
	 NeighborsCountView view_neighborsCount = this->particles->getNeighborsCountList().getView();
	 NeighborsView view_neighbors = this->particles->getNeighborsList().getView();

   auto compareTwoParticles = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j  ) //mutable
	 {
    	if( ( l2Norm( view_points[ i ] - view_points[ j ] ) < searchRadius ) && ( i != j ) )
			{
  			 view_neighbors[ ( ParticleConfig::maxOfNeigborsPerParticle )*i + view_neighborsCount[ i ] ] = j;
  			 view_neighborsCount[ i ]++;
			}
	 };

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i  ) //mutable
	 {
		 this->neighborParticleLoop( i, compareTwoParticles );
	 };
	 Algorithms::ParallelFor< DeviceType >::exec( 0, particles->getNumberOfParticles(), particleLoop );
}

} // ParticleSystem
} // TNL

