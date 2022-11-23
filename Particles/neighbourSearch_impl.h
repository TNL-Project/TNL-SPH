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
NeighborSearch< ParticleConfig, ParticleSystem >::particlesToCells
()
{
   GlobalIndexType nPtcs = particles->getNumberOfParticles();

   auto view_firstCellParticle = this->firstCellParticle.getView();
   auto view_lastCellParticle = this->lastCellParticle.getView();
   //auto view_particleCellIndex = this->particles->getParticleCellIndices().getData(); //works for loop
   auto view_particleCellIndex = this->particles->getParticleCellIndices().getView();

   //: //resolve first particle by hand
   //view_firstCellParticle[  view_particleCellIndex[ 0 ] ] = 0;
   view_firstCellParticle.setElement( view_particleCellIndex.getElement( 0 ) ,  0 );
   //view_lastCellParticle[  view_particleCellIndex[ 0 ] ] = ( view_particleCellIndex[ 0 ] != view_particleCellIndex[ 0+1 ] ) ? 0 : -1;
	 view_lastCellParticle.setElement( view_particleCellIndex.getElement( 0 ), ( view_particleCellIndex.getElement( 0 ) != view_particleCellIndex.getElement( 0+1 ) ) ? 0 : -1 );

   auto init = [=] __cuda_callable__ ( int i ) mutable
   {
   		if( view_particleCellIndex[ i ] != view_particleCellIndex[ i-1 ] )
         view_firstCellParticle[  view_particleCellIndex[ i ] ] = i ;
     	if( view_particleCellIndex[ i ] != view_particleCellIndex[ i+1 ] )
        view_lastCellParticle[  view_particleCellIndex[ i ] ] =  i ;
   };
   Algorithms::ParallelFor< DeviceType >::exec( 1, this->particles->getNumberOfParticles() -1, init );

   //: //resolve last particle by hand
   //view_firstCellParticle[  view_particleCellIndex[ nPtcs-1 ] ] = ( view_particleCellIndex[ nPtcs -1 ] != view_particleCellIndex[ nPtcs-2 ] ) ? nPtcs-1 : -1;
   view_firstCellParticle.setElement( view_particleCellIndex.getElement( nPtcs - 1 ), ( view_particleCellIndex.getElement( nPtcs -1 ) != view_particleCellIndex.getElement( nPtcs-2 ) ) ? nPtcs-1 : -1 );
   //view_lastCellParticle[  view_particleCellIndex[ this->particles->getNumberOfParticles()-1 ] ] = this->particles->getNumberOfParticles()-1;
   view_lastCellParticle.setElement(  view_particleCellIndex.getElement( nPtcs-1 ), nPtcs - 1 );

}


template< typename ParticleConfig, typename ParticleSystem >
__cuda_callable__
void
NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells
//(Cell centralCell, Cell neighborCell)
(LocalIndexType centralCell, LocalIndexType neighborCell)
{

   auto view_firstCellParticle = this->firstCellParticle.getView();
   auto view_lastCellParticle = this->lastCellParticle.getView();

   if(view_firstCellParticle[ neighborCell ] != -1)
   {
    auto f = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j ) mutable
    {
      if((l2Norm(particles->getPoint(i) - particles->getPoint(j)) < this->particles->getSearchRadius()) && (i != j))
      this->particles->setNeighbor(static_cast< LocalIndexType>( j ), static_cast< LocalIndexType> ( i ) ); //i is nbs, j is central!
			printf("HIT!");

      printf("Particle i: %f, %f particle j: %f, %f, l2Norm: %f\n", particles->getPoint(i)[0], particles->getPoint(i)[1], particles->getPoint(j)[0], particles->getPoint(j)[1], l2Norm(particles->getPoint(i) - particles->getPoint(j)));

    };

    Algorithms::ParallelFor2D< DeviceType >::exec(
        ( LocalIndexType ) firstCellParticle[ neighborCell ],
        ( LocalIndexType ) firstCellParticle[ centralCell ],
        ( LocalIndexType ) lastCellParticle[ neighborCell ] + 1,
        ( LocalIndexType ) lastCellParticle[ centralCell ] + 1,
        f );
    }
}

template< typename ParticleConfig, typename ParticleSystem >
__cuda_callable__
void
NeighborSearch< ParticleConfig, ParticleSystem >::runCycleOverGrid()
{
   auto view_firstCellParticle = this->firstCellParticle.getView();
   auto view_lastCellParticle = this->lastCellParticle.getView();

	 /* old, ready to remove
 	const int dim = particles->grid->getEntitiesCount( 2 );
   auto init = [=] __cuda_callable__ ( int i ) mutable
   {
     if( view_firstCellParticle[ i ] != -1 )
     {
        //const myCell centralCell = particles->grid->template getEntity< GridCell >( i );
   						//const int dim = particles->getCountOfGridCells( );
        //old: centralCell.refresh();
        //old: const typename myCell::template NeighborEntities< 2 >& neighborEntities = centralCell.getNeighborEntities();
        //#include "somethingNotNice.h"
      }
   };
   Algorithms::ParallelFor< DeviceType >::exec( 0, ParticleConfig::gridXsize*ParticleConfig::gridYsize, init );
	 */

   auto initF = [=] __cuda_callable__ ( myCell centralCell  ) mutable
	 {
		 const GlobalIndexType i = centralCell.getIndex();
		 if( view_firstCellParticle[ i ] != -1 )
		 {
		  	//#include "somethingNotNice.h"
				//if( centralCell.isBoundary() == false )

				const LocalIndexType mp = centralCell.template getNeighbourEntity< 2, -1, 1 >().getIndex();
				NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mp);

				const LocalIndexType zp = centralCell.template getNeighbourEntity< 2, 0, 1 >().getIndex();
				NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zp);

				const LocalIndexType pp = centralCell.template getNeighbourEntity< 2, 1, 1 >().getIndex();
				NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pp);

				const LocalIndexType mz = centralCell.template getNeighbourEntity< 2, -1, 0 >().getIndex();
				NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mz);

				const LocalIndexType zz = centralCell.template getNeighbourEntity< 2, 0, 0 >().getIndex();
				NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zz);

				const LocalIndexType pz = centralCell.template getNeighbourEntity< 2, 1, 0 >().getIndex();
				NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pz);

				const LocalIndexType mm = centralCell.template getNeighbourEntity< 2, -1, -1 >().getIndex();
				NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, mm);

				const LocalIndexType zm = centralCell.template getNeighbourEntity< 2, 0, -1 >().getIndex();
				NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, zm);

				const LocalIndexType pm = centralCell.template getNeighbourEntity< 2, 1, -1 >().getIndex();
				NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells(i, pm);
				printf("[%d, %d, %d, %d, %d, %d, %d, %d, %d]\n", mp, zp, pp, mz, zz, pz, mm, zm, pm);
		 }
	 };
	 particles->grid->template forInteriorEntities< 2 >( initF );

}

template< typename ParticleConfig, typename ParticleSystem >
__cuda_callable__
void
NeighborSearch< ParticleConfig, ParticleSystem >::searchForNeighbors()
{
   NeighborSearch< ParticleConfig, ParticleSystem >::particlesToCells();
   NeighborSearch< ParticleConfig, ParticleSystem >::runCycleOverGrid();
}

} // ParticleSystem
} // TNL

