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
   GlobalIndexType nPtcs = particles.getNumberOfParticles();

   auto view_firstCellParticle = this->firstCellParticle.getView();
   auto view_lastCellParticle = this->lastCellParticle.getView();

   //resolve first particle by hand
   view_firstCellParticle[  particles.getParticleCellIndex( 0 ) ] = 0;
   view_lastCellParticle[  particles.getParticleCellIndex( 0 ) ] = (particles.getParticleCellIndex( 0 ) != particles.getParticleCellIndex( 0+1 )) ? 0 : -1;

   auto init = [=] __cuda_callable__ ( int i ) mutable
   {
     if(particles.getParticleCellIndex( i ) != particles.getParticleCellIndex( i-1 ))
       view_firstCellParticle[  particles.getParticleCellIndex( i ) ] = i ;
     if(particles.getParticleCellIndex( i ) != particles.getParticleCellIndex( i+1 ))
       view_lastCellParticle[  particles.getParticleCellIndex( i ) ] =  i ;
   };
   Algorithms::ParallelFor< DeviceType >::exec( 1, this->particles.getNumberOfParticles() -1, init );

   //resolve last particle by hand
   view_firstCellParticle[  particles.getParticleCellIndex( nPtcs-1 ) ] = (particles.getParticleCellIndex( nPtcs -1 ) != particles.getParticleCellIndex( nPtcs-2 )) ? nPtcs-1 : -1;
   view_lastCellParticle[  particles.getParticleCellIndex( this->particles.getNumberOfParticles()-1 ) ] = this->particles.getNumberOfParticles()-1;

}


template< typename ParticleConfig, typename ParticleSystem >
__cuda_callable__
void
NeighborSearch< ParticleConfig, ParticleSystem >::getNeighborsFromTwoCells
//(Cell centralCell, Cell neighborCell)
(LocalIndexType centralCell, LocalIndexType neighborCell)
{

   auto f = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j ) mutable
   {
     //if(l2Norm(particles.getPoint(i) - particles.getPoint(j)))
     //add to neighbor list
   };

   Algorithms::ParallelFor2D< DeviceType, Algorithms::AsynchronousMode >::exec(
       ( LocalIndexType ) firstCellParticle[ neighborCell ],
       ( LocalIndexType ) firstCellParticle[ centralCell ],
       ( LocalIndexType ) lastCellParticle[ neighborCell ],
       ( LocalIndexType ) lastCellParticle[ centralCell ],
       f );

}

template< typename ParticleConfig, typename ParticleSystem >
__cuda_callable__
void
NeighborSearch< ParticleConfig, ParticleSystem >::runCycleOverGrid()
{

   auto init = [=] __cuda_callable__ ( int i ) mutable
   {

   };
   Algorithms::ParallelFor< DeviceType >::exec( 0, 100, init );

}

template< typename ParticleConfig, typename ParticleSystem >
__cuda_callable__
void
NeighborSearch< ParticleConfig, ParticleSystem >::searchForNeighbors()
{

}

} // ParticleSystem
} // TNL

