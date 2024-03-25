#include "GhostZone.h"
#include <TNL/Meshes/Grid.h>

namespace TNL {
namespace ParticleSystem {

template< typename ParticleSystem >
class DistributedParticleSystem
{
public:

   using ParticleSystemType = ParticleSystem;
   using DeviceType = typename ParticleSystem::Device;
   using GlobalIndexType = typename ParticleSystem::GlobalIndexType;
   using RealType = typename ParticleSystem::RealType;
   using PointType = typename ParticleSystem::PointType;
   using IndexVectorType = typename ParticleSystem::IndexVectorType;
   using IndexArrayType = typename ParticleSystem::IndexArrayType;
   using ParticleZone = ParticleZone< typename ParticleSystem::Config >;
   using GridType = TNL::Meshes::Grid< typename ParticleSystem::Config::spaceDimension, RealType, DeviceType, GlobalIndexType >;
   using DistributedGridType = DistributedMesh< GridType >;

   //initialize
   void
   initialize( unsigned int numberOfParticles,
               unsigned int numberOfAllocatedParticles,
               RealType searchRadius,
               IndexVectorType gridSize,
               PointType gridOrigin,
               IndexVectorType domainDecomposition )
   {
      //set local particle system
      this->localParticles->setSize( numberOfAllocatedParticles );
      this->localParticles->setSearchRadius( searchRadius );
      this->localParticles->setGridSize( gridSize );
      this->localParticles->setGridOrigin( gridOrigin );
      this->localParticles->setNumberOfParticles( numberOfParticles );
      this->localParticles->setFirstActiveParticle( 0 );
      this->localParticles->setLastActiveParticle( numberOfParticles - 1 );
      //set distributed grid
      this->distributedGrid.setDomainDecomposition( domainDecomposition );

   }

   void
   seDistributedGrid( const DistributedGridType& distributedGrid ) {}

   [[nodiscard]] DistributedGridType&
   getDistributedGrid()
   {
      return this->distributedGrid;
   }

   [[nodiscard]] const DistributedGridType&
   getDistributedGrid() const
   {
      return this->distributedGrid;
   }

   void
   setDistributedGridParameters( const IndexVectorType& globalGridSize,
                                 const PointType& globalGridOrigin,
                                 const IndexVectorType& domainDecomposition )
   {
      //TODO: Grid should be pobably pointer...

      //int size = 10;
      int rank = TNL::MPI::GetRank();
      int nproc = TNL::MPI::GetSize();

      PointType globalOrigin;
      PointType globalProportions;
      GridType globalGrid;

      //globalGrid.setDimensions( size, size );
      globalGrid.setDimensions( globalGridSize );
      globalGrid.setDomain( globalOrigin, globalGridSize );

      //distributedGrid = new DistributedGridType();
      //distributedGrid->setDomainDecomposition( typename DistributedGridType::CoordinatesType( 3, 3 ) );

      distributedGrid.setDomainDecomposition( domainDecomposition );
      distributedGrid.setGlobalGrid( globalGrid );

      typename DistributedGridType::SubdomainOverlapsType lowerOverlap, upperOverlap;
      SubdomainOverlapsGetter< GridType >::getOverlaps( distributedGrid, lowerOverlap, upperOverlap, 1 );
      distributedGrid.setOverlaps( lowerOverlap, upperOverlap );
   }

   //initialize innerOverlpas
   void
   initializeInnerOverlaps( const GlobalIndexType& gridSize, const PointType& gridOrigin )
   {

   }

   //collect particles to innerOverlaps
   void
   collectParticlesInInnerOverlaps()
   {
      const int* neighbors = this->distributedGrid->getNeighbors();
      for( int i = 0; i < this->distributedGrid->getNeighborsCount(); i++ ) {
         //TODO: We shoud limit ourselves only to filled zones to save the call time
         innerOverlaps[ i ].updateParticlesInZone();
      }
   }

   void
   mergeTwoOverlaps( const GlobalIndexType& recvOverlapIdx, const GlobalIndexType& sendOverlapIdx )
   {
      auto recvOverlapParticlesInZone_view = innerOverlaps[ recvOverlapIdx ].getParticlesInZone().getView();
      const int recvOverlapNumberOfZoneParticles = innerOverlaps[ recvOverlapIdx ].getNumberOfParticles();
      const auto sendOverlapParticlesInZone_view = innerOverlaps[ sendOverlapIdx ].getParticlesInZone().getConstView();
      const int sendOverlapNumberOfZoneParticles = innerOverlaps[ sendOverlapIdx ].getNumberOfParticles();

      auto copyIndices = [=] __cuda_callable__ ( int i ) mutable
      {
         recvOverlapParticlesInZone_view[ recvOverlapNumberOfZoneParticles + i ] = sendOverlapParticlesInZone_view[ i ];

      };
      Algorithms::parallelFor< DeviceType >( 0, sendOverlapNumberOfZoneParticles, copyIndices );
      //TODO: Add function allowing to increase number of particles in zone
      innerOverlaps[ recvOverlapIdx ].updateNumberOfParticlesInZone(
            recvOverlapNumberOfZoneParticles + sendOverlapNumberOfZoneParticles );

   }

   //linearize indicse inside innerOverlaps
   void
   linearize()
   {
      // 6 5 7
      // 0 * 1
      // 2 3 4
      const int* neighbors = this->distributedGrid->getNeighbors();

      //add up to segments
      //TODO: Add condition that currently only 2D distribution is allowd
      //TODO: Rewrite this in terms of dimension marks
      //add 2, 4 to 3
      mergeTwoOverlaps( 3, 2 );
      mergeTwoOverlaps( 3, 4 );
      //add 4, 7 to 1
      mergeTwoOverlaps( 1, 4 );
      mergeTwoOverlaps( 1, 7 );
      //add 7, 6 to 5
      mergeTwoOverlaps( 5, 6 );
      mergeTwoOverlaps( 5, 7 );
      //add 6, 2 to 0
      mergeTwoOverlaps( 0, 6 );
      mergeTwoOverlaps( 0, 6 );

      //copy all to innerOverlapsLinearized
      for( int i = 0; i < this->distributedGrid->getNeighborsCount(); i++ ) {
         const auto zoneParticleIndices_view = innerOverlaps.getParticlesInZone().getConstView();
         const int numberOfZoneParticles = innerOverlaps[ i ].getNumberOfParticles();
         if( numberOfZoneParticles == 0 )
            continue;
         auto innerOverlapsLinearized_view = innerOverlapsLinearized.getView();
         const int offset = this->numberOfParticlesInOverlaps;

         auto copyIndices = [=] __cuda_callable__ ( int i ) mutable
         {
            innerOverlapsLinearized_view[ offset + i ] = zoneParticleIndices_view[ i ];
         };
         Algorithms::parallelFor< DeviceType >( 0, numberOfZoneParticles, copyIndices );
         this->numberOfParticlesInOverlaps += numberOfZoneParticles;
      }
   }

protected:

   ParticleSystem localParticles;
   DistributedGridType distributedGrid;
   Containers::Array< GlobalIndexType, Devices::Host, int > innerOverlaps;

   GlobalIndexType numberOfParticlesInOverlaps;
   IndexArrayType innerOverlapsLinearized;


};

}  //namespace ParticleSystem
}  //namespace TNL

