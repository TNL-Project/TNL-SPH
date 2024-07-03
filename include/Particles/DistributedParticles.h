#include "GhostZone.h"
#include <TNL/Meshes/Grid.h>
#include <TNL/Meshes/DistributedMeshes/DistributedGrid.h>
#include <TNL/Meshes/DistributedMeshes/SubdomainOverlapsGetter.h>
#include <TNL/Meshes/DistributedMeshes/Directions.h>
#include <cmath>

namespace TNL {
namespace ParticleSystem {

using namespace TNL::Meshes::DistributedMeshes;

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
   using ParticleZoneType = ParticleZone< typename ParticleSystem::Config >;
   using ParticleZonePointerType = typename Pointers::SharedPointer< ParticleZoneType >;
   using GridType = TNL::Meshes::Grid< ParticleSystem::Config::spaceDimension, RealType, DeviceType, GlobalIndexType >;
   using DistributedGridType = TNL::Meshes::DistributedMeshes::DistributedMesh< GridType >;
   using SubdomainCoordinates = Containers::StaticVector< 2, int >;

   [[nodiscard]] static constexpr int
   getNeighborsCount()
   {
      return DistributedGridType::getNeighborsCount();
   }

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
   setDistributedGrid( const DistributedGridType& distributedGrid ) {}

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

   [[nodiscard]] const Containers::Array< ParticleZoneType, Devices::Host, int >&
   getInnerOverlaps() const
   {
      return this->innerOverlaps;
   }

   [[nodiscard]] Containers::StaticArray< DistributedGridType::getNeighborsCount(), GlobalIndexType >&
   getSubdomainsParticlesCount()
   {
      return this->subdomainsParticlesCount;
   }

   [[nodiscard]] Containers::StaticArray< DistributedGridType::getNeighborsCount(), float >&
   getSubdomainsCompTime()
   {
      return this->subdomainsCompTime;
   }

   const GlobalIndexType
   getNumberOfParticlesForLoadBalancing() const
   {
      return this->subdomainParticlesCount;
   }

   void
   setNumberOfParticlesForLoadBalancing( const GlobalIndexType numberOfParticles )
   {
      this->subdomainParticlesCount = numberOfParticles;
   }

   const GlobalIndexType
   getCompTime() const
   {
      return this->subdomainCompTime;
   }

   void
   setCompTimeForLoadBalancing( const RealType compTime )
   {
      this->subdomainsCompTime = compTime;
   }

   void
   setParticlesCountResizeTrashold( const GlobalIndexType trashold )
   {
      this->particlesCountResizeTrashold = trashold;
   }

   void
   setCompTimeResizePercetnageTrashold( const RealType trashold )
   {
      this->computationalTimeResizeTrashold = trashold;
   }

   //[[nodiscard]] const MPI::Comm&
   //getCommunicator() const
   //{
   //   return communicator;
   //}

   [[nodiscard]] MPI::Comm&
   getCommunicator()
   {
      return communicator;
   }


   //TODO: 1D and 2D decompositions should be separated, currently we asume only 1D
   //initialize innerOverlpas
   void
   initializeInnerOverlaps( const GlobalIndexType numberOfOverlapsLayers, const int numberOfParticlesPerCell = 25 )
   {
      const IndexVectorType localGridDimensions = distributedGrid.getLocalMesh().getDimensions();
      const IndexVectorType increaseLocalGridSizeDueToOverlaps = 2 * numberOfOverlapsLayers;
      const IndexVectorType localGridDimensionsWithOverlap = distributedGrid.getLocalMesh().getDimensions() + increaseLocalGridSizeDueToOverlaps;
      const GlobalIndexType zoneWidth = 1 + numberOfOverlapsLayers;

      const int* neighbors = this->getDistributedGrid().getNeighbors();
      for( int i = 0; i < this->getDistributedGrid().getNeighborsCount(); i++ ) {
         if( neighbors[ i ] != -1 ){
            Containers::StaticVector< 2, int > direction = Directions::template getXYZ< 2 >( i );
            // decode overlaps dimensions
            IndexVectorType zoneOriginIdx;
            IndexVectorType zoneDimensions;
            for( int j = 0; j < 2; j ++ ){
               // assign zone origin
               if( direction[ j ] == 1 )
                  zoneOriginIdx[ j ] = localGridDimensions[ j ];
               else
                  zoneOriginIdx[ j ] = 0;

               // assign zone dimensions
               if( direction[ j ] == 0 )
                  zoneDimensions[ j ] = localGridDimensions[ j ] + zoneWidth;
               else
                  zoneDimensions[ j ] = zoneWidth;
            }
            // set overlaps
            innerOverlaps[ i ].setNumberOfParticlesPerCell( numberOfParticlesPerCell );
            innerOverlaps[ i ].assignCells( zoneOriginIdx, zoneDimensions, localGridDimensionsWithOverlap );
         }
      }
   }

   void
   initializeInnerOverlaps3D( const GlobalIndexType numberOfOverlapsLayers, const int numberOfParticlesPerCell = 125 )
   {
      const IndexVectorType localGridDimensions = distributedGrid.getLocalMesh().getDimensions();
      const IndexVectorType increaseLocalGridSizeDueToOverlaps = 2 * numberOfOverlapsLayers;
      const IndexVectorType localGridDimensionsWithOverlap = distributedGrid.getLocalMesh().getDimensions() + increaseLocalGridSizeDueToOverlaps;
      const GlobalIndexType zoneWidth = 1 + numberOfOverlapsLayers;

      const int* neighbors = this->getDistributedGrid().getNeighbors();
      for( int i = 0; i < this->getDistributedGrid().getNeighborsCount(); i++ ) {
         innerOverlaps[ i ].resetParticles();
         innerOverlaps[ i ].resetZoneCells();
         if( neighbors[ i ] != -1 ){
            //Containers::StaticVector< 3, int > direction = Directions::template getXYZ< 3 >( i );
            Containers::StaticVector< 2, int > direction = Directions::template getXYZ< 2 >( i );
            // decode overlaps dimensions
            IndexVectorType zoneOriginIdx;
            IndexVectorType zoneDimensions;
            for( int j = 0; j < 2; j ++ ){
               // assign zone origin
               if( direction[ j ] == 1 )
                  zoneOriginIdx[ j ] = localGridDimensions[ j ]; //added - 1
                  //zoneOriginIdx[ j ] = localGridDimensionsWithOverlap[ j ] - 2; //this is added
               else
                  zoneOriginIdx[ j ] = 0;

               // assign zone dimensions
               if( direction[ j ] == 0 )
                  //zoneDimensions[ j ] = localGridDimensions[ j ] + zoneWidth;
                  zoneDimensions[ j ] = localGridDimensionsWithOverlap[ j ];
               else
                  zoneDimensions[ j ] = zoneWidth;
            }
            zoneOriginIdx[ 2 ] = 0;
            zoneDimensions[ 2 ] = localGridDimensionsWithOverlap[ 2 ];
            // set overlaps
            innerOverlaps[ i ].setNumberOfParticlesPerCell( numberOfParticlesPerCell );
            innerOverlaps[ i ].assignCells( zoneOriginIdx, zoneDimensions, localGridDimensionsWithOverlap ); //FIXME: Why there is no localGridDimensionWithOverlap??????
         }
      }
   }

   void
   setDistributedGridParameters( const IndexVectorType& globalGridSize,
                                 const PointType& globalGridOrigin,
                                 const IndexVectorType& localGridSize,
                                 const PointType& localGridOrigin,
                                 const GlobalIndexType& numberOfOverlapsLayers,
                                 const RealType& searchRadius,
                                 const SubdomainCoordinates& domainDecomposition,
                                 MPI::Comm& comm )
   {
      this->communicator = comm;

      std::cout << "rank: <" << TNL::MPI::GetRank() << " globalGridSize: " << globalGridSize << \
                                                       " globalGridOrigin: " << globalGridOrigin << \
                                                       " localGridSize: " << localGridSize << \
                                                       " localGridOrigin: " << localGridOrigin << \
                                                       " localGridOriginINCELLS: " << TNL::ceil( ( localGridOrigin - globalGridOrigin ) / searchRadius ) << \
                                                       " numberOfOverlapsLayers " << numberOfOverlapsLayers << \
                                                       " searchRadius " << searchRadius << std::endl;

      //TODO: Pass as globalGrid
      GridType globalGrid;
      globalGrid.setDimensions( globalGridSize );
      globalGrid.setDomain( globalGridOrigin, globalGridSize );
      const PointType spaceStepsVector = searchRadius;
      globalGrid.setSpaceSteps( spaceStepsVector );


      //FIXME: Ugly workaround
      const IndexVectorType domainDecompositionVect = { domainDecomposition[ 0 ], domainDecomposition[ 1 ], 1 };
      //distributedGrid.setDomainDecomposition( domainDecomposition );
      distributedGrid.setDomainDecomposition( domainDecompositionVect );
      distributedGrid.setGlobalGrid( globalGrid );

      typename DistributedGridType::SubdomainOverlapsType lowerOverlap, upperOverlap;
      Meshes::DistributedMeshes::SubdomainOverlapsGetter< GridType >::getOverlaps( &distributedGrid, lowerOverlap, upperOverlap, 1 );
      distributedGrid.setOverlaps( lowerOverlap, upperOverlap );


      // Since the global is not distributed unifromly, we need to update parameters of local grid
      distributedGrid.localGrid.setOrigin( localGridOrigin );
      distributedGrid.localGrid.setDimensions( localGridSize );
      distributedGrid.localGrid.setSpaceSteps( distributedGrid.globalGrid.getSpaceSteps() );

      //NOTE: This is probably not necessay unless its used in initializeInnerOverlaps
      using CoordinatesType = typename DistributedGridType::CoordinatesType;
      CoordinatesType interiorBegin = this->distributedGrid.lowerOverlap;
      CoordinatesType interiorEnd = distributedGrid.localGrid.getDimensions() - this->distributedGrid.upperOverlap;
      const int* neighbors = distributedGrid.getNeighbors();
      if( neighbors[ ZzYzXm ] == -1 )
         interiorBegin[ 0 ] += 1;
      if( neighbors[ ZzYzXp ] == -1 )
         interiorEnd[ 0 ] -= 1;
      if( ZzYmXz < distributedGrid.getNeighborsCount() && neighbors[ ZzYmXz ] == -1 )
         interiorBegin[ 1 ] += 1;
      if( ZzYpXz < distributedGrid.getNeighborsCount() && neighbors[ ZzYpXz ] == -1 )
         interiorEnd[ 1 ] -= 1;
      if( ZmYzXz < distributedGrid.getNeighborsCount() && neighbors[ ZmYzXz ] == -1 )
         interiorBegin[ 2 ] += 1;
      if( ZpYzXz < distributedGrid.getNeighborsCount() && neighbors[ ZpYzXz ] == -1 )
         interiorEnd[ 2 ] -= 1;
      distributedGrid.localGrid.setInteriorBegin( interiorBegin );
      distributedGrid.localGrid.setInteriorEnd( interiorEnd );

      //Initialize inner particle zones to collect particles
      //TODO: Define inner overlaps with given size directly.
      innerOverlaps.resize( getNeighborsCount() );
      if constexpr ( ParticleSystem::spaceDimension == 2 )
         initializeInnerOverlaps( numberOfOverlapsLayers );
      if constexpr ( ParticleSystem::spaceDimension == 3 )
         initializeInnerOverlaps3D( numberOfOverlapsLayers );

      //initialize load balancing measures
      subdomainsParticlesCount = 0;
      subdomainsCompTime = 0.;
   }

   void
   updateDistriutedGridParameters( const IndexVectorType& updatedGridDimensions, const PointType& updatedGridOrigin, const GlobalIndexType& numberOfOverlapsLayers, const RealType& searchRadius )
   {
      distributedGrid.localGrid.setOrigin( updatedGridOrigin );
      distributedGrid.localGrid.setDimensions( updatedGridDimensions );
      //there is some wierd grid size, but search radius is not necessay here
      const PointType spaceStepsVector = searchRadius;
      distributedGrid.localGrid.setSpaceSteps( spaceStepsVector );

      if constexpr ( ParticleSystem::spaceDimension == 2 )
         initializeInnerOverlaps( numberOfOverlapsLayers );
      if constexpr ( ParticleSystem::spaceDimension == 3 )
         initializeInnerOverlaps3D( numberOfOverlapsLayers );
   }

   //collect particles to innerOverlaps
   template< typename ParticlePointer >
   void
   collectParticlesInInnerOverlaps( ParticlePointer& particles )
   {
      //NOTE: This was original idea, but the overlap size is much smaller.
      const int* neighbors = this->getDistributedGrid().getNeighbors();
      for( int i = 0; i < this->getDistributedGrid().getNeighborsCount(); i++ ) {
         //TODO: We shoud limit ourselves only to filled zones to save the call time
         if( neighbors[ i ] != -1 ){
            //innerOverlaps[ i ].resetParticles(); //this is not required
            innerOverlaps[ i ].updateParticlesInZone( particles );
            if( TNL::MPI::GetRank() == 2 )
            std::cout << "NUMBER OF PARTICLES IN ZONE: " << innerOverlaps[ i ].getNumberOfParticles() << " NUMBER OF CELLS: " << innerOverlaps[ i ].getNumberOfCells() << std::endl;
            //std::cout << innerOverlaps[ i ].getCellsInZone() << std::endl;
         }
      }
   }

   std::pair< IndexVectorType, IndexVectorType >
   loadBalancingDomainAdjustment()
   {
      IndexVectorType gridDimensionsAdjustment = 0;
      IndexVectorType gridOriginAdjustment = 0.;
      const int* neighbors = this->getDistributedGrid().getNeighbors();
      std::cout << "(RANK: " << TNL::MPI::GetRank()  << ") subdomainsParticlesCount = " <<  subdomainsParticlesCount << ", subdomainParticlesCount: " << subdomainParticlesCount  << std::endl;

      if( neighbors[ ZzYzXm ] != -1 ){
         const GlobalIndexType particlesCountDifference = subdomainParticlesCount - subdomainsParticlesCount[ ZzYzXm ];
         if( particlesCountDifference > this->particlesCountResizeTrashold  ){
            //pCD > pCRT -> interface mm => gridDimension.x++, gridOrigin--
            gridDimensionsAdjustment[ 0 ]--;
            gridOriginAdjustment[ 0 ]++;
            std::cout << "(RANK: " << TNL::MPI::GetRank()  << ") Balancing option: ( subdomainParticlesCount - subdomainsParticlesCount[ ZzYzXm ] ) = " << particlesCountDifference << " > pCRT."  << std::endl;
         }
         if( particlesCountDifference < ( ( -1 ) * this->particlesCountResizeTrashold ) ){
            //pCD > pCRT -> interface mm => gridDimension.x--, gridOrigin++
            gridDimensionsAdjustment[ 0 ]++;
            gridOriginAdjustment[ 0 ]--;
            std::cout << "(RANK: " << TNL::MPI::GetRank()  << ") Balancing option: ( subdomainParticlesCount - subdomainsParticlesCount[ ZzYzXm ] ) = " << particlesCountDifference << " < ( -1 ) * pCRT."  << std::endl;
         }

      }

      if( neighbors[ ZzYzXp ] != -1 ){
         const GlobalIndexType particlesCountDifference = subdomainParticlesCount - subdomainsParticlesCount[ ZzYzXp ];
         if( particlesCountDifference > this->particlesCountResizeTrashold  ){
            //pCD > pCRT -> interface pp => gridDimensions.x++, gridOring unchanged
            gridDimensionsAdjustment[ 0 ]--;
            std::cout << "(RANK: " << TNL::MPI::GetRank()  << ") Balancing option: ( subdomainParticlesCount - subdomainsParticlesCount[ ZzYzXp ] ) = " <<  particlesCountDifference << " > pCRT."  << std::endl;
         }
         if( particlesCountDifference < ( ( -1 ) * this->particlesCountResizeTrashold ) ){
            //pCD > pCRT -> interface pm => gridDimensions.x--, gridOrigin unchanged
            gridDimensionsAdjustment[ 0 ]++;
            std::cout << "(RANK: " << TNL::MPI::GetRank()  << ") Balancing option: ( subdomainParticlesCount - subdomainsParticlesCount[ ZzYzXp ] ) = " <<  particlesCountDifference << " < ( -1 ) * pCRT."  << std::endl;
         }
      }

      return std::make_pair( gridDimensionsAdjustment, gridOriginAdjustment );
   }



   //collect particles to innerOverlaps
   template< typename ParticlePointer >
   void
   DEBUG_printParticlesInOverlap( ParticlePointer& particles )
   {
      auto points_view = particles->getPoints().getConstView();
      auto cellIndex_view = particles->getParticleCellIndices().getConstView();

      using CellIndexer = typename ParticleSystem::CellIndexer;

      //NOTE: This was original idea, but the overlap size is much smaller.
      const int* neighbors = this->getDistributedGrid().getNeighbors();
      for( int i = 0; i < this->getDistributedGrid().getNeighborsCount(); i++ ) {
         //TODO: We shoud limit ourselves only to filled zones to save the call time
         if( neighbors[ i ] != -1 ){
            const auto ghostZoneView = innerOverlaps[ i ].getParticlesInZone().getConstView();
            const GlobalIndexType numberOfParticlesInZone = innerOverlaps[ i ].getNumberOfParticles();

            const RealType searchRadius = particles->getSearchRadius();
            const PointType gridOrigin = particles->getGridOrigin();
            const PointType gridDimension = particles->getGridDimensions();

            const float scaleFactor = 1.f ;

            auto init = [=] __cuda_callable__ ( GlobalIndexType i ) mutable
            {
               const GlobalIndexType p = ghostZoneView[ i ];
               //printf( "[ %f, %f, (%d), (%d), <%f>, <%f>, {%f} ]",
               printf( "[ %.12f, %f, (%d), (%d), <%f>, <%f>, {%.12f} ]",
                     points_view[ p ][ 0 ],
                     points_view[ p ][ 1 ],
                     cellIndex_view[ p ],
                     CellIndexer::EvaluateCellIndex( points_view[ p ], gridOrigin, gridDimension, searchRadius ),
                     ( points_view[ p ][ 0 ] * scaleFactor - gridOrigin[ 0 ] * scaleFactor ) / ( searchRadius * scaleFactor ) ,
                     TNL::floor( ( points_view[ p ][ 0 ] - gridOrigin[ 0 ] ) / searchRadius ),
                     points_view[ p ][ 0 ] - gridOrigin[ 0 ] );
            };
            Algorithms::parallelFor< DeviceType >( 0, numberOfParticlesInZone, init );
         }

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

   ////TODO: This is desinged currently only for 1D decomposition
   //void
   //collectParticlesToTransfer()
   //{
   //   //2D version
   //   //const int* neighbors = this->getDistributedGrid().getNeighbors();
   //   //for( int i = 0; i < this->getNeighborsCount(); i++ ){
   //   //   if( neighbors[ i ] != -1 ){

   //   //      auto particlesInZone_view = innerOverlaps[ i ].getParticlesInZone().getView();
   //   //      //compa


   //   //   }
   //   //}

   //   //1D version
   //   if( distributedGrid.isThereNeighbor( Directions::template getXYZ< 2 >( ZzYzXm ) ) ){

   //      auto particlesInZone_view = innerOverlaps[ ZzYzXm ].getParticlesInZone().getView();

   //   }
   //   if( distributedGrid.isThereNeighbor( Directions::template getXYZ< 2 >( ZzYzXp ) ) ){
   //
   //   }

   //}

   //linearize indicse inside innerOverlaps
   void
   linearize()
   {
      // 6 5 7
      // 0 * 1
      // 2 3 4
      const int* neighbors = this->distributedGrid->getNeighbors();

      /*
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
      */

      //copy all to innerOverlapsLinearized
      for( int i = 0; i < this->distributedGrid->getNeighborsCount(); i++ ) {
         if( neighbors[ i ] != -1 ){
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
   }

   void
   writeProlog( TNL::Logger& logger ) //const noexcept
   {
      logger.writeParameter( "Distributed particles:", "" );
      distributedGrid.writeProlog( logger );

      //debug info from local subdomains
      logger.writeParameter( "Lower overlap:", distributedGrid.getLowerOverlap() );
      logger.writeParameter( "Upper overlap:", distributedGrid.getUpperOverlap() );
      logger.writeParameter( "Local grid origin:", distributedGrid.localGrid.getOrigin() );
      logger.writeParameter( "Local grid end:", distributedGrid.localGrid.getOrigin() + distributedGrid.localGrid.getSpaceSteps() * distributedGrid.localGrid.getDimensions() );
      logger.writeParameter( "Local grid size:", distributedGrid.localGrid.getSpaceSteps() * distributedGrid.localGrid.getDimensions() );
      logger.writeParameter( "Local grid space steps:", distributedGrid.localGrid.getSpaceSteps() );
      logger.writeParameter( "Local grid dimensions:", distributedGrid.localGrid.getDimensions() );
      logger.writeParameter( "Local grid interior begin:", distributedGrid.localGrid.getInteriorBegin() );
      logger.writeParameter( "Local grid interior end:", distributedGrid.localGrid.getInteriorEnd() );
      //logger.writeParameter( "Total neighbors count:", distributedGrid.getNeighborsCount() );
      //zones
      for( int i = 0; i < innerOverlaps.getSize(); i++ )
         innerOverlaps[ i ].writeProlog( logger );



   }

protected:

   ParticleSystem localParticles;
   DistributedGridType distributedGrid;
   //Containers::Array< GlobalIndexType, Devices::Host, int > innerOverlaps; //TODO: What was this idee?

   Containers::Array< ParticleZoneType, Devices::Host, int > innerOverlaps; //TODO: What was this idee?

   GlobalIndexType numberOfParticlesInOverlaps;

   Containers::StaticArray< DistributedGridType::getNeighborsCount(), int > innerOverlapsOffests;
   IndexArrayType innerOverlapsLinearized;

   Containers::StaticArray< DistributedGridType::getNeighborsCount(), int > subdomainsParticlesCount;
   GlobalIndexType subdomainParticlesCount;
   int particlesCountResizeTrashold;

   Containers::StaticArray< DistributedGridType::getNeighborsCount(), float > subdomainsCompTime;
   RealType subdomainCompTime;
   RealType computationalTimeResizeTrashold;

   MPI::Comm communicator = MPI_COMM_WORLD;


};

}  //namespace ParticleSystem
}  //namespace TNL

