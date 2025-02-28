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

   void
   //setCommunicator( const MPI::Comm& communicator )
   setCommunicator( MPI::Comm& communicator )
   {
      this->communicator = communicator;
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
   setDistributedGridParameters( const RealType& searchRadius,
                                 const IndexVectorType& globalGridDimension,
                                 const PointType& globalGridOrigin,
                                 const IndexVectorType& localGridDimensions,
                                 const IndexVectorType& localGridOriginCoords,
                                 const int& numberOfOverlapsLayers,
                                 const Containers::StaticVector< 2, int >& numberOfSubdomains )
   {
      // The topology of grid is handled by DistributedGrid.h class, which is set here.
      // setup global grid (NOTE: This is not nice due to distributed grid interface)
      GridType globalGrid;
      globalGrid.setDimensions( globalGridDimension );
      globalGrid.setDomain( globalGridOrigin, globalGridDimension );
      const PointType spaceStepsVector = searchRadius;
      globalGrid.setSpaceSteps( spaceStepsVector );

      // setup distributed grid
      const IndexVectorType domainDecompositionVect = { numberOfSubdomains[ 0 ], numberOfSubdomains[ 1 ], 1 }; //FIXME
      distributedGrid.setDomainDecomposition( domainDecompositionVect );
      distributedGrid.setGlobalGrid( globalGrid );
      typename DistributedGridType::SubdomainOverlapsType lowerOverlap, upperOverlap;
      Meshes::DistributedMeshes::SubdomainOverlapsGetter< GridType >::getOverlaps( &distributedGrid, lowerOverlap, upperOverlap, numberOfOverlapsLayers );
      distributedGrid.setOverlaps( lowerOverlap, upperOverlap );

      // since the global is not distributed unifromly, we need to update parameters of local grid
      const PointType localGridOrigin = globalGridOrigin + searchRadius * localGridOriginCoords;
      distributedGrid.localGrid.setOrigin( localGridOrigin );
      distributedGrid.localGrid.setDimensions( localGridDimensions );
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
      const int* neighbors = this->getDistributedGrid().getNeighbors();
      for( int i = 0; i < this->getDistributedGrid().getNeighborsCount(); i++ ) {
         if( neighbors[ i ] != -1 )
            innerOverlaps[ i ].updateParticlesInZone( particles );
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

   void
   writeProlog( TNL::Logger& logger ) //const noexcept
   {
      logger.writeParameter( "Distributed particles:", "" );
      distributedGrid.writeProlog( logger );
      logger.writeParameter( "Lower overlap:", distributedGrid.getLowerOverlap() );
      logger.writeParameter( "Upper overlap:", distributedGrid.getUpperOverlap() );
      logger.writeParameter( "Local grid origin:", distributedGrid.localGrid.getOrigin() );
      logger.writeParameter( "Local grid end:", distributedGrid.localGrid.getOrigin() + distributedGrid.localGrid.getSpaceSteps() * distributedGrid.localGrid.getDimensions() );
      logger.writeParameter( "Local grid size:", distributedGrid.localGrid.getSpaceSteps() * distributedGrid.localGrid.getDimensions() );
      logger.writeParameter( "Local grid space steps:", distributedGrid.localGrid.getSpaceSteps() );
      logger.writeParameter( "Local grid dimensions:", distributedGrid.localGrid.getDimensions() );
      logger.writeParameter( "Local grid interior begin:", distributedGrid.localGrid.getInteriorBegin() );
      logger.writeParameter( "Local grid interior end:", distributedGrid.localGrid.getInteriorEnd() );
      // print zones details
      //logger.writeParameter( "Total neighbors count:", distributedGrid.getNeighborsCount() );
      //for( int i = 0; i < innerOverlaps.getSize(); i++ )
      //   innerOverlaps[ i ].writeProlog( logger );
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

