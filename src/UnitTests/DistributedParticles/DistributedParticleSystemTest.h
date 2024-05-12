#pragma once

#include "../main_mpi.h"

#include <complex>
#include <gtest/gtest.h>
#include "gtest/gtest.h"
#include "gmock/gmock.h" //test vectors

#include <Particles/ParticlesLinkedList.h>
#include <Particles/DistributedParticles.h>
#include <Particles/DistributedParticlesSynchronizer.h>
#include <type_traits>


using namespace TNL;
using namespace ParticleSystem;

template< typename Device >
class Particles2DConfig
{
   public:
   using DeviceType = Device;

   using GlobalIndexType = int;
   using LocalIndexType = int;
   using CellIndexType = int;
   using RealType = float;

   static constexpr int spaceDimension = 2;

   using UseWithDomainDecomposition = std::true_type;
   using CoordinatesType = Containers::StaticVector< spaceDimension, int >;
   using CellIndexerType = SimpleCellIndex< spaceDimension, Particles2DConfig, std::index_sequence< 1, 0 > >;
};

template< typename Device >
class Particles2DSetup
{
   public:
   using ParticlesConfig = Particles2DConfig< Device >;
   using ParticlesTraitsType = ParticlesTraits< ParticlesConfig, Device >;
   using IndexVectorType = typename ParticlesTraitsType::IndexVectorType;
   using PointType = typename ParticlesTraitsType::PointType;

   const IndexVectorType numberOfSubdomains = { 3 , 1 };

   const int numberOfParticles = 20;
   const int numberOfAllocatedParticles = 28;

   const float searchRadius = 0.5;
   const PointType gridOrigin = { -searchRadius, -searchRadius };
   const IndexVectorType overlapSize = 1;
   const PointType overlapWidth = searchRadius;

   const IndexVectorType gridDimension = { 6, 5 };
   const int numberOfGridCells = gridDimension[ 0 ] * gridDimension[ 1 ];
};

template< typename Device >
struct SubdomainsParametersX3Y1
{
   using ParticlesConfig = Particles2DConfig< Device >;
   using ParticlesTraitsType = ParticlesTraits< ParticlesConfig, Device >;
   using IndexVectorType = typename ParticlesTraitsType::IndexVectorType;
   using PointType = typename ParticlesTraitsType::PointType;

   const IndexVectorType subdomains = { 3 , 1 };
   const int subdomainsCount = subdomains[ 0 ] * subdomains[ 1 ];
   Containers::StaticVector< 3, int > numberOfParticles = { 20, 20, 20 };
   Containers::StaticVector< 3, int > numberOfAllocatedParticles = { 30, 30, 30 };
   const float searchRadius = 0.5;

   const PointType subdomain1Origin = { -searchRadius, -searchRadius };
   const PointType subdomain2Origin = { 2 * searchRadius, 2 * searchRadius };
   const PointType subdomain3Origin = { 7 * searchRadius, 7 * searchRadius };
   Containers::StaticVector< 3, PointType > gridOrigin = { subdomain1Origin, subdomain2Origin, subdomain3Origin };
   const IndexVectorType subdomain1Dimension = { 3, 3 };
   const IndexVectorType subdomain2Dimension = { 5, 5 };
   const IndexVectorType subdomain3Dimension = { 4, 4 };
   Containers::StaticVector< 3, IndexVectorType > gridDimension = { subdomain1Dimension,
                                                                    subdomain2Dimension,
                                                                    subdomain3Dimension };
   const IndexVectorType overlapSize = 1;
   const PointType overlapWidth = searchRadius;

};

template< typename ParticlePointer >
bool assignPoints2D( ParticlePointer& particles )
{
   //setup points:
   auto points = particles->getPoints().getView();

   //[ 1, 1 ]
   points.setElement( 10, { 0.32, 0.2 } );
   //[ 2, 1 ]
   points.setElement( 4, { 0.78, 0.23 } );
   points.setElement( 16, { 0.62, 0.39 } );
   //[ 3, 1 ]
   points.setElement( 5, { 1.22, 0.25 } );
   //[ 4, 1 ]
   points.setElement( 11, { 1.6, 0.1 } );

   //[ 1, 2 ]
   points.setElement( 0, { 0.17, 0.6 } );
   points.setElement( 18, { 0.3, 0.82 } );
   //[ 2, 2 ]
   points.setElement( 19, { 0.65, 0.7 } );
   //[ 3, 2 ]
   points.setElement( 6, { 1.25, 0.58 } );
   points.setElement( 9, { 1.2, 0.92 } );
   //[ 4, 2 ]
   points.setElement( 2, { 1.67, 0.69 } );
   points.setElement( 14, { 1.9, 0.58 } );
   points.setElement( 15, { 1.8, 0.89 } );

   //[ 1, 3 ]
   points.setElement( 1, { 0.18, 1.3 } );
   points.setElement( 12, { 0.4, 1.15 } );
   //[ 2, 3 ]
   points.setElement( 7, { 0.53, 1.33 } );
   points.setElement( 8, { 0.8, 1.3 } );
   points.setElement( 17, { 0.7, 1.18 } );
   //[ 3, 3 ]
   points.setElement( 13, { 1.3, 1.21 } );
   //[ 4, 3 ]
   points.setElement( 3, { 1.7, 1.33 } );

  return true;
}



TEST( DistributedParticles1DSplittingTest, DistributedParticlesInitialization )
{
   using Device = TNL::Devices::Host;
   using ParticlesSetup = Particles2DSetup< Device >;
   using SubdomainsSetup = SubdomainsParametersX3Y1< Device >;

   using Particles = TNL::ParticleSystem::ParticlesLinkedList< ParticlesSetup::ParticlesConfig, Device >;
   using ParticlesPointer = typename Pointers::SharedPointer< Particles, Device >;
   using DistributedParticles = TNL::ParticleSystem::DistributedParticleSystem< Particles >;
   using DistributedParticlesPointer = typename Pointers ::SharedPointer< DistributedParticles, Device >;
   using Synchronizer = TNL::ParticleSystem::DistributedParticlesSynchronizer< DistributedParticles >;

   ParticlesSetup setup;
   SubdomainsSetup subdomainsSetup;

   ParticlesPointer particles;
   DistributedParticlesPointer distributedParticles;
   Synchronizer synchronizer;
   MPI::Comm communicator = MPI_COMM_WORLD;

   const int rank = TNL::MPI::GetRank();

   // setup particles
   particles->setSize( subdomainsSetup.numberOfAllocatedParticles[ rank ] );
   particles->setSearchRadius( subdomainsSetup.searchRadius );
   particles->setGridSize( subdomainsSetup.gridDimension[ rank ] + subdomainsSetup.overlapSize );
   particles->setGridOrigin( subdomainsSetup.gridOrigin[ rank ] - subdomainsSetup.overlapWidth );
   particles->setNumberOfParticles( subdomainsSetup.numberOfParticles[ rank ] );
   particles->setFirstActiveParticle( 0 );
   particles->setLastActiveParticle( subdomainsSetup.numberOfParticles[ rank ] - 1 );

   // setup distributed particles
   std::ofstream logFile( "testLog" );
   TNL::Logger logger( 0, logFile );
   distributedParticles->setDistributedGridParameters( subdomainsSetup.gridDimension[ rank ] + subdomainsSetup.overlapSize,
                                                       subdomainsSetup.gridOrigin[ rank ] - subdomainsSetup.overlapWidth,
                                                       subdomainsSetup.gridDimension[ rank ],
                                                       subdomainsSetup.gridOrigin[ rank ],
                                                       1,
                                                       subdomainsSetup.searchRadius,
                                                       subdomainsSetup.subdomains,
                                                       communicator );

   //assgn particles
   ASSERT_TRUE( assignPoints2D( particles ) );

   //test number of points
   EXPECT_EQ( particles->getNumberOfParticles(), 20 );
   EXPECT_EQ( particles->getNumberOfAllocatedParticles(), 30 );

   using Grid = typename DistributedParticles::GridType;
   using Point = typename Grid::PointType;
   using CoordinatesType = typename Grid::CoordinatesType;

   Grid localGrid = distributedParticles->getDistributedGrid().getLocalMesh();
   if( rank == 0 ){
      Point s1_localGridOrigin = { -0.5, -0.5 };
      EXPECT_EQ( localGrid.getOrigin(), s1_localGridOrigin );
      CoordinatesType s1_localGridDimension = { 3, 3 };
      EXPECT_EQ( localGrid.getDimensions(), s1_localGridDimension );
   }
   if( rank == 1 ){
      Point s1_localGridOrigin = { 1.0, 1.0 };
      EXPECT_EQ( localGrid.getOrigin(), s1_localGridOrigin );
      CoordinatesType s1_localGridDimension = { 5, 5 };
      EXPECT_EQ( localGrid.getDimensions(), s1_localGridDimension );
   }
   if( rank == 2 ){
      Point s1_localGridOrigin = { 3.5, 3.5 };
      EXPECT_EQ( localGrid.getOrigin(), s1_localGridOrigin );
      CoordinatesType s1_localGridDimension = { 4, 4 };
      EXPECT_EQ( localGrid.getDimensions(), s1_localGridDimension );
   }

}

