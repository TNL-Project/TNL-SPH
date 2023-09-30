#pragma once

#include <gtest/gtest.h>
#include <Particles/ParticlesLinkedListFloating.h>
#include <TNL/Algorithms/Segments/Ellpack.h> //TODO: Reconsider.
#include "gtest/gtest.h"
#include "gmock/gmock.h" //test vectors


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

   using CoordinatesType = Containers::StaticVector< spaceDimension, int >;
   using CellIndexerType = SimpleCellIndex< spaceDimension, Particles2DConfig, std::index_sequence< 0, 1 > >;
   using NeighborListType = typename Algorithms::Segments::Ellpack< Device, int >; //deprecated
};

template< typename Device >
class Particles2DSetup
{
   public:
   using ParticlesConfig = Particles2DConfig< Device >;
   using ParticlesTraitsType = ParticlesTraits< ParticlesConfig, Device >;
   using IndexVectorType = typename ParticlesTraitsType::IndexVectorType;
   using PointType = typename ParticlesTraitsType::PointType;

   const int numberOfParticles = 20;
   const int numberOfAllocatedParticles = 28;

   const float searchRadius = 0.5;
   const int gridXsize = 6;
   const int gridYsize = 5;
   const PointType gridOrigin = { -searchRadius, -searchRadius };

   const IndexVectorType gridSize = { gridXsize, gridYsize };
   const int numberOfGridCells = gridXsize * gridYsize;
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

TEST( SearchForNeighbors2DTest, ParticlesPropertiesHost )
{
   using Device = TNL::Devices::Host;
   using ParticlesSetup = Particles2DSetup< Device >;
   using Particles = TNL::ParticleSystem::ParticlesLinkedList< ParticlesSetup::ParticlesConfig, Device >;
   using ParticlesPointer = typename Pointers::SharedPointer< Particles, Device >;

   ParticlesSetup setup;
   ParticlesPointer particles( setup.numberOfParticles,
                               setup.numberOfAllocatedParticles,
                               setup.searchRadius,
                               setup.numberOfGridCells );
   particles->setGridSize( setup.gridSize );
   particles->setGridOrigin( setup.gridOrigin );

   //assgn particles
   ASSERT_TRUE( assignPoints2D( particles ) );

   //test number of points
   EXPECT_EQ( particles->getNumberOfParticles(), 20 );
   EXPECT_EQ( particles->getNumberOfAllocatedParticles(), 28 );

}

TEST( SearchForNeighbors2DTest, ParticlesPropertiesCuda )
{
   using Device = TNL::Devices::Cuda;
   using ParticlesSetup = Particles2DSetup< Device >;
   using Particles = TNL::ParticleSystem::ParticlesLinkedList< ParticlesSetup::ParticlesConfig, Device >;
   using ParticlesPointer = typename Pointers::SharedPointer< Particles, Device >;

   ParticlesSetup setup;
   ParticlesPointer particles( setup.numberOfParticles,
                               setup.numberOfAllocatedParticles,
                               setup.searchRadius,
                               setup.numberOfGridCells );
   particles->setGridSize( setup.gridSize );
   particles->setGridOrigin( setup.gridOrigin );

   //assgn particles
   ASSERT_TRUE( assignPoints2D( particles ) );

   //test number of points
   EXPECT_EQ( particles->getNumberOfParticles(), 20 );
   EXPECT_EQ( particles->getNumberOfAllocatedParticles(), 28 );

}

TEST( SearchForNeighbors2DTest, ComputeParticleCellIndicesCuda )
{
   using Device = TNL::Devices::Cuda;
   using ParticlesSetup = Particles2DSetup< Device >;
   using Particles = TNL::ParticleSystem::ParticlesLinkedList< ParticlesSetup::ParticlesConfig, Device >;
   using ParticlesPointer = typename Pointers::SharedPointer< Particles, Device >;

   ParticlesSetup setup;
   ParticlesPointer particles( setup.numberOfParticles,
                               setup.numberOfAllocatedParticles,
                               setup.searchRadius,
                               setup.numberOfGridCells );
   particles->setGridSize( setup.gridSize );
   particles->setGridOrigin( setup.gridOrigin );

   //assgn particles
   ASSERT_TRUE( assignPoints2D( particles ) );

   //test number of points
   particles->computeParticleCellIndices();

   const auto cellIndices = particles->getParticleCellIndices().getConstView();

   //[ 1, 1 ]
   EXPECT_EQ( cellIndices.getElement( 10 ), 7 );
   //[ 2, 1 ]
   EXPECT_EQ( cellIndices.getElement( 4 ), 8 );
   EXPECT_EQ( cellIndices.getElement( 16 ), 8 );
   //[ 3, 1 ]
   EXPECT_EQ( cellIndices.getElement( 5 ), 9 );
   //[ 4, 1 ]
   EXPECT_EQ( cellIndices.getElement( 11 ), 10 );

   //[ 1, 2 ]
   EXPECT_EQ( cellIndices.getElement( 0 ), 13 );
   EXPECT_EQ( cellIndices.getElement( 18 ), 13 );
   //[ 2, 2 ]
   EXPECT_EQ( cellIndices.getElement( 19 ), 14 );
   //[ 3, 2 ]
   EXPECT_EQ( cellIndices.getElement( 6 ), 15 );
   EXPECT_EQ( cellIndices.getElement( 9 ), 15 );
   //[ 4, 2 ]
   EXPECT_EQ( cellIndices.getElement( 2 ), 16 );
   EXPECT_EQ( cellIndices.getElement( 14 ), 16 );
   EXPECT_EQ( cellIndices.getElement( 15 ), 16 );

   //[ 1, 3 ]
   EXPECT_EQ( cellIndices.getElement( 1 ), 19 );
   EXPECT_EQ( cellIndices.getElement( 12 ), 19 );
   //[ 2, 3 ]
   EXPECT_EQ( cellIndices.getElement( 17 ), 20 );
   EXPECT_EQ( cellIndices.getElement( 7 ), 20 );
   EXPECT_EQ( cellIndices.getElement( 8 ), 20 );
   //[ 3, 3 ]
   EXPECT_EQ( cellIndices.getElement( 13 ), 21 );
   //[ 4, 3 ]
   EXPECT_EQ( cellIndices.getElement( 3 ), 22 );

}

TEST( SearchForNeighbors2DTest, SortParticlesCuda )
{
   using Device = TNL::Devices::Cuda;
   using ParticlesSetup = Particles2DSetup< Device >;
   using Particles = TNL::ParticleSystem::ParticlesLinkedList< ParticlesSetup::ParticlesConfig, Device >;
   using ParticlesPointer = typename Pointers::SharedPointer< Particles, Device >;
   using Point = typename Particles::PointType;

   ParticlesSetup setup;
   ParticlesPointer particles( setup.numberOfParticles,
                               setup.numberOfAllocatedParticles,
                               setup.searchRadius,
                               setup.numberOfGridCells );
   particles->setGridSize( setup.gridSize );
   particles->setGridOrigin( setup.gridOrigin );

   //assgn particles
   ASSERT_TRUE( assignPoints2D( particles ) );

   //test number of points
   particles->computeParticleCellIndices();
   particles->sortParticles();

   const auto points = particles->getPoints().getConstView();
   const auto cellIndices = particles->getParticleCellIndices().getConstView();
   const auto permutations = particles->getSortPermutations()->getConstView(); //FIXME: Is necessary to have map as pointer?

   //[ 1, 1 ]
   EXPECT_EQ( cellIndices.getElement( 0 ), 7 );
   Point p0 = { 0.32, 0.2 };
   EXPECT_EQ( points.getElement( 0 ), p0 );
   EXPECT_EQ( permutations.getElement( 0 ), 10 );
   //[ 2, 1 ]
   EXPECT_EQ( cellIndices.getElement( 1 ), 8 );
   Point p1 = { 0.78, 0.23 };
   EXPECT_EQ( points.getElement( 1 ), p1 );
   EXPECT_EQ( permutations.getElement( 1 ), 4 );
   EXPECT_EQ( cellIndices.getElement( 2 ), 8 );
   Point p2 = { 0.62, 0.39 };
   EXPECT_EQ( points.getElement( 2 ), p2 );
   EXPECT_EQ( permutations.getElement( 2 ), 16 );
   //[ 3, 1 ]
   EXPECT_EQ( cellIndices.getElement( 3 ), 9 );
   Point p3 = { 1.22, 0.25 };
   EXPECT_EQ( points.getElement( 3 ), p3 );
   EXPECT_EQ( permutations.getElement( 3 ), 5 );
   //[ 4, 1 ]
   EXPECT_EQ( cellIndices.getElement( 4 ), 10 );
   Point p4 = { 1.6, 0.1 };
   EXPECT_EQ( points.getElement( 4 ), p4 );
   EXPECT_EQ( permutations.getElement( 4 ), 11 );

   //[ 1, 2 ]
   EXPECT_EQ( cellIndices.getElement( 5 ), 13 );
   Point p5 = { 0.17, 0.6 };
   EXPECT_EQ( points.getElement( 5 ), p5 );
   EXPECT_EQ( permutations.getElement( 5 ), 0 );
   EXPECT_EQ( cellIndices.getElement( 6 ), 13 );
   Point p6 = { 0.3, 0.82 };
   EXPECT_EQ( points.getElement( 6 ), p6 );
   EXPECT_EQ( permutations.getElement( 6 ), 18 );
   //[ 2, 2 ]
   EXPECT_EQ( cellIndices.getElement( 7 ), 14 );
   Point p7 = { 0.65, 0.7 };
   EXPECT_EQ( points.getElement( 7 ), p7 );
   EXPECT_EQ( permutations.getElement( 7 ), 19 );
   //[ 3, 2 ]
   EXPECT_EQ( cellIndices.getElement( 8 ), 15 );
   Point p8 = { 1.25, 0.58 };
   EXPECT_EQ( points.getElement( 8 ), p8 );
   EXPECT_EQ( permutations.getElement( 8 ), 6 );
   EXPECT_EQ( cellIndices.getElement( 9 ), 15 );
   Point p9 = { 1.2, 0.92 };
   EXPECT_EQ( points.getElement( 9 ), p9 );
   EXPECT_EQ( permutations.getElement( 9 ), 9 );
   //[ 4, 2 ]
   EXPECT_EQ( cellIndices.getElement( 10 ), 16 );
   Point p10 = { 1.67, 0.69 };
   EXPECT_EQ( points.getElement( 10 ), p10 );
   EXPECT_EQ( permutations.getElement( 10 ), 2 );
   EXPECT_EQ( cellIndices.getElement( 11 ), 16 );
   Point p11 = { 1.9, 0.58 };
   EXPECT_EQ( points.getElement( 11 ), p11 );
   EXPECT_EQ( permutations.getElement( 11 ), 14 );
   EXPECT_EQ( cellIndices.getElement( 12 ), 16 );
   Point p12 = { 1.8, 0.89 };
   EXPECT_EQ( points.getElement( 12 ), p12 );
   EXPECT_EQ( permutations.getElement( 12 ), 15 );

   //[ 1, 3 ]
   EXPECT_EQ( cellIndices.getElement( 13 ), 19 );
   Point p13 = { 0.18, 1.3 };
   EXPECT_EQ( points.getElement( 13 ), p13 );
   EXPECT_EQ( permutations.getElement( 13 ), 1 );
   EXPECT_EQ( cellIndices.getElement( 14 ), 19 );
   Point p14 = { 0.4, 1.15 };
   EXPECT_EQ( points.getElement( 14 ), p14 );
   EXPECT_EQ( permutations.getElement( 14 ), 12 );
   //[ 2, 3 ]
   EXPECT_EQ( cellIndices.getElement( 15 ), 20 );
   Point p15 = { 0.53, 1.33 };
   EXPECT_EQ( points.getElement( 15 ), p15 );
   EXPECT_EQ( permutations.getElement( 15 ), 7 );
   EXPECT_EQ( cellIndices.getElement( 16 ), 20 );
   Point p16 = { 0.8, 1.3 };
   EXPECT_EQ( points.getElement( 16 ), p16 );
   EXPECT_EQ( permutations.getElement( 16 ), 8 );
   EXPECT_EQ( cellIndices.getElement( 17 ), 20 );
   Point p17 = { 0.7, 1.18 };
   EXPECT_EQ( points.getElement( 17 ), p17 );
   EXPECT_EQ( permutations.getElement( 17 ), 17 );
   //[ 3, 3 ]
   EXPECT_EQ( cellIndices.getElement( 18 ), 21 );
   Point p18 = { 1.3, 1.21 };
   EXPECT_EQ( points.getElement( 18 ), p18 );
   EXPECT_EQ( permutations.getElement( 18 ), 13 );
   //[ 4, 3 ]
   EXPECT_EQ( cellIndices.getElement( 19 ), 22 );
   Point p19 = { 1.7, 1.33 };
   EXPECT_EQ( points.getElement( 19 ), p19 );
   EXPECT_EQ( permutations.getElement( 19 ), 3 );

}

