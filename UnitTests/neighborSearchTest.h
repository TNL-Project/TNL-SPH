#pragma once

// #ifdef HAVE_GTEST
#include <gtest/gtest.h>

#include <TNL/Algorithms/Segments/Ellpack.h>

#include "../Particles/Particles.h"
#include "../Particles/neighbourSearch.h"

#include "gtest/gtest.h"

using namespace TNL;
using namespace ParticleSystem;

// build particle system

template< typename Device >
class ParticlesConfig
{
   public:

   using GlobalIndexType = int;
   using LocalIndexType = int;
   using CellIndexType = int;
   using RealType = float;

   using DeviceType = Device;
   using CellIndexerType = SimpleCellIndex<ParticlesConfig< DeviceType >, DeviceType>;

   static constexpr int spaceDimension = 2;
   static constexpr int maxOfNeigborsPerParticle = 9;

   static constexpr float searchRadius = 0.75;
   static constexpr int gridXsize = 2;
   static constexpr int gridYsize = 2;

   static constexpr int gridXbegin = 0;
   static constexpr int gridYbegin = 0;

   using NeighborListType = typename Algorithms::Segments::Ellpack< DeviceType, GlobalIndexType >;
};

template< typename Device >
class ParticlesConfigCentric
{
   public:

   using GlobalIndexType = int;
   using LocalIndexType = int;
   using CellIndexType = int;
   using RealType = float;

   using DeviceType = Device;
   using CellIndexerType = SimpleCellIndex<ParticlesConfigCentric< DeviceType >, DeviceType>;

   static constexpr int spaceDimension = 2;
   static constexpr int maxOfNeigborsPerParticle = 9;

   static constexpr float searchRadius = 0.75;
   static constexpr int gridXsize = 2;
   static constexpr int gridYsize = 2;


   static constexpr RealType gridXbegin = -0.5;
   static constexpr RealType gridYbegin = -0.5;

   using NeighborListType = typename Algorithms::Segments::Ellpack< DeviceType, GlobalIndexType >;
};

template< typename Device >
bool generate2DParticleSystem( Particles< ParticlesConfig< Device >, Device >& particles)
{

  //
  //particles = Particles< ParticlesConfig< Device >, Device >( 9, 0.75 );
  particles = Particles< ParticlesConfig< Device >, Device >( 9, 0.75 );

  //setup points
  particles.setPoint( 0, { 0., 0. } );
  particles.setPoint( 1, { 0., 0.5 } );
  particles.setPoint( 2, { 0., 1. } );

  particles.setPoint( 3, { 0.5, 0. } );
  particles.setPoint( 4, { 0.5, 0.5 } );
  particles.setPoint( 5, { 0.5, 1. } );

  particles.setPoint( 6, { 1., 0. } );
  particles.setPoint( 7, { 1., 0.5 } );
  particles.setPoint( 8, { 1., 1. } );


  return true;
}

template< typename Device >
bool generate2DParticleSystemCentric( Particles< ParticlesConfigCentric< Device >, Device >& particles)
{

  //
  //particles = Particles< ParticlesConfig< Device >, Device >( 9, 0.75 );
  particles = Particles< ParticlesConfigCentric< Device >, Device >( 9, 0.75 );

  //setup points
  particles.setPoint( 0, { -0.5, -0.5 } );
  particles.setPoint( 1, { -0.5, 0. } );
  particles.setPoint( 2, { -0.5, 0.5 } );

  particles.setPoint( 3, { 0., -0.5 } );
  particles.setPoint( 4, { 0., 0. } );
  particles.setPoint( 5, { 0., 0.5 } );

  particles.setPoint( 6, { 0.5, -0.5 } );
  particles.setPoint( 7, { 0.5, 0. } );
  particles.setPoint( 8, { 0.5, 0.5 } );


  return true;
}

TEST( GetTwoTest, Two ) {
    EXPECT_EQ( 2, 2 );
}

TEST( SearchForNeighborsTest, NeighborSearchHost )
{
   using ParticlesConfig =  ParticlesConfig< TNL::Devices::Host > ;
   using ParticlesHost = Particles< ParticlesConfig, Devices::Host >;
   using neighborSearchHost = NeighborSearch< ParticlesConfig, ParticlesHost >;

   using PointType = typename ParticlesHost::PointType;

   ParticlesHost particles;

   ASSERT_TRUE( generate2DParticleSystem( particles ) );

   //test number of points
   EXPECT_EQ( particles.getNumberOfParticles(), 9 );

   //test particle position
   PointType point0( 0., 0. ),  point1( 0., 0.5 ),  point2( 0., 1. ),
             point3( 0.5, 0. ), point4( 0.5, 0.5 ), point5( 0.5, 1. ),
             point6( 1., 0. ),  point7( 1., 0.5 ),  point8( 1., 1. );

   EXPECT_EQ( particles.getPoint( 0 ), point0 );
   EXPECT_EQ( particles.getPoint( 1 ), point1 );
   EXPECT_EQ( particles.getPoint( 2 ), point2 );
   EXPECT_EQ( particles.getPoint( 3 ), point3 );
   EXPECT_EQ( particles.getPoint( 4 ), point4 );
   EXPECT_EQ( particles.getPoint( 5 ), point5 );
   EXPECT_EQ( particles.getPoint( 6 ), point6 );
   EXPECT_EQ( particles.getPoint( 7 ), point7 );
   EXPECT_EQ( particles.getPoint( 8 ), point8 );

   //test grid index
   particles.computeGridCellIndices();

   //test particle grid index
   particles.computeParticleCellIndices();

   particles.sortParticles();

   //prep the neighbor search
   neighborSearchHost nbs( particles, 4  );
   nbs.particlesToCells();
   nbs.runCycleOverGrid();


   //for(unsigned int i = 0; i < particles.getNumberOfParticles(); i++)
   //{
   //  std::cout << "Particle id: " << i << " has number of nbs: ";
   //  std::cout << particles.getNeighborsCount( i );
   //  std::cout << " with id: [ ";
   //  for(int j = 0; j < particles.getNeighborsCount( i ); j++)
   //    std::cout << particles.getNeighbor(i, j) << " ";
   //  std::cout << "]"<< std::endl;
   //}


   //test neighbor list
   // - test number of neigbors
   EXPECT_EQ( particles.getNeighborsCount( 0 ), 8 );
   EXPECT_EQ( particles.getNeighborsCount( 1 ), 5 );
   EXPECT_EQ( particles.getNeighborsCount( 2 ), 3 );
   EXPECT_EQ( particles.getNeighborsCount( 3 ), 5 );
   EXPECT_EQ( particles.getNeighborsCount( 4 ), 5 );
   EXPECT_EQ( particles.getNeighborsCount( 5 ), 3 );
   EXPECT_EQ( particles.getNeighborsCount( 6 ), 3 );
   EXPECT_EQ( particles.getNeighborsCount( 7 ), 5 );
   EXPECT_EQ( particles.getNeighborsCount( 8 ), 3 );

   // - test neighbors
   EXPECT_EQ( particles.getNeighbor( 0 , 0 ), 6 );
   EXPECT_EQ( particles.getNeighbor( 0 , 1 ), 7 );
   EXPECT_EQ( particles.getNeighbor( 0 , 2 ), 8 );
   EXPECT_EQ( particles.getNeighbor( 0 , 3 ), 1 );
   EXPECT_EQ( particles.getNeighbor( 0 , 4 ), 2 );
   EXPECT_EQ( particles.getNeighbor( 0 , 5 ), 3 );
   EXPECT_EQ( particles.getNeighbor( 0 , 6 ), 4 );
   EXPECT_EQ( particles.getNeighbor( 0 , 7 ), 5 );

   EXPECT_EQ( particles.getNeighbor( 1 , 0 ), 6 );
   EXPECT_EQ( particles.getNeighbor( 1 , 1 ), 7 );
   EXPECT_EQ( particles.getNeighbor( 1 , 2 ), 0 );
   EXPECT_EQ( particles.getNeighbor( 1 , 3 ), 2 );
   EXPECT_EQ( particles.getNeighbor( 1 , 4 ), 3 );

   EXPECT_EQ( particles.getNeighbor( 2 , 0 ), 0 );
   EXPECT_EQ( particles.getNeighbor( 2 , 1 ), 1 );
   EXPECT_EQ( particles.getNeighbor( 2 , 2 ), 3 );

   EXPECT_EQ( particles.getNeighbor( 3 , 0 ), 0 );
   EXPECT_EQ( particles.getNeighbor( 3 , 1 ), 1 );
   EXPECT_EQ( particles.getNeighbor( 3 , 2 ), 2 );
   EXPECT_EQ( particles.getNeighbor( 3 , 3 ), 4 );
   EXPECT_EQ( particles.getNeighbor( 3 , 4 ), 5 );

   EXPECT_EQ( particles.getNeighbor( 4 , 0 ), 7 );
   EXPECT_EQ( particles.getNeighbor( 4 , 1 ), 8 );
   EXPECT_EQ( particles.getNeighbor( 4 , 2 ), 0 );
   EXPECT_EQ( particles.getNeighbor( 4 , 3 ), 3 );
   EXPECT_EQ( particles.getNeighbor( 4 , 4 ), 5 );

   EXPECT_EQ( particles.getNeighbor( 5 , 0 ), 0 );
   EXPECT_EQ( particles.getNeighbor( 5 , 1 ), 3 );
   EXPECT_EQ( particles.getNeighbor( 5 , 2 ), 4 );

   EXPECT_EQ( particles.getNeighbor( 6 , 0 ), 7 );
   EXPECT_EQ( particles.getNeighbor( 6 , 1 ), 0 );
   EXPECT_EQ( particles.getNeighbor( 6 , 2 ), 1 );

   EXPECT_EQ( particles.getNeighbor( 7 , 0 ), 6 );
   EXPECT_EQ( particles.getNeighbor( 7 , 1 ), 8 );
   EXPECT_EQ( particles.getNeighbor( 7 , 2 ), 0 );
   EXPECT_EQ( particles.getNeighbor( 7 , 3 ), 1 );
   EXPECT_EQ( particles.getNeighbor( 7 , 4 ), 4 );

   EXPECT_EQ( particles.getNeighbor( 8 , 0 ), 7 );
   EXPECT_EQ( particles.getNeighbor( 8 , 1 ), 0 );
   EXPECT_EQ( particles.getNeighbor( 8 , 2 ), 4 );

}

TEST( SearchForNeighborsTest, NeighborSearchHostCentric )
{
   using ParticlesConfigCentric =  ParticlesConfigCentric< TNL::Devices::Host > ;
   using ParticlesHost = Particles< ParticlesConfigCentric, Devices::Host >;
   using neighborSearchHost = NeighborSearch< ParticlesConfigCentric, ParticlesHost >;

   using PointType = typename ParticlesHost::PointType;

   ParticlesHost particles;

   //test scuccessful creation of particle object
   ASSERT_TRUE( generate2DParticleSystemCentric( particles ) );

   //test number of points
   EXPECT_EQ( particles.getNumberOfParticles(), 9 );

   //test particle position
   PointType point0( -0.5, -0.5 ),  point1( -0.5, 0. ),  point2( -0.5, 0.5 ),
             point3( 0., -0.5 ), point4( 0., 0.0 ), point5( 0., 0.5 ),
             point6( 0.5, -0.5 ),  point7( 0.5, 0. ),  point8( 0.5, 0.5 );

   EXPECT_EQ( particles.getPoint( 0 ), point0 );
   EXPECT_EQ( particles.getPoint( 1 ), point1 );
   EXPECT_EQ( particles.getPoint( 2 ), point2 );
   EXPECT_EQ( particles.getPoint( 3 ), point3 );
   EXPECT_EQ( particles.getPoint( 4 ), point4 );
   EXPECT_EQ( particles.getPoint( 5 ), point5 );
   EXPECT_EQ( particles.getPoint( 6 ), point6 );
   EXPECT_EQ( particles.getPoint( 7 ), point7 );
   EXPECT_EQ( particles.getPoint( 8 ), point8 );

   //test grid index
   particles.computeGridCellIndices();
   //maybe even this: std::cout << "Grid indices: " << particles.getGridCellIndices() << std::endl;

   //test particle grid index
   particles.computeParticleCellIndices();

   particles.sortParticles();
   //test this: std::cout << "\nPrint sorted particle arrays. " << std::endl;
   //test this: std::cout << "Particle points: " << particles.getPoints() << std::endl;
   //test this: std::cout << "Grid indices: " << particles.getGridCellIndices() << std::endl;
   //test this: std::cout << "Particle indices: " << particles.getParticleCellIndices() << std::endl;

   //prep the neighbor search
   neighborSearchHost nbs( particles, 4  );
   nbs.particlesToCells();
   nbs.runCycleOverGrid();

   //test neighbor list
   // - test number of neigbors
   EXPECT_EQ( particles.getNeighborsCount( 0 ), 8 );
   EXPECT_EQ( particles.getNeighborsCount( 1 ), 5 );
   EXPECT_EQ( particles.getNeighborsCount( 2 ), 3 );
   EXPECT_EQ( particles.getNeighborsCount( 3 ), 5 );
   EXPECT_EQ( particles.getNeighborsCount( 4 ), 5 );
   EXPECT_EQ( particles.getNeighborsCount( 5 ), 3 );
   EXPECT_EQ( particles.getNeighborsCount( 6 ), 3 );
   EXPECT_EQ( particles.getNeighborsCount( 7 ), 5 );
   EXPECT_EQ( particles.getNeighborsCount( 8 ), 3 );

   // - test neighbors
   EXPECT_EQ( particles.getNeighbor( 0 , 0 ), 6 );
   EXPECT_EQ( particles.getNeighbor( 0 , 1 ), 7 );
   EXPECT_EQ( particles.getNeighbor( 0 , 2 ), 8 );
   EXPECT_EQ( particles.getNeighbor( 0 , 3 ), 1 );
   EXPECT_EQ( particles.getNeighbor( 0 , 4 ), 2 );
   EXPECT_EQ( particles.getNeighbor( 0 , 5 ), 3 );
   EXPECT_EQ( particles.getNeighbor( 0 , 6 ), 4 );
   EXPECT_EQ( particles.getNeighbor( 0 , 7 ), 5 );

   EXPECT_EQ( particles.getNeighbor( 1 , 0 ), 6 );
   EXPECT_EQ( particles.getNeighbor( 1 , 1 ), 7 );
   EXPECT_EQ( particles.getNeighbor( 1 , 2 ), 0 );
   EXPECT_EQ( particles.getNeighbor( 1 , 3 ), 2 );
   EXPECT_EQ( particles.getNeighbor( 1 , 4 ), 3 );

   EXPECT_EQ( particles.getNeighbor( 2 , 0 ), 0 );
   EXPECT_EQ( particles.getNeighbor( 2 , 1 ), 1 );
   EXPECT_EQ( particles.getNeighbor( 2 , 2 ), 3 );

   EXPECT_EQ( particles.getNeighbor( 3 , 0 ), 0 );
   EXPECT_EQ( particles.getNeighbor( 3 , 1 ), 1 );
   EXPECT_EQ( particles.getNeighbor( 3 , 2 ), 2 );
   EXPECT_EQ( particles.getNeighbor( 3 , 3 ), 4 );
   EXPECT_EQ( particles.getNeighbor( 3 , 4 ), 5 );

   EXPECT_EQ( particles.getNeighbor( 4 , 0 ), 7 );
   EXPECT_EQ( particles.getNeighbor( 4 , 1 ), 8 );
   EXPECT_EQ( particles.getNeighbor( 4 , 2 ), 0 );
   EXPECT_EQ( particles.getNeighbor( 4 , 3 ), 3 );
   EXPECT_EQ( particles.getNeighbor( 4 , 4 ), 5 );

   EXPECT_EQ( particles.getNeighbor( 5 , 0 ), 0 );
   EXPECT_EQ( particles.getNeighbor( 5 , 1 ), 3 );
   EXPECT_EQ( particles.getNeighbor( 5 , 2 ), 4 );

   EXPECT_EQ( particles.getNeighbor( 6 , 0 ), 7 );
   EXPECT_EQ( particles.getNeighbor( 6 , 1 ), 0 );
   EXPECT_EQ( particles.getNeighbor( 6 , 2 ), 1 );

   EXPECT_EQ( particles.getNeighbor( 7 , 0 ), 6 );
   EXPECT_EQ( particles.getNeighbor( 7 , 1 ), 8 );
   EXPECT_EQ( particles.getNeighbor( 7 , 2 ), 0 );
   EXPECT_EQ( particles.getNeighbor( 7 , 3 ), 1 );
   EXPECT_EQ( particles.getNeighbor( 7 , 4 ), 4 );

   EXPECT_EQ( particles.getNeighbor( 8 , 0 ), 7 );
   EXPECT_EQ( particles.getNeighbor( 8 , 1 ), 0 );
   EXPECT_EQ( particles.getNeighbor( 8 , 2 ), 4 );

}


/*
TEST( SearchForNeighborsTest, NeighborSearchCuda )
{

  using ParticlesCuda = Particles< ParticlesConfig, Deivces::Cuda >;

}
*/


