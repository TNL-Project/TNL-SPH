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

   static constexpr float searchRadius = 1;
   static constexpr int gridXsize = 2;
   static constexpr int gridYsize = 2;

   using NeighborListType = typename Algorithms::Segments::Ellpack< DeviceType, GlobalIndexType >;
};

template< typename Device >
bool generate2DParticleSystem( Particles< ParticlesConfig< Device >, Device > particles)
{

  //

  //setup points
  particles.setPoint(0, {0., 0.});
  particles.setPoint(1, {0., 0.5});
  particles.setPoint(2, {0., 1.});

  particles.setPoint(3, {0.5, 0.});
  particles.setPoint(4, {0.5, 0.5});
  particles.setPoint(5, {0.5, 1.});

  particles.setPoint(6, {1., 0.});
  particles.setPoint(7, {1., 0.5});
  particles.setPoint(8, {1., 1.});


  return true;
}

TEST(GetTwoTest, Two) {
    EXPECT_EQ(2, 2);
}

TEST( SearchForNeighborsTest, NeighborSearchHost )
{
   //using ParticlesHost = Particles< ParticlesConfig< TNL::Devices::Host>, Devices::Host >;

   //ParticlesHost particles;
   //ASSERT_TRUE( generate2DParticleSystem( particles ) );


}

/*
TEST( SearchForNeighborsTest, NeighborSearchCuda )
{

  using ParticlesCuda = Particles< ParticlesConfig, Deivces::Cuda >;

}
*/


