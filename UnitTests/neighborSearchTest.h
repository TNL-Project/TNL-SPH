#pragma once

// #ifdef HAVE_GTEST
#include <gtest/gtest.h>

#include <TNL/Algorithms/Segments/Ellpack.h>

#include "../Particles/Particles.h"
#include "../Particles/neighbourSearch.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h" //test vectors

using namespace TNL;
using namespace ParticleSystem;

/**
 * Particle config to test particle set of 9 particles in firt quadrant.
 * ( i.e. positions of all particles are positive)
 * NeighborSearch: SimpleCellIndex
 */
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
   //using CellIndexerType = ZOrderCurve<ParticlesConfig< DeviceType >, DeviceType>;

   static constexpr int spaceDimension = 2;
   static constexpr int maxOfNeigborsPerParticle = 9;

   static constexpr float searchRadius = 0.75;
   static constexpr int gridXsize = 2 + 2;
   static constexpr int gridYsize = 2 + 2;
   static constexpr RealType gridXbegin = 0 - searchRadius * 1;
   static constexpr RealType gridYbegin = 0 - searchRadius * 1;

   using NeighborListType = typename Algorithms::Segments::Ellpack< DeviceType, GlobalIndexType >;
};

/**
 * Particle config to test particle set of 16 particles placed in the center
 * of coordinate system.
 * ( i.e. positions of particles have all signs)
 * NeighborSearch: SimpleCellIndex
 */
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
   static constexpr int maxOfNeigborsPerParticle = 12;

   static constexpr float searchRadius = 0.6;
   static constexpr int gridXsize = 2 + 2;
   static constexpr int gridYsize = 3 + 2;
   static constexpr RealType gridXbegin = -0.5 - searchRadius * 1;
   static constexpr RealType gridYbegin = -0.5 - searchRadius * 1;

   using NeighborListType = typename Algorithms::Segments::Ellpack< DeviceType, GlobalIndexType >;
};

/**
 * Generate particle system of 9 particles in first quadrant.
 */
//template< typename Device >
template< typename Device, typename ParticlePointer >
//bool generate2DParticleSystem( Particles< ParticlesConfig< Device >, Device >& particles)
bool generate2DParticleSystem( ParticlePointer particles )
{
  //constructor:
  //particles = Particles< ParticlesConfig< Device >, Device >( 9, 0.75 );

  //setup points:
  particles->setPoint( 0, { 0., 0. } );
  particles->setPoint( 1, { 0., 0.5 } );
  particles->setPoint( 2, { 0., 1. } );

  particles->setPoint( 3, { 0.5, 0. } );
  particles->setPoint( 4, { 0.5, 0.5 } );
  particles->setPoint( 5, { 0.5, 1. } );

  particles->setPoint( 6, { 1., 0. } );
  particles->setPoint( 7, { 1., 0.5 } );
  particles->setPoint( 8, { 1., 1. } );

  return true;
}

template< typename Device, typename ParticlePointer >
//bool generate2DParticleSystemCentric( Particles< ParticlesConfigCentric< Device >, Device >& particles)
bool generate2DParticleSystemCentric( ParticlePointer particles )
{

  //constructor:
  //particles = Particles< ParticlesConfigCentric< Device >, Device >( 12, 0.6 );

  //setup points:
  particles->setPoint( 0, { -0.5, -0.5 } );
  particles->setPoint( 1, { -0.5, 0. } );
  particles->setPoint( 2, { -0.5, 0.5 } );
  particles->setPoint( 3, { -0.5, 1. } );

  particles->setPoint( 4, { 0., -0.5 } );
  particles->setPoint( 5, { 0., 0. } );
  particles->setPoint( 6, { 0., 0.5 } );
  particles->setPoint( 7, { 0., 1. } );

  particles->setPoint( 8, { 0.5, -0.5 } );
  particles->setPoint( 9, { 0.5, 0. } );
  particles->setPoint( 10, { 0.5, 0.5 } );
  particles->setPoint( 11, { 0.5, 1. } );

  return true;
}


TEST( SearchForNeighborsTest, NeighborSearchHost )
{
   using ParticlesConfig =  ParticlesConfig< TNL::Devices::Host > ;
   using ParticlesHost = Particles< ParticlesConfig, Devices::Host >;
   using neighborSearchHost = NeighborSearch< ParticlesConfig, ParticlesHost >;
	 using NeighborSearchPointer = typename Pointers::SharedPointer< neighborSearchHost, Devices::Host >;


	 using PointType = typename ParticlesHost::PointType;
	 using ParticlePointer = typename Pointers::SharedPointer< ParticlesHost, TNL::Devices::Host >;

   //ParticlesHost particles;
   ParticlePointer particles( 9, 0.75 );

   ASSERT_TRUE( ( generate2DParticleSystem< TNL::Devices::Host, ParticlePointer >( particles ) )  );

   //test number of points
   EXPECT_EQ( particles->getNumberOfParticles(), 9 );

   //test particle position
   PointType point0( 0., 0. ),  point1( 0., 0.5 ),  point2( 0., 1. ),
             point3( 0.5, 0. ), point4( 0.5, 0.5 ), point5( 0.5, 1. ),
             point6( 1., 0. ),  point7( 1., 0.5 ),  point8( 1., 1. );

   EXPECT_EQ( particles->getPoint( 0 ), point0 );
   EXPECT_EQ( particles->getPoint( 1 ), point1 );
   EXPECT_EQ( particles->getPoint( 2 ), point2 );
   EXPECT_EQ( particles->getPoint( 3 ), point3 );
   EXPECT_EQ( particles->getPoint( 4 ), point4 );
   EXPECT_EQ( particles->getPoint( 5 ), point5 );
   EXPECT_EQ( particles->getPoint( 6 ), point6 );
   EXPECT_EQ( particles->getPoint( 7 ), point7 );
   EXPECT_EQ( particles->getPoint( 8 ), point8 );

   //compute grid and particle cell index, sort particles:
   // ... this is something that would be also good to test
   particles->computeGridCellIndices();
   particles->computeParticleCellIndices();
   particles->sortParticles();

   //create neighbor search and search for neighbors:
   NeighborSearchPointer nbs( particles, ParticlesConfig::gridXsize * ParticlesConfig::gridYsize  );
   nbs->searchForNeighbors();
   //: //nbs->searchForNeighborsClean();

	//: //std::cout << "Nbsearch getCellFirstParticleList: " << nbs->getCellFirstParticleList() << std::endl;
	//: //std::cout << "Nbsearch getCellLastParticleList: " << nbs->getCellLastParticleList() << std::endl;

	//: //std::cout << "Particle positions: " << particles->getPoints() << std::endl;
	std::cout << "Particle cell indices: " << particles->getParticleCellIndices() << std::endl;
	std::cout << "Grid cell indices: " << particles->getGridCellIndices() << std::endl;
	//: //:	particles->GetParticlesInformations();

   //: //test neighbor list
   //: // - test number of neigbors
	//:  for ( unsigned int p = 0; p < 9; p++ )
	//: 	  if ( particles->getPoint( p ) == point4 )
   //: 		 	 EXPECT_EQ( particles->getNeighborsCount( p ), 8 );
	//:  		else if ( particles->getPoint( p ) == point1 )
   //: 			 EXPECT_EQ( particles->getNeighborsCount( p ), 5 );
	//:  		else if ( particles->getPoint( p ) == point3 )
   //: 			 EXPECT_EQ( particles->getNeighborsCount( p ), 5 );
	//:  		else if ( particles->getPoint( p ) == point5 )
   //: 			 EXPECT_EQ( particles->getNeighborsCount( p ), 5 );
	//:  		else if ( particles->getPoint( p ) == point7 )
   //: 			 EXPECT_EQ( particles->getNeighborsCount( p ), 5 );
	//:  		else if ( particles->getPoint( p ) == point0 )
   //: 			 EXPECT_EQ( particles->getNeighborsCount( p ), 3 );
	//:  		else if ( particles->getPoint( p ) == point2 )
   //: 			 EXPECT_EQ( particles->getNeighborsCount( p ), 3 );
	//:  		else if ( particles->getPoint( p ) == point6 )
   //: 			 EXPECT_EQ( particles->getNeighborsCount( p ), 3 );
	//:  		else if ( particles->getPoint( p ) == point8 )
   //: 			 EXPECT_EQ( particles->getNeighborsCount( p ), 3 );
	//:  		else
	//: 			std::cout << "Unexpected error." << std::endl;

  // // - test neighbors
  // EXPECT_EQ( particles->getNeighbor( 0 , 0 ), 6 );
  // EXPECT_EQ( particles->getNeighbor( 0 , 1 ), 7 );
  // EXPECT_EQ( particles->getNeighbor( 0 , 2 ), 8 );
  // EXPECT_EQ( particles->getNeighbor( 0 , 3 ), 1 );
  // EXPECT_EQ( particles->getNeighbor( 0 , 4 ), 2 );
  // EXPECT_EQ( particles->getNeighbor( 0 , 5 ), 3 );
  // EXPECT_EQ( particles->getNeighbor( 0 , 6 ), 4 );
  // EXPECT_EQ( particles->getNeighbor( 0 , 7 ), 5 );

	 std::vector< int > v = {1, 2, 3};
	 //ASSERT_THAT(v, ElementsAre(1, 2, 3));
	 //ASSERT({1, 2, 3}, ElementsAre(1, 2, 3));
		//ASSERT_THAT({2, 3, 4}, ElementsAre(5, 10, 15));
	 //EXPECT_THAT( v , ElementsMapOnto( 1 , 2 , 3 ) );
	 //using testing::ElementsAre;
	 //using testing::ElementsMapOnto;
	 ////EXPECT_THAT( v , ElementsAre( 2 , 1 , 3 ) );
	 //EXPECT_THAT( v , ElementsMapOnto( 1 , 2 , 3 ) );

	 /*
   EXPECT_EQ( particles->getNeighbor( 1 , 0 ), 6 );
   EXPECT_EQ( particles->getNeighbor( 1 , 1 ), 7 );
   EXPECT_EQ( particles->getNeighbor( 1 , 2 ), 0 );
   EXPECT_EQ( particles->getNeighbor( 1 , 3 ), 2 );
   EXPECT_EQ( particles->getNeighbor( 1 , 4 ), 3 );

   EXPECT_EQ( particles->getNeighbor( 2 , 0 ), 0 );
   EXPECT_EQ( particles->getNeighbor( 2 , 1 ), 1 );
   EXPECT_EQ( particles->getNeighbor( 2 , 2 ), 3 );

   EXPECT_EQ( particles->getNeighbor( 3 , 0 ), 0 );
   EXPECT_EQ( particles->getNeighbor( 3 , 1 ), 1 );
   EXPECT_EQ( particles->getNeighbor( 3 , 2 ), 2 );
   EXPECT_EQ( particles->getNeighbor( 3 , 3 ), 4 );
   EXPECT_EQ( particles->getNeighbor( 3 , 4 ), 5 );

   EXPECT_EQ( particles->getNeighbor( 4 , 0 ), 7 );
   EXPECT_EQ( particles->getNeighbor( 4 , 1 ), 8 );
   EXPECT_EQ( particles->getNeighbor( 4 , 2 ), 0 );
   EXPECT_EQ( particles->getNeighbor( 4 , 3 ), 3 );
   EXPECT_EQ( particles->getNeighbor( 4 , 4 ), 5 );

   EXPECT_EQ( particles->getNeighbor( 5 , 0 ), 0 );
   EXPECT_EQ( particles->getNeighbor( 5 , 1 ), 3 );
   EXPECT_EQ( particles->getNeighbor( 5 , 2 ), 4 );

   EXPECT_EQ( particles->getNeighbor( 6 , 0 ), 7 );
   EXPECT_EQ( particles->getNeighbor( 6 , 1 ), 0 );
   EXPECT_EQ( particles->getNeighbor( 6 , 2 ), 1 );

   EXPECT_EQ( particles->getNeighbor( 7 , 0 ), 6 );
   EXPECT_EQ( particles->getNeighbor( 7 , 1 ), 8 );
   EXPECT_EQ( particles->getNeighbor( 7 , 2 ), 0 );
   EXPECT_EQ( particles->getNeighbor( 7 , 3 ), 1 );
   EXPECT_EQ( particles->getNeighbor( 7 , 4 ), 4 );

   EXPECT_EQ( particles->getNeighbor( 8 , 0 ), 7 );
   EXPECT_EQ( particles->getNeighbor( 8 , 1 ), 0 );
   EXPECT_EQ( particles->getNeighbor( 8 , 2 ), 4 );
	 */
}

//: #ifdef HAVE_CUDA
//: template< typename ParticlePointer >
//: __global__ void testSetPointKernel( ParticlePointer particles )
//: {
//:   particles->setPoint( 0, { 0., 0. } );
//:   particles->setPoint( 1, { 0., 0.5 } );
//:   particles->setPoint( 2, { 0., 1. } );
//:
//:   particles->setPoint( 3, { 0.5, 0. } );
//:   particles->setPoint( 4, { 0.5, 0.5 } );
//:   particles->setPoint( 5, { 0.5, 1. } );
//:
//:   particles->setPoint( 6, { 1., 0. } );
//:   particles->setPoint( 7, { 1., 0.5 } );
//:   particles->setPoint( 8, { 1., 1. } );
//: }
//:
//:
//: template< typename ParticlePointer >
//: __global__ void temp_printNeighbors( ParticlePointer particles )
//: {
//: 	if( threadIdx.x < 9 )
//: 		printf(" %d: %d\n", threadIdx.x , particles->getNeighborsCount( threadIdx.x ) );
//: }
//: #endif /* HAVE_CUDA */
//:
//: #ifdef HAVE_CUDA
//: TEST( SearchForNeighborsTest, NeighborSearchCuda )
//: {
//:    using ParticlesConfig =  ParticlesConfig< TNL::Devices::Cuda > ;
//:    using ParticlesCuda = Particles< ParticlesConfig, Devices::Cuda >;
//:    using neighborSearchCuda = NeighborSearch< ParticlesConfig, ParticlesCuda >;
//: 	 using NeighborSearchPointer = typename Pointers::SharedPointer< neighborSearchCuda, Devices::Cuda >;
//:
//:    using PointType = typename ParticlesCuda::PointType;
//: 	 using ParticlePointer = typename Pointers::SharedPointer< ParticlesCuda, TNL::Devices::Cuda >;
//:
//:    //ParticlesCuda particles;
//:    ParticlePointer particles( 9, 0.75 );
//:
//:    ////ASSERT_TRUE( generate2DParticleSystem( particles ) );
//: 	 testSetPointKernel< ParticlePointer ><<< 1, 1 >>>( particles );
//: 	 cudaDeviceSynchronize();
//: 	 TNL_CHECK_CUDA_DEVICE;
//:
//:    //test number of points
//:    EXPECT_EQ( particles->getNumberOfParticles(), 9 );
//:
//:    //test particle position
//:    PointType point0( 0., 0. ),  point1( 0., 0.5 ),  point2( 0., 1. ),
//:              point3( 0.5, 0. ), point4( 0.5, 0.5 ), point5( 0.5, 1. ),
//:              point6( 1., 0. ),  point7( 1., 0.5 ),  point8( 1., 1. );
//:
//:
//:
//:    EXPECT_EQ( particles->getPoints().getElement( 0 ), point0 ); //getPoint works in kernel
//:    EXPECT_EQ( particles->getPoints().getElement( 1 ), point1 );
//:    EXPECT_EQ( particles->getPoints().getElement( 2 ), point2 );
//:    EXPECT_EQ( particles->getPoints().getElement( 3 ), point3 );
//:    EXPECT_EQ( particles->getPoints().getElement( 4 ), point4 );
//:    EXPECT_EQ( particles->getPoints().getElement( 5 ), point5 );
//:    EXPECT_EQ( particles->getPoints().getElement( 6 ), point6 );
//:    EXPECT_EQ( particles->getPoints().getElement( 7 ), point7 );
//:    EXPECT_EQ( particles->getPoints().getElement( 8 ), point8 );
//:
//:    //compute grid and particle cell index, sort particles:
//:    // ... this is something that would be also good to test
//:    particles->computeGridCellIndices();
//:    particles->computeParticleCellIndices();
//:    particles->sortParticles();
//:
//:    //create neighbor search and search for neighbors:
//:    NeighborSearchPointer nbs( particles, ParticlesConfig::gridXsize * ParticlesConfig::gridYsize  );
//:    nbs->searchForNeighbors();
//: 	 cudaDeviceSynchronize();
//: 	 TNL_CHECK_CUDA_DEVICE;
//:
//: 	 //std::cout << " Poitns: " << particles->getPoints( ) << std::endl;
//: 	 //std::cout << " CellIdx: " << particles->getParticleCellIndices( ) << std::endl;
//: 	 //std::cout << " GridIdx: " << particles->getGridCellIndices( ) << std::endl;
//: 	 //std::cout << " nbs first: " << nbs->getCellFirstParticleList( ) << std::endl;
//: 	 //std::cout << " nbs flast: " << nbs->getCellLastParticleList( ) << std::endl;
//: 	 //temp_printNeighbors< ParticlePointer ><<< 1, 9 >>>( particles );
//:
//: 	 for ( unsigned int p = 0; p < 9; p++ )
//: 		  if ( particles->getPoints().getElement( p ) == point4 )
//:    		 	 EXPECT_EQ( particles->getNeighborsCountList().getElement( p ), 8 );
//: 	 		else if ( particles->getPoints().getElement( p ) == point1 )
//:    			 EXPECT_EQ( particles->getNeighborsCountList().getElement( p ), 5 );
//: 	 		else if ( particles->getPoints().getElement( p ) == point3 )
//:    			 EXPECT_EQ( particles->getNeighborsCountList().getElement( p ), 5 );
//: 	 		else if ( particles->getPoints().getElement( p ) == point5 )
//:    			 EXPECT_EQ( particles->getNeighborsCountList().getElement( p ), 5 );
//: 	 		else if ( particles->getPoints().getElement( p ) == point7 )
//:    			 EXPECT_EQ( particles->getNeighborsCountList().getElement( p ), 5 );
//: 	 		else if ( particles->getPoints().getElement( p ) == point0 )
//:    			 EXPECT_EQ( particles->getNeighborsCountList().getElement( p ), 3 );
//: 	 		else if ( particles->getPoints().getElement( p ) == point2 )
//:    			 EXPECT_EQ( particles->getNeighborsCountList().getElement( p ), 3 );
//: 	 		else if ( particles->getPoints().getElement( p ) == point6 )
//:    			 EXPECT_EQ( particles->getNeighborsCountList().getElement( p ), 3 );
//: 	 		else if ( particles->getPoints().getElement( p ) == point8 )
//:    			 EXPECT_EQ( particles->getNeighborsCountList().getElement( p ), 3 );
//: 	 		else
//: 				std::cout << "Unexpected error." << std::endl;
//:
//:    //: //test neighbor list
//:    //: // - test number of neigbors
//:    //: EXPECT_EQ( particles->getNeighborsCount( 0 ), 8 );
//:    //: EXPECT_EQ( particles->getNeighborsCount( 1 ), 5 );
//:    //: EXPECT_EQ( particles.getNeighborsCount( 2 ), 3 );
//:    //: EXPECT_EQ( particles.getNeighborsCount( 3 ), 5 );
//:    //: EXPECT_EQ( particles.getNeighborsCount( 4 ), 5 );
//:    //: EXPECT_EQ( particles.getNeighborsCount( 5 ), 3 );
//:    //: EXPECT_EQ( particles.getNeighborsCount( 6 ), 3 );
//:    //: EXPECT_EQ( particles.getNeighborsCount( 7 ), 5 );
//:    //: EXPECT_EQ( particles.getNeighborsCount( 8 ), 3 );
//:
//:    //: // - test neighbors
//:    //: EXPECT_EQ( particles.getNeighbor( 0 , 0 ), 6 );
//:    //: EXPECT_EQ( particles.getNeighbor( 0 , 1 ), 7 );
//:    //: EXPECT_EQ( particles.getNeighbor( 0 , 2 ), 8 );
//:    //: EXPECT_EQ( particles.getNeighbor( 0 , 3 ), 1 );
//:    //: EXPECT_EQ( particles.getNeighbor( 0 , 4 ), 2 );
//:    //: EXPECT_EQ( particles.getNeighbor( 0 , 5 ), 3 );
//:    //: EXPECT_EQ( particles.getNeighbor( 0 , 6 ), 4 );
//:    //: EXPECT_EQ( particles.getNeighbor( 0 , 7 ), 5 );
//:
//:    //: EXPECT_EQ( particles.getNeighbor( 1 , 0 ), 6 );
//:    //: EXPECT_EQ( particles.getNeighbor( 1 , 1 ), 7 );
//:    //: EXPECT_EQ( particles.getNeighbor( 1 , 2 ), 0 );
//:    //: EXPECT_EQ( particles.getNeighbor( 1 , 3 ), 2 );
//:    //: EXPECT_EQ( particles.getNeighbor( 1 , 4 ), 3 );
//:
//:    //: EXPECT_EQ( particles.getNeighbor( 2 , 0 ), 0 );
//:    //: EXPECT_EQ( particles.getNeighbor( 2 , 1 ), 1 );
//:    //: EXPECT_EQ( particles.getNeighbor( 2 , 2 ), 3 );
//:
//:    //: EXPECT_EQ( particles.getNeighbor( 3 , 0 ), 0 );
//:    //: EXPECT_EQ( particles.getNeighbor( 3 , 1 ), 1 );
//:    //: EXPECT_EQ( particles.getNeighbor( 3 , 2 ), 2 );
//:    //: EXPECT_EQ( particles.getNeighbor( 3 , 3 ), 4 );
//:    //: EXPECT_EQ( particles.getNeighbor( 3 , 4 ), 5 );
//:
//:    //: EXPECT_EQ( particles.getNeighbor( 4 , 0 ), 7 );
//:    //: EXPECT_EQ( particles.getNeighbor( 4 , 1 ), 8 );
//:    //: EXPECT_EQ( particles.getNeighbor( 4 , 2 ), 0 );
//:    //: EXPECT_EQ( particles.getNeighbor( 4 , 3 ), 3 );
//:    //: EXPECT_EQ( particles.getNeighbor( 4 , 4 ), 5 );
//:
//:    //: EXPECT_EQ( particles.getNeighbor( 5 , 0 ), 0 );
//:    //: EXPECT_EQ( particles.getNeighbor( 5 , 1 ), 3 );
//:    //: EXPECT_EQ( particles.getNeighbor( 5 , 2 ), 4 );
//:
//:    //: EXPECT_EQ( particles.getNeighbor( 6 , 0 ), 7 );
//:    //: EXPECT_EQ( particles.getNeighbor( 6 , 1 ), 0 );
//:    //: EXPECT_EQ( particles.getNeighbor( 6 , 2 ), 1 );
//:
//:    //: EXPECT_EQ( particles.getNeighbor( 7 , 0 ), 6 );
//:    //: EXPECT_EQ( particles.getNeighbor( 7 , 1 ), 8 );
//:    //: EXPECT_EQ( particles.getNeighbor( 7 , 2 ), 0 );
//:    //: EXPECT_EQ( particles.getNeighbor( 7 , 3 ), 1 );
//:    //: EXPECT_EQ( particles.getNeighbor( 7 , 4 ), 4 );
//:
//:    //: EXPECT_EQ( particles.getNeighbor( 8 , 0 ), 7 );
//:    //: EXPECT_EQ( particles.getNeighbor( 8 , 1 ), 0 );
//:    //: EXPECT_EQ( particles.getNeighbor( 8 , 2 ), 4 );
//: }
//: #endif /* HAVE_CUDA */
//:
//: TEST( SearchForNeighborsTest, NeighborSearchHostCentric )
//: {
//:    using ParticlesConfigCentric =  ParticlesConfigCentric< TNL::Devices::Host > ;
//:    using ParticlesHost = Particles< ParticlesConfigCentric, Devices::Host >;
//:    using neighborSearchHost = NeighborSearch< ParticlesConfigCentric, ParticlesHost >;
//: 	 using NeighborSearchPointer = typename Pointers::SharedPointer< neighborSearchHost, Devices::Cuda >;
//:
//:    using PointType = typename ParticlesHost::PointType;
//: 	 using ParticlePointer = typename Pointers::SharedPointer< ParticlesHost, TNL::Devices::Host >;
//:
//:    ParticlePointer particles( 12, 0.6 );
//:
//:    //test scuccessful creation of particle object
//:    ASSERT_TRUE( ( generate2DParticleSystemCentric< TNL::Devices::Host, ParticlePointer >( particles ) ) );
//:
//:    //test number of points
//:    EXPECT_EQ( particles->getNumberOfParticles(), 12 );
//:
//:    //test particle position
//:    PointType point0( -0.5, -0.5 ), point1( -0.5, 0. ), point2( -0.5, 0.5 ), point3( -0.5, 1. ),
//:              point4( 0., -0.5 ),   point5( 0., 0.0 ),  point6( 0., 0.5 ),   point7( 0, 1. ),
//:              point8( 0.5, -0.5 ),  point9( 0.5, 0. ),  point10( 0.5, 0.5 ), point11( 0.5, 1. );
//:
//:    EXPECT_EQ( particles->getPoint( 0 ), point0 );
//:    EXPECT_EQ( particles->getPoint( 1 ), point1 );
//:    EXPECT_EQ( particles->getPoint( 2 ), point2 );
//:    EXPECT_EQ( particles->getPoint( 3 ), point3 );
//:    EXPECT_EQ( particles->getPoint( 4 ), point4 );
//:    EXPECT_EQ( particles->getPoint( 5 ), point5 );
//:    EXPECT_EQ( particles->getPoint( 6 ), point6 );
//:    EXPECT_EQ( particles->getPoint( 7 ), point7 );
//:    EXPECT_EQ( particles->getPoint( 8 ), point8 );
//:    EXPECT_EQ( particles->getPoint( 9 ), point9 );
//:    EXPECT_EQ( particles->getPoint( 10 ), point10 );
//:    EXPECT_EQ( particles->getPoint( 11 ), point11 );
//:
//:    //compute grid and particle cell index, sort particles:
//:    // ... this is something that would be also good to test
//:    particles->computeGridCellIndices();
//:    particles->computeParticleCellIndices();
//:    particles->sortParticles();
//:
//:    //create neighbor search and search for neighbors:
//:    NeighborSearchPointer nbs( particles, ParticlesConfigCentric::gridXsize * ParticlesConfigCentric::gridYsize  );
//:    nbs->searchForNeighbors();
//:
//: 	 for ( unsigned int p = 0; p < 12; p++ )
//: 		  if ( particles->getPoint( p ) == point0 )
//:    		 	 EXPECT_EQ( particles->getNeighborsCount( p ), 2 );
//: 	 		else if ( particles->getPoint( p ) == point1 )
//:    			 EXPECT_EQ( particles->getNeighborsCount( p ), 3 );
//: 	 		else if ( particles->getPoint( p ) == point2 )
//:    			 EXPECT_EQ( particles->getNeighborsCount( p ), 3 );
//: 	 		else if ( particles->getPoint( p ) == point3 )
//:    			 EXPECT_EQ( particles->getNeighborsCount( p ), 2 );
//: 	 		else if ( particles->getPoint( p ) == point4 )
//:    			 EXPECT_EQ( particles->getNeighborsCount( p ), 3 );
//: 	 		else if ( particles->getPoint( p ) == point5 )
//:    			 EXPECT_EQ( particles->getNeighborsCount( p ), 4 );
//: 	 		else if ( particles->getPoint( p ) == point6 )
//:    			 EXPECT_EQ( particles->getNeighborsCount( p ), 4 );
//: 	 		else if ( particles->getPoint( p ) == point7 )
//:    			 EXPECT_EQ( particles->getNeighborsCount( p ), 3 );
//: 	 		else if ( particles->getPoint( p ) == point8 )
//:    			 EXPECT_EQ( particles->getNeighborsCount( p ), 2 );
//: 	 		else if ( particles->getPoint( p ) == point9 )
//:    			 EXPECT_EQ( particles->getNeighborsCount( p ), 3 );
//: 	 		else if ( particles->getPoint( p ) == point10 )
//:    			 EXPECT_EQ( particles->getNeighborsCount( p ), 3 );
//: 	 		else if ( particles->getPoint( p ) == point11 )
//:    			 EXPECT_EQ( particles->getNeighborsCount( p ), 2 );
//: 	 		else
//: 				std::cout << "Unexpected error." << std::endl;
//:
//: 	 /*
//:    //test neighbor list
//:    // - test number of neigbors
//:    EXPECT_EQ( particles->getNeighborsCount( 0 ), 4 );
//:    EXPECT_EQ( particles->getNeighborsCount( 1 ), 3 );
//:    EXPECT_EQ( particles->getNeighborsCount( 2 ), 2 );
//:    EXPECT_EQ( particles->getNeighborsCount( 3 ), 3 );
//:    EXPECT_EQ( particles->getNeighborsCount( 4 ), 2 );
//:    EXPECT_EQ( particles->getNeighborsCount( 5 ), 3 );
//:    EXPECT_EQ( particles->getNeighborsCount( 6 ), 3 );
//:    EXPECT_EQ( particles->getNeighborsCount( 7 ), 4 );
//:    EXPECT_EQ( particles->getNeighborsCount( 8 ), 3 );
//:    EXPECT_EQ( particles->getNeighborsCount( 9 ), 2 );
//:    EXPECT_EQ( particles->getNeighborsCount( 10 ), 3 );
//:    EXPECT_EQ( particles->getNeighborsCount( 11 ), 2 );
//:
//:    // - test neighbors
//:    EXPECT_EQ( particles->getNeighbor( 0 , 0 ), 7 );
//:    EXPECT_EQ( particles->getNeighbor( 0 , 1 ), 1 );
//:    EXPECT_EQ( particles->getNeighbor( 0 , 2 ), 3 );
//:    EXPECT_EQ( particles->getNeighbor( 0 , 3 ), 5 );
//:
//:    EXPECT_EQ( particles->getNeighbor( 1 , 0 ), 0 );
//:    EXPECT_EQ( particles->getNeighbor( 1 , 1 ), 2 );
//:    EXPECT_EQ( particles->getNeighbor( 1 , 2 ), 4 );
//:
//:    EXPECT_EQ( particles->getNeighbor( 2 , 0 ), 1 );
//:    EXPECT_EQ( particles->getNeighbor( 2 , 1 ), 3 );
//:
//:    EXPECT_EQ( particles->getNeighbor( 3 , 0 ), 6 );
//:    EXPECT_EQ( particles->getNeighbor( 3 , 1 ), 0 );
//:    EXPECT_EQ( particles->getNeighbor( 3 , 2 ), 2 );
//:
//:    EXPECT_EQ( particles->getNeighbor( 4 , 0 ), 1 );
//:    EXPECT_EQ( particles->getNeighbor( 4 , 1 ), 5 );
//:
//:    EXPECT_EQ( particles->getNeighbor( 5 , 0 ), 8 );
//:    EXPECT_EQ( particles->getNeighbor( 5 , 1 ), 0 );
//:    EXPECT_EQ( particles->getNeighbor( 5 , 2 ), 4 );
//:
//:    EXPECT_EQ( particles->getNeighbor( 6 , 0 ), 9 );
//:    EXPECT_EQ( particles->getNeighbor( 6 , 1 ), 7 );
//:    EXPECT_EQ( particles->getNeighbor( 6 , 2 ), 3 );
//:
//:    EXPECT_EQ( particles->getNeighbor( 7 , 0 ), 10 );
//:    EXPECT_EQ( particles->getNeighbor( 7 , 1 ), 6 );
//:    EXPECT_EQ( particles->getNeighbor( 7 , 2 ), 8 );
//:    EXPECT_EQ( particles->getNeighbor( 7 , 3 ), 0 );
//:
//:    EXPECT_EQ( particles->getNeighbor( 8 , 0 ), 11 );
//:    EXPECT_EQ( particles->getNeighbor( 8 , 1 ), 7 );
//:    EXPECT_EQ( particles->getNeighbor( 8 , 2 ), 5 );
//:
//:    EXPECT_EQ( particles->getNeighbor( 9 , 0 ), 10 );
//:    EXPECT_EQ( particles->getNeighbor( 9 , 1 ), 6 );
//:
//:    EXPECT_EQ( particles->getNeighbor( 10 , 0 ), 9 );
//:    EXPECT_EQ( particles->getNeighbor( 10 , 1 ), 11 );
//:    EXPECT_EQ( particles->getNeighbor( 10 , 2 ), 7 );
//:
//:    EXPECT_EQ( particles->getNeighbor( 11 , 0 ), 10 );
//:    EXPECT_EQ( particles->getNeighbor( 11 , 1 ), 8 );
//: 	 */
//: }
//:
//: /*
//: TEST( SearchForNeighborsTest, NeighborSearchCuda )
//: {
//:
//:   using ParticlesCuda = Particles< ParticlesConfig, Deivces::Cuda >;
//:
//: }
//: */

