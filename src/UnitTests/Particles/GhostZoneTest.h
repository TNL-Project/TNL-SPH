#pragma once

#include <gtest/gtest.h>
#include <Particles/ParticlesLinkedListFloating.h>
#include <Particles/GhostZone.h>
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

TEST( GhostZonesConstruction2DTest, CollectParticlesInZoneCuda )
{
   using Device = TNL::Devices::Cuda;
   using ParticlesSetup = Particles2DSetup< Device >;
   using Particles = TNL::ParticleSystem::ParticlesLinkedList< ParticlesSetup::ParticlesConfig, Device >;
   using ParticlesPointer = typename Pointers::SharedPointer< Particles, Device >;

   using GhostZone = ParticleZone< ParticlesSetup::ParticlesConfig >;
   using IndexVectorType = typename GhostZone::IndexVectorType;

   ParticlesSetup setup;
   ParticlesPointer particles( setup.numberOfParticles,
                               setup.numberOfAllocatedParticles,
                               setup.searchRadius,
                               setup.numberOfGridCells );
   particles->setGridSize( setup.gridSize );
   particles->setGridOrigin( setup.gridOrigin );

   //assgn particles
   ASSERT_TRUE( assignPoints2D( particles ) );

   GhostZone zone_A( 3, 5 );
   IndexVectorType zone_A_begin = { 1, 1 };
   IndexVectorType zone_A_size = { 0, 2 };

   //:GhostZone zone_B( 3, 5 );
   //:GhostZone zone_C( 4, 5 );
   //:GhostZone zone_D( 3, 5 );

   zone_A.template assignCells< Particles::CellIndexer >( zone_A_begin, zone_A_size, particles->getGridSize() );
   std::cout << zone_A.getCellsInZone() << std::endl;
   std::cout << zone_A.getParticlesInZone() << std::endl;

}

