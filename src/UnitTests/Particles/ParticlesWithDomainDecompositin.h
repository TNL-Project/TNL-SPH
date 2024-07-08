#pragma once

#include <gtest/gtest.h>
#include "gtest/gtest.h"
#include "gmock/gmock.h" //test vectors

#include <Particles/ParticlesLinkedList.h>

using namespace TNL;
using namespace ParticleSystem;

// Config for 2D particle system with decomposition
template< typename Device >
class ParticlesWithDecomposition2DConfig
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
   using CellIndexerType = SimpleCellIndex< spaceDimension,
                                            ParticlesWithDecomposition2DConfig,
                                            std::index_sequence< 0, 1 > >;
};

// Setup (grouping structure) with input parameters for tested particle system
template< typename Device >
struct ParticlesWithDecomposition2DSetup
{
   using Config = ParticlesWithDecomposition2DConfig< Device >;
   using ParticlesTraitsType = ParticlesTraits< Config, Device >;
   using IndexVectorType = typename ParticlesTraitsType::IndexVectorType;
   using PointType = typename ParticlesTraitsType::PointType;

   // input parameters
   const int numberOfParticles = 20;
   const int numberOfAllocatedParticles = 28;

   const float searchRadius = 0.5;
   const PointType gridOrigin = { -searchRadius, -searchRadius };
   const IndexVectorType gridDimensions = { 6, 5 };

   // input parameters for enhanced decomposition paritlces
   const int overlapWidth = 1;
   const PointType gridReferentialOrigin = { -3 * searchRadius, -3 * searchRadius };
   const IndexVectorType gridOriginGlobalCoords = { 1, 1 };
};

TEST( ParticlesWithDecomposition2DTest, ParticlesPropertiesHost )
{
   using Device = TNL::Devices::Host;
   using Setup = ParticlesWithDecomposition2DSetup< Device >;
   using Particles = TNL::ParticleSystem::ParticlesLinkedList< Setup::Config, Device >;
   using ParticlesPointer = typename Pointers::SharedPointer< Particles, Device >;
   using PointType = typename Setup::PointType;
   using IndexVectorType = typename Setup::IndexVectorType;

   Setup setup;
   ParticlesPointer particles;
   // input parameters
   particles->setSize( setup.numberOfAllocatedParticles );
   particles->setNumberOfParticles( setup.numberOfParticles );

   particles->setSearchRadius( setup.searchRadius );
   particles->setGridOrigin( setup.gridOrigin );
   particles->setGridDimensions( setup.gridDimensions );

   // input parameters for enhanced decomposition paritlces
   particles->setOverlapWidth( setup.overlapWidth );
   particles->setGridReferentialOrigin( setup.gridReferentialOrigin );
   particles->setGridOriginGlobalCoords( setup.gridOriginGlobalCoords );

   // assign particles
   ASSERT_TRUE( assignPoints2D( particles ) );

   // test particle system default properties
   EXPECT_EQ( particles->getNumberOfParticles(), 20 );
   EXPECT_EQ( particles->getNumberOfAllocatedParticles(), 28 );
   EXPECT_EQ( particles->getPoints().getSize(), 28 );
   EXPECT_EQ( particles->getSortPermutations()->getSize(), 28 );
   EXPECT_EQ( particles->getParticleCellIndices().getSize(), 28 ); //cellList related

   EXPECT_EQ( particles->getSearchRadius(), 0.5 );
   const PointType gridOrigin = { -0.5, -0.5 };
   EXPECT_EQ( particles->getGridOrigin(), gridOrigin );
   const IndexVectorType gridDimensions = { 6, 5 };
   EXPECT_EQ( particles->getGridDimensions(), gridDimensions );
   EXPECT_EQ( particles->getCellFirstLastParticleList().getSize(), 56 );

   EXPECT_EQ( particles->getOverlapWidth(), 1 );
   const PointType gridReferentialOrigin = { -1.5, -1.5 };
   EXPECT_EQ( particles->getGridReferentialOrigin(), gridReferentialOrigin );
   const IndexVectorType gridOriginGlobalCoords = { 1, 1 };
   EXPECT_EQ( particles->getGridOriginGlobalCoords(), gridOriginGlobalCoords );

   // test particle system computed properties
   const PointType gridOriginWithOverlap = { -1.0, -1.0 };
   EXPECT_EQ( particles->getGridOriginWithOverlap(), gridOriginWithOverlap );
   const IndexVectorType gridDimensionsWithOverlap = { 8, 7 };
   EXPECT_EQ( particles->getGridDimensionsWithOverlap(), gridDimensionsWithOverlap );

   // change overlap, recompute sizes
   particles->setOverlapWidth( 3 );
   const PointType gridOriginWithOverlap_scaledUp = { -2.0, -2.0 };
   EXPECT_EQ( particles->getGridOriginWithOverlap(), gridOriginWithOverlap_scaledUp );
   const IndexVectorType gridDimensionsWithOverlap_scaledUp = { 12, 11 };
   EXPECT_EQ( particles->getGridDimensionsWithOverlap(), gridDimensionsWithOverlap_scaledUp );
   EXPECT_EQ( particles->getCellFirstLastParticleList().getSize(), 132 );

   particles->setOverlapWidth( 2 );
   const PointType gridOriginWithOverlap_scaledDown = { -1.5, -1.5 };
   EXPECT_EQ( particles->getGridOriginWithOverlap(), gridOriginWithOverlap_scaledDown );
   const IndexVectorType gridDimensionsWithOverlap_scaledDown = { 10, 9 };
   EXPECT_EQ( particles->getGridDimensionsWithOverlap(), gridDimensionsWithOverlap_scaledDown );
   EXPECT_EQ( particles->getCellFirstLastParticleList().getSize(), 90 );

   // change grid size
   particles->setGridDimensions( { 3, 4 } );
   const IndexVectorType gridDimensions_updated = { 3, 4 };
   EXPECT_EQ( particles->getGridDimensions(), gridDimensions_updated );
   const IndexVectorType gridDimensionsWithOverlap_updated = { 7, 8 };
   EXPECT_EQ( particles->getGridDimensionsWithOverlap(), gridDimensionsWithOverlap_updated );
   EXPECT_EQ( particles->getCellFirstLastParticleList().getSize(), 56 );
}

TEST( ParticlesWithDecomposition2DTest, ParticlesPropertiesDevice )
{
   using Device = TNL::Devices::Cuda;
   using Setup = ParticlesWithDecomposition2DSetup< Device >;
   using Particles = TNL::ParticleSystem::ParticlesLinkedList< Setup::Config, Device >;
   using ParticlesPointer = typename Pointers::SharedPointer< Particles, Device >;
   using PointType = typename Setup::PointType;
   using IndexVectorType = typename Setup::IndexVectorType;

   Setup setup;
   ParticlesPointer particles;
   // input parameters
   particles->setSize( setup.numberOfAllocatedParticles );
   particles->setNumberOfParticles( setup.numberOfParticles );

   particles->setSearchRadius( setup.searchRadius );
   particles->setGridOrigin( setup.gridOrigin );
   particles->setGridDimensions( setup.gridDimensions );

   // input parameters for enhanced decomposition paritlces
   particles->setOverlapWidth( setup.overlapWidth );
   particles->setGridReferentialOrigin( setup.gridReferentialOrigin );
   particles->setGridOriginGlobalCoords( setup.gridOriginGlobalCoords );

   // assign particles
   ASSERT_TRUE( assignPoints2D( particles ) );

   // test particle system default properties
   EXPECT_EQ( particles->getNumberOfParticles(), 20 );
   EXPECT_EQ( particles->getNumberOfAllocatedParticles(), 28 );
   EXPECT_EQ( particles->getPoints().getSize(), 28 );
   EXPECT_EQ( particles->getSortPermutations()->getSize(), 28 );
   EXPECT_EQ( particles->getParticleCellIndices().getSize(), 28 ); //cellList related

   EXPECT_EQ( particles->getSearchRadius(), 0.5 );
   const PointType gridOrigin = { -0.5, -0.5 };
   EXPECT_EQ( particles->getGridOrigin(), gridOrigin );
   const IndexVectorType gridDimensions = { 6, 5 };
   EXPECT_EQ( particles->getGridDimensions(), gridDimensions );
   EXPECT_EQ( particles->getCellFirstLastParticleList().getSize(), 56 );

   EXPECT_EQ( particles->getOverlapWidth(), 1 );
   const PointType gridReferentialOrigin = { -1.5, -1.5 };
   EXPECT_EQ( particles->getGridReferentialOrigin(), gridReferentialOrigin );
   const IndexVectorType gridOriginGlobalCoords = { 1, 1 };
   EXPECT_EQ( particles->getGridOriginGlobalCoords(), gridOriginGlobalCoords );

   // test particle system computed properties
   const PointType gridOriginWithOverlap = { -1.0, -1.0 };
   EXPECT_EQ( particles->getGridOriginWithOverlap(), gridOriginWithOverlap );
   const IndexVectorType gridDimensionsWithOverlap = { 8, 7 };
   EXPECT_EQ( particles->getGridDimensionsWithOverlap(), gridDimensionsWithOverlap );

   // change overlap, recompute sizes
   particles->setOverlapWidth( 3 );
   const PointType gridOriginWithOverlap_scaledUp = { -2.0, -2.0 };
   EXPECT_EQ( particles->getGridOriginWithOverlap(), gridOriginWithOverlap_scaledUp );
   const IndexVectorType gridDimensionsWithOverlap_scaledUp = { 12, 11 };
   EXPECT_EQ( particles->getGridDimensionsWithOverlap(), gridDimensionsWithOverlap_scaledUp );
   EXPECT_EQ( particles->getCellFirstLastParticleList().getSize(), 132 );

   particles->setOverlapWidth( 2 );
   const PointType gridOriginWithOverlap_scaledDown = { -1.5, -1.5 };
   EXPECT_EQ( particles->getGridOriginWithOverlap(), gridOriginWithOverlap_scaledDown );
   const IndexVectorType gridDimensionsWithOverlap_scaledDown = { 10, 9 };
   EXPECT_EQ( particles->getGridDimensionsWithOverlap(), gridDimensionsWithOverlap_scaledDown );
   EXPECT_EQ( particles->getCellFirstLastParticleList().getSize(), 90 );

   // change grid size
   particles->setGridDimensions( { 3, 4 } );
   const IndexVectorType gridDimensions_updated = { 3, 4 };
   EXPECT_EQ( particles->getGridDimensions(), gridDimensions_updated );
   const IndexVectorType gridDimensionsWithOverlap_updated = { 7, 8 };
   EXPECT_EQ( particles->getGridDimensionsWithOverlap(), gridDimensionsWithOverlap_updated );
   EXPECT_EQ( particles->getCellFirstLastParticleList().getSize(), 56 );
}

TEST( ParticlesWithDecomposition2DTest, ComputeParticleCellIndicesCuda )
{
   using Device = TNL::Devices::Cuda;
   using Setup = ParticlesWithDecomposition2DSetup< Device >;
   using Particles = TNL::ParticleSystem::ParticlesLinkedList< Setup::Config, Device >;
   using ParticlesPointer = typename Pointers::SharedPointer< Particles, Device >;
   using PointType = typename Setup::PointType;
   using IndexVectorType = typename Setup::IndexVectorType;

   Setup setup;
   ParticlesPointer particles;
   // input parameters
   particles->setSize( setup.numberOfAllocatedParticles );
   particles->setNumberOfParticles( setup.numberOfParticles );

   particles->setSearchRadius( setup.searchRadius );
   particles->setGridOrigin( setup.gridOrigin );
   particles->setGridDimensions( setup.gridDimensions );

   // input parameters for enhanced decomposition paritlces
   particles->setOverlapWidth( setup.overlapWidth );
   particles->setGridReferentialOrigin( setup.gridReferentialOrigin );
   particles->setGridOriginGlobalCoords( setup.gridOriginGlobalCoords );

   // assign particles
   ASSERT_TRUE( assignPoints2D( particles ) );

   // compute and test particle cell indices
   particles->computeParticleCellIndices();
   const auto cellIndices = particles->getParticleCellIndices().getConstView();

   //[ 1, 1 ]
   EXPECT_EQ( cellIndices.getElement( 10 ), 18 );
   //[ 2, 1 ]
   EXPECT_EQ( cellIndices.getElement( 4 ), 19 );
   EXPECT_EQ( cellIndices.getElement( 16 ), 19 );
   //[ 3, 1 ]
   EXPECT_EQ( cellIndices.getElement( 5 ), 20 );
   //[ 4, 1 ]
   EXPECT_EQ( cellIndices.getElement( 11 ), 21 );

   //[ 1, 2 ]
   EXPECT_EQ( cellIndices.getElement( 0 ), 26 );
   EXPECT_EQ( cellIndices.getElement( 18 ), 26 );
   //[ 2, 2 ]
   EXPECT_EQ( cellIndices.getElement( 19 ), 27 );
   //[ 3, 2 ]
   EXPECT_EQ( cellIndices.getElement( 6 ), 28 );
   EXPECT_EQ( cellIndices.getElement( 9 ), 28 );
   //[ 4, 2 ]
   EXPECT_EQ( cellIndices.getElement( 2 ), 29 );
   EXPECT_EQ( cellIndices.getElement( 14 ), 29 );
   EXPECT_EQ( cellIndices.getElement( 15 ), 29 );

   //[ 1, 3 ]
   EXPECT_EQ( cellIndices.getElement( 1 ), 34 );
   EXPECT_EQ( cellIndices.getElement( 12 ), 34 );
   //[ 2, 3 ]
   EXPECT_EQ( cellIndices.getElement( 17 ), 35 );
   EXPECT_EQ( cellIndices.getElement( 7 ), 35 );
   EXPECT_EQ( cellIndices.getElement( 8 ), 35 );
   //[ 3, 3 ]
   EXPECT_EQ( cellIndices.getElement( 13 ), 36 );
   //[ 4, 3 ]
   EXPECT_EQ( cellIndices.getElement( 3 ), 37 );


}

