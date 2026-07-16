#pragma once

#include <gtest/gtest.h>
#include "gtest/gtest.h"
#include "gmock/gmock.h" //test vectors

#include <Particles/ParticlesLinkedList.h>

using namespace TNL;
using namespace TNL::Particles;

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
   using CoordinatesType = typename ParticlesTraitsType::CoordinatesType;
   using PointType = typename ParticlesTraitsType::PointType;

   // input parameters
   const int numberOfParticles = 20;
   const int numberOfAllocatedParticles = 28;

   const float searchRadius = 0.5;
   const PointType origin = { -searchRadius, -searchRadius };
   const CoordinatesType dimensions = { 6, 5 };

   // input parameters for enhanced decomposition paritlces
   const int overlapWidth = 1;
   const PointType referentialOrigin = { origin[ 0 ] - overlapWidth * searchRadius, origin[ 1 ]  - overlapWidth * searchRadius };
   const CoordinatesType globalOriginCoordinates = { overlapWidth, overlapWidth };
};

TEST( ParticlesWithDecomposition2DTest, ParticlesPropertiesHost )
{
   using Device = TNL::Devices::Host;
   using Setup = ParticlesWithDecomposition2DSetup< Device >;
   using Particles = TNL::Particles::ParticlesLinkedList< Setup::Config, Device >;
   using ParticlesPointer = typename Pointers::SharedPointer< Particles, Device >;
   using PointType = typename Setup::PointType;
   using CoordinatesType = typename Setup::CoordinatesType;

   Setup setup;
   ParticlesPointer particles;
   // input parameters
   particles->setSize( setup.numberOfAllocatedParticles );
   particles->setNumberOfParticles( setup.numberOfParticles );

   particles->setSearchRadius( setup.searchRadius );
   particles->setOrigin( setup.origin );
   particles->setDimensions( setup.dimensions );

   // input parameters for enhanced decomposition paritlces
   particles->setOverlapWidth( setup.overlapWidth );
   particles->setReferentialOrigin( setup.referentialOrigin );
   particles->setGlobalOriginCoordinates( setup.globalOriginCoordinates );

   // assign particles
   ASSERT_TRUE( assignPoints2D( particles ) );

   // test particle system default properties
   EXPECT_EQ( particles->getNumberOfParticles(), 20 );
   EXPECT_EQ( particles->getNumberOfAllocatedParticles(), 28 );
   EXPECT_EQ( particles->getPoints().getSize(), 28 );
   EXPECT_EQ( particles->getSortPermutations()->getSize(), 28 );
   EXPECT_EQ( particles->getParticleCellIndices().getSize(), 28 ); //cellList related

   EXPECT_EQ( particles->getSearchRadius(), 0.5 );
   const PointType origin = { -0.5, -0.5 };
   EXPECT_EQ( particles->getOrigin(), origin );
   const CoordinatesType dimensions = { 6, 5 };
   EXPECT_EQ( particles->getDimensions(), dimensions );
   EXPECT_EQ( particles->getCellFirstLastParticleList().getSize(), 56 );

   EXPECT_EQ( particles->getOverlapWidth(), 1 );
   const PointType referentialOrigin = { -1.0, -1.0 };
   EXPECT_EQ( particles->getReferentialOrigin(), referentialOrigin );
   const CoordinatesType globalOriginCoordinates = { 1, 1 };
   EXPECT_EQ( particles->getGlobalOriginCoordinates(), globalOriginCoordinates );

   // NOTE: With changing overlap width, referentialOrigin should be changed aswell
   // test particle system computed properties
   const PointType originWithOverlap = { -1.0, -1.0 };
   EXPECT_EQ( particles->getOriginWithOverlap(), originWithOverlap );
   const CoordinatesType dimensionsWithOverlap = { 8, 7 };
   EXPECT_EQ( particles->getDimensionsWithOverlap(), dimensionsWithOverlap );

   // change overlap, recompute sizes
   particles->setOverlapWidth( 3 );
   const PointType originWithOverlap_scaledUp = { -2.0, -2.0 };
   EXPECT_EQ( particles->getOriginWithOverlap(), originWithOverlap_scaledUp );
   const CoordinatesType dimensionsWithOverlap_scaledUp = { 12, 11 };
   EXPECT_EQ( particles->getDimensionsWithOverlap(), dimensionsWithOverlap_scaledUp );
   EXPECT_EQ( particles->getCellFirstLastParticleList().getSize(), 132 );

   particles->setOverlapWidth( 2 );
   const PointType originWithOverlap_scaledDown = { -1.5, -1.5 };
   EXPECT_EQ( particles->getOriginWithOverlap(), originWithOverlap_scaledDown );
   const CoordinatesType dimensionsWithOverlap_scaledDown = { 10, 9 };
   EXPECT_EQ( particles->getDimensionsWithOverlap(), dimensionsWithOverlap_scaledDown );
   EXPECT_EQ( particles->getCellFirstLastParticleList().getSize(), 90 );

   // change grid size
   particles->setDimensions( { 3, 4 } );
   const CoordinatesType dimensions_updated = { 3, 4 };
   EXPECT_EQ( particles->getDimensions(), dimensions_updated );
   const CoordinatesType dimensionsWithOverlap_updated = { 7, 8 };
   EXPECT_EQ( particles->getDimensionsWithOverlap(), dimensionsWithOverlap_updated );
   EXPECT_EQ( particles->getCellFirstLastParticleList().getSize(), 56 );
}

TEST( ParticlesWithDecomposition2DTest, ParticlesPropertiesDevice )
{
   using Device = TNL::Devices::Cuda;
   using Setup = ParticlesWithDecomposition2DSetup< Device >;
   using Particles = TNL::Particles::ParticlesLinkedList< Setup::Config, Device >;
   using ParticlesPointer = typename Pointers::SharedPointer< Particles, Device >;
   using PointType = typename Setup::PointType;
   using CoordinatesType = typename Setup::CoordinatesType;

   Setup setup;
   ParticlesPointer particles;
   // input parameters
   particles->setSize( setup.numberOfAllocatedParticles );
   particles->setNumberOfParticles( setup.numberOfParticles );

   particles->setSearchRadius( setup.searchRadius );
   particles->setOrigin( setup.origin );
   particles->setDimensions( setup.dimensions );

   // input parameters for enhanced decomposition paritlces
   particles->setOverlapWidth( setup.overlapWidth );
   particles->setReferentialOrigin( setup.referentialOrigin );
   particles->setGlobalOriginCoordinates( setup.globalOriginCoordinates );

   // assign particles
   ASSERT_TRUE( assignPoints2D( particles ) );

   // test particle system default properties
   EXPECT_EQ( particles->getNumberOfParticles(), 20 );
   EXPECT_EQ( particles->getNumberOfAllocatedParticles(), 28 );
   EXPECT_EQ( particles->getPoints().getSize(), 28 );
   EXPECT_EQ( particles->getSortPermutations()->getSize(), 28 );
   EXPECT_EQ( particles->getParticleCellIndices().getSize(), 28 ); //cellList related

   EXPECT_EQ( particles->getSearchRadius(), 0.5 );
   const PointType origin = { -0.5, -0.5 };
   EXPECT_EQ( particles->getOrigin(), origin );
   const CoordinatesType dimensions = { 6, 5 };
   EXPECT_EQ( particles->getDimensions(), dimensions );
   EXPECT_EQ( particles->getCellFirstLastParticleList().getSize(), 56 );

   EXPECT_EQ( particles->getOverlapWidth(), 1 );
   const PointType referentialOrigin = { -1.0, -1.0 };
   EXPECT_EQ( particles->getReferentialOrigin(), referentialOrigin );
   const CoordinatesType globalOriginCoordinates = { 1, 1 };
   EXPECT_EQ( particles->getGlobalOriginCoordinates(), globalOriginCoordinates );

   // test particle system computed properties
   const PointType originWithOverlap = { -1.0, -1.0 };
   EXPECT_EQ( particles->getOriginWithOverlap(), originWithOverlap );
   const CoordinatesType dimensionsWithOverlap = { 8, 7 };
   EXPECT_EQ( particles->getDimensionsWithOverlap(), dimensionsWithOverlap );

   // NOTE: With changing overlap width, referentialOrigin should be changed aswell
   // change overlap, recompute sizes
   particles->setOverlapWidth( 3 );
   const PointType originWithOverlap_scaledUp = { -2.0, -2.0 };
   EXPECT_EQ( particles->getOriginWithOverlap(), originWithOverlap_scaledUp );
   const CoordinatesType dimensionsWithOverlap_scaledUp = { 12, 11 };
   EXPECT_EQ( particles->getDimensionsWithOverlap(), dimensionsWithOverlap_scaledUp );
   EXPECT_EQ( particles->getCellFirstLastParticleList().getSize(), 132 );

   particles->setOverlapWidth( 2 );
   const PointType originWithOverlap_scaledDown = { -1.5, -1.5 };
   EXPECT_EQ( particles->getOriginWithOverlap(), originWithOverlap_scaledDown );
   const CoordinatesType dimensionsWithOverlap_scaledDown = { 10, 9 };
   EXPECT_EQ( particles->getDimensionsWithOverlap(), dimensionsWithOverlap_scaledDown );
   EXPECT_EQ( particles->getCellFirstLastParticleList().getSize(), 90 );

   // change grid size
   particles->setDimensions( { 3, 4 } );
   const CoordinatesType dimensions_updated = { 3, 4 };
   EXPECT_EQ( particles->getDimensions(), dimensions_updated );
   const CoordinatesType dimensionsWithOverlap_updated = { 7, 8 };
   EXPECT_EQ( particles->getDimensionsWithOverlap(), dimensionsWithOverlap_updated );
   EXPECT_EQ( particles->getCellFirstLastParticleList().getSize(), 56 );
}

TEST( ParticlesWithDecomposition2DTest, ComputeParticleCellIndicesCuda )
{
   using Device = TNL::Devices::Cuda;
   using Setup = ParticlesWithDecomposition2DSetup< Device >;
   using Particles = TNL::Particles::ParticlesLinkedList< Setup::Config, Device >;
   using ParticlesPointer = typename Pointers::SharedPointer< Particles, Device >;
   using PointType = typename Setup::PointType;
   using CoordinatesType = typename Setup::CoordinatesType;

   Setup setup;
   ParticlesPointer particles;
   // input parameters
   particles->setSize( setup.numberOfAllocatedParticles );
   particles->setNumberOfParticles( setup.numberOfParticles );

   particles->setSearchRadius( setup.searchRadius );
   particles->setOrigin( setup.origin );
   particles->setDimensions( setup.dimensions );

   // input parameters for enhanced decomposition paritlces
   particles->setOverlapWidth( setup.overlapWidth );
   particles->setReferentialOrigin( setup.referentialOrigin );
   particles->setGlobalOriginCoordinates( setup.globalOriginCoordinates );

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

   // update  overlap size
   particles->setOverlapWidth( 2 );
   // NOTE: With changing overlap width, referentialOrigin should be changed aswell
   const PointType shiftReferentialOrigin_overlap2 = particles->getSearchRadius() * particles->getOverlapWidth();
   particles->setReferentialOrigin( particles->getOrigin() - shiftReferentialOrigin_overlap2 );

   particles->computeParticleCellIndices();
   const auto cellIndices_resizedOverlap = particles->getParticleCellIndices().getConstView();

   //[ 1, 1 ]
   EXPECT_EQ( cellIndices_resizedOverlap.getElement( 10 ), 44 );
   //[ 2, 1 ]
   EXPECT_EQ( cellIndices_resizedOverlap.getElement( 4 ), 45 );
   EXPECT_EQ( cellIndices_resizedOverlap.getElement( 16 ), 45 );
   //[ 3, 1 ]
   EXPECT_EQ( cellIndices_resizedOverlap.getElement( 5 ), 46 );
   //[ 4, 1 ]
   EXPECT_EQ( cellIndices_resizedOverlap.getElement( 11 ), 47 );

   //[ 1, 2 ]
   EXPECT_EQ( cellIndices_resizedOverlap.getElement( 0 ), 54 );
   EXPECT_EQ( cellIndices_resizedOverlap.getElement( 18 ), 54 );
   //[ 2, 2 ]
   EXPECT_EQ( cellIndices_resizedOverlap.getElement( 19 ), 55 );
   //[ 3, 2 ]
   EXPECT_EQ( cellIndices_resizedOverlap.getElement( 6 ), 56 );
   EXPECT_EQ( cellIndices_resizedOverlap.getElement( 9 ), 56 );
   //[ 4, 2 ]
   EXPECT_EQ( cellIndices_resizedOverlap.getElement( 2 ), 57 );
   EXPECT_EQ( cellIndices_resizedOverlap.getElement( 14 ), 57 );
   EXPECT_EQ( cellIndices_resizedOverlap.getElement( 15 ), 57 );

   //[ 1, 3 ]
   EXPECT_EQ( cellIndices_resizedOverlap.getElement( 1 ), 64 );
   EXPECT_EQ( cellIndices_resizedOverlap.getElement( 12 ), 64 );
   //[ 2, 3 ]
   EXPECT_EQ( cellIndices_resizedOverlap.getElement( 17 ), 65 );
   EXPECT_EQ( cellIndices_resizedOverlap.getElement( 7 ), 65 );
   EXPECT_EQ( cellIndices_resizedOverlap.getElement( 8 ), 65 );
   //[ 3, 3 ]
   EXPECT_EQ( cellIndices_resizedOverlap.getElement( 13 ), 66 );
   //[ 4, 3 ]
   EXPECT_EQ( cellIndices_resizedOverlap.getElement( 3 ), 67 );

}

