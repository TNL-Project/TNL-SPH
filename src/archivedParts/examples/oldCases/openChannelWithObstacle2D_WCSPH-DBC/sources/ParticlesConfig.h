#pragma once

#include <TNL/Containers/StaticVector.h>
#include <TNL/Algorithms/Segments/CSR.h>
#include <TNL/Algorithms/Segments/Ellpack.h>

#include <Particles/GenerateCellIndexFloating.h>
#include <Particles/ParticlesTraits.h>

namespace TNL {
namespace ParticleSystem {
namespace ParticleSystemConfig {

/**
 * TYPES OF PARTICLE SYSTEM AND NEIGHBOR SEARCH (necessary)
 *
 * This class is used to store parameters necessary for particle system,
 * i.e. data types for quantities and indices. It also defines dimension
 * and attributes of particle system.
 *
 * It is necessary to enter TYPES for:
 * - GlobalIndexType
 * - LocalIndexType
 * - CellIndexType
 * - RealType
 *
 * - CoordinatesType
 * - CellIndexerType - algorithm to index particle for neighborSearch
 * - NeighborListType - NOT USED ATM
 */
template< typename Device >
class ParticleSystemConfig
{
   public:
   using DeviceType = Device;

   using GlobalIndexType = int;
   using LocalIndexType = int;
   using CellIndexType = int;
   using RealType = float;

   static constexpr int spaceDimension = 2;

   using CoordinatesType = Containers::StaticVector< spaceDimension, int >;
   using CellIndexerType = SimpleCellIndex< spaceDimension, ParticleSystemConfig, std::index_sequence< 1, 0 > >;
   using NeighborListType = typename Algorithms::Segments::Ellpack< DeviceType, int >; //deprecated
};

/**
 * PARAMETERS OF PARTICLE SYSTEM AND NEIGHBOR SEARCH (necessary)
 *
 * This class is used to store parameters necessary for particle system,
 * neighbor search and SPH simulation initialization. This is basically
 * just placeholder to load the data.
 *
 * It is necessary to enter:
 * - Number of loaded particles and nuber of allocated particles for
 *   all "objects" (i.e. fluid, boundary, inlet buffer, etc.) of the
 *   SPH simulation.
 * -- numberOfParticles - number of loaded fluid particles
 * -- numberOfAllocatedParticles - number of initially allocated fluid
 * -- numberOfBoundaryParticles - number of loaded bound particles
 * -- numberOfAllocatedBoundaryParticles - number of initially allocated bound
 *
 * - searchRadius - define the size of grid cell
 * - gridXsize - define the size of grid in dim X
 * - gridYsize - define the size of grid in dim Y
 * - gridXorigin - define the origin of grid in dim X
 * - gridYorigin - define the origin of grid in dim Y
 */
template< typename Device >
class ParticleInitialSetup
{
   public:
   using ParticlesConfig = ParticleSystemConfig< Device >;
   using ParticlesTraitsType = ParticlesTraits< ParticlesConfig, Device >;
   using IndexVectorType = typename ParticlesTraitsType::IndexVectorType;
   using PointType = typename ParticlesTraitsType::PointType;

   const int numberOfParticles = 15000;
   const int numberOfAllocatedParticles = 50000;
   const int numberOfBoundaryParticles = 1512;
   const int numberOfAllocatedBoundaryParticles = 1512;

   int numberOfOpenBoundaryParticles[ 2 ] = { 200,  200 };
   int numberOfAllocatedOpenBoundaryParticles[ 2 ] = { 200, 600 };

   const float searchRadius = 0.0056568 * 1.001;
   const int gridXsize = 181;
   const int gridYsize = 12;
   const PointType gridOrigin = { -0.0096768f, -0.0096768f };

   const IndexVectorType gridSize = { gridXsize, gridYsize };
   const int numberOfGridCells = gridXsize * gridYsize;
};


} //namespace ParticleSystemConfig
} //namespace ParticleSystem
} //namespace TNL

