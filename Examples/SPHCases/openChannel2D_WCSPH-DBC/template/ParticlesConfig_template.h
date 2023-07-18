#pragma once

#include <TNL/Containers/StaticVector.h>
#include <TNL/Algorithms/Segments/CSR.h>
#include <TNL/Algorithms/Segments/Ellpack.h>

#include "../../../Particles/GenerateCellIndex.h"
#include "../../../Particles/ParticlesTraits.h"

namespace TNL {
namespace ParticleSystem {
namespace ParticleSystemConfig {

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
class ParticleInitialSetup
{
   public:

   static constexpr int numberOfParticles = placeholderFluidParticles;
   static constexpr int numberOfAllocatedParticles = placeholderAllocatedFluidParticles;
   static constexpr int numberOfBoundaryParticles = placeholderBoundaryParticles;
   static constexpr int numberOfAllocatedBoundaryParticles = placeholderBoundaryParticles;

   static constexpr int numberOfInletParticles = placeholderInletParticles;
   static constexpr int numberOfAllocatedInletParticles = placeholderAllocatedInletParticles;
   static constexpr int numberOfOutletParticles = placeholderOutletParticles;
   static constexpr int numberOfAllocatedOutletParticles = placeholderAllocatedOutletParticles;

   static constexpr float searchRadius = placeholderSearchRadius * 1.001;
   static constexpr int gridXsize = placeholderGridXSize;
   static constexpr int gridYsize = placeholderGridYSize;
   static constexpr float gridXorigin = placeholderGridXBeginf;
   static constexpr float gridYorigin = placeholderGridYBeginf;
};

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

   static constexpr int spaceDimension = placeholderDimension;

   using CoordinatesType = Containers::StaticVector< spaceDimension, int >;
   using CellIndexerType = SimpleCellIndex< spaceDimension, ParticleSystemConfig >;
   using NeighborListType = typename Algorithms::Segments::Ellpack< DeviceType, int >; //deprecated
};


//template< typename Device >
//class ParticleSystemConfig
//{
//   public:
//   using DeviceType = Device;
//
//   using GlobalIndexType = int;
//   using LocalIndexType = int;
//   using CellIndexType = int;
//   using RealType = float;
//
//   static constexpr int spaceDimension = placeholderDimension;
//   static constexpr int numberOfParticles = placeholderFluidParticles;
//   static constexpr int numberOfAllocatedParticles = placeholderAllocatedFluidParticles;
//   static constexpr float reallocationCoef = 1.5f;
//   static constexpr int maxOfNeigborsPerParticle = 70;
//
//   static constexpr RealType searchRadius = placeholderSearchRadius*1.001;
//   static constexpr int gridXsize = placeholderGridXSize;
//   static constexpr int gridYsize = placeholderGridYSize;
//
//   //static constexpr CoordinatesType origin = {0, 0}; //.. I would like something like this
//   static constexpr RealType gridXbegin = placeholderGridXBegin;
//   static constexpr RealType gridYbegin = placeholderGridYBegin;
//
//   using CoordinatesType = Containers::StaticVector< spaceDimension, int >;
//   using CellIndexerType = SimpleCellIndex< spaceDimension, ParticleSystemConfig >; //?
//   using NeighborListType = typename Algorithms::Segments::Ellpack< DeviceType, int >;
//};
//
//template< typename Device >
//class ParticleSystemConfig_boundary
//{
//   public:
//   using DeviceType = Device;
//
//   using GlobalIndexType = int;
//   using LocalIndexType = int;
//   using CellIndexType = int;
//   using RealType = float;
//
//   static constexpr int spaceDimension = placeholderDimension;
//   static constexpr int numberOfParticles = placeholderBoundaryParticles;
//   static constexpr int numberOfAllocatedParticles = placeholderBoundaryParticles;
//   static constexpr int maxOfNeigborsPerParticle = 70;
//
//   static constexpr RealType searchRadius = placeholderSearchRadius*1.001;
//   static constexpr int gridXsize = placeholderGridXSize;
//   static constexpr int gridYsize = placeholderGridYSize;
//
//   //static constexpr CoordinatesType origin = {0, 0}; //.. I would like something like this
//   static constexpr RealType gridXbegin = placeholderGridXBegin;
//   static constexpr RealType gridYbegin = placeholderGridYBegin;
//
//   using CoordinatesType = Containers::StaticVector< spaceDimension, int >;
//   using CellIndexerType = SimpleCellIndex< spaceDimension, ParticleSystemConfig_boundary >; //?
//   using NeighborListType = typename Algorithms::Segments::Ellpack< DeviceType, int >;
//};
//
//template< typename Device >
//class ParticleSystemConfig_inletBuffer
//{
//   public:
//   using DeviceType = Device;
//
//   using GlobalIndexType = int;
//   using LocalIndexType = int;
//   using CellIndexType = int;
//   using RealType = float;
//
//   static constexpr int spaceDimension = placeholderDimension;
//   static constexpr int numberOfParticles = placeholderBufferParticles;
//   static constexpr int numberOfAllocatedParticles = placeholderAllocatedBufferParticles;
//   static constexpr int maxOfNeigborsPerParticle = 70;
//
//   static constexpr RealType searchRadius = placeholderSearchRadius*1.001;
//   static constexpr int gridXsize = placeholderGridXSize;
//   static constexpr int gridYsize = placeholderGridYSize;
//
//   //static constexpr CoordinatesType origin = {0, 0}; //.. I would like something like this
//   static constexpr RealType gridXbegin = placeholderGridXBegin;
//   static constexpr RealType gridYbegin = placeholderGridYBegin;
//
//   using CoordinatesType = Containers::StaticVector< spaceDimension, int >;
//   using CellIndexerType = SimpleCellIndex< spaceDimension, ParticleSystemConfig_inletBuffer >;
//   using NeighborListType = typename Algorithms::Segments::Ellpack< DeviceType, int >;
//};
//
//template< typename Device >
//class ParticleSystemConfig_inlet2Buffer
//{
//   public:
//   using DeviceType = Device;
//
//   using GlobalIndexType = int;
//   using LocalIndexType = int;
//   using CellIndexType = int;
//   using RealType = float;
//
//   static constexpr int spaceDimension = placeholderDimension;
//   static constexpr int numberOfParticles = placeholderBuffer2Particles;
//   static constexpr int numberOfAllocatedParticles = placeholderAllocatedBuffer2Particles;
//   static constexpr int maxOfNeigborsPerParticle = 70;
//
//   static constexpr RealType searchRadius = placeholderSearchRadius*1.001;
//   static constexpr int gridXsize = placeholderGridXSize;
//   static constexpr int gridYsize = placeholderGridYSize;
//
//   //static constexpr CoordinatesType origin = {0, 0}; //.. I would like something like this
//   static constexpr RealType gridXbegin = placeholderGridXBegin;
//   static constexpr RealType gridYbegin = placeholderGridYBegin;
//
//   using CoordinatesType = Containers::StaticVector< spaceDimension, int >;
//   using CellIndexerType = SimpleCellIndex< spaceDimension, ParticleSystemConfig_inlet2Buffer >;
//   using NeighborListType = typename Algorithms::Segments::Ellpack< DeviceType, int >;
//};

} //namespace ParticleSystemConfig
} //namespace ParticleSystem
} //namespace TNL

