#pragma once

#include <TNL/Containers/StaticVector.h>
#include <TNL/Algorithms/Segments/CSR.h>
#include <TNL/Algorithms/Segments/Ellpack.h>

#include "../../../Particles/GenerateCellIndex.h"

namespace TNL {
namespace ParticleSystem {

template< typename Device >
class ParticleSystemConfig
{
   public:
   using DeviceType = Device;

   using GlobalIndexType = int;
   using LocalIndexType = int;
   using CellIndexType = int;
   using RealType = float;

   static constexpr int spaceDimension = 3;
   static constexpr int numberOfParticles = 79380;
   static constexpr int numberOfAllocatedParticles = 79380;
   static constexpr int maxOfNeigborsPerParticle = 70;

   static constexpr RealType searchRadius = 0.08*1.001;
   static constexpr int gridXsize = 44;
   static constexpr int gridYsize = 16;
   static constexpr int gridZsize = 21;

   //static constexpr CoordinatesType origin = {0, 0}; //.. I would like something like this
   static constexpr RealType gridXbegin = -0.1202f;
   static constexpr RealType gridYbegin = -0.1202f;
   static constexpr RealType gridZbegin = -0.1202f;

   using CoordinatesType = Containers::StaticVector< spaceDimension, int >;
   using CellIndexerType = SimpleCellIndex< spaceDimension, ParticleSystemConfig >; //?
   using NeighborListType = typename Algorithms::Segments::Ellpack< DeviceType, int >;
};

template< typename Device >
class ParticleSystemConfig_boundary
{
   public:
   using DeviceType = Device;

   using GlobalIndexType = int;
   using LocalIndexType = int;
   using CellIndexType = int;
   using RealType = float;

   static constexpr int spaceDimension = 3;
   static constexpr int numberOfParticles = 93042;
   static constexpr int numberOfAllocatedParticles = 93042;
   static constexpr int maxOfNeigborsPerParticle = 70;

   static constexpr RealType searchRadius = 0.08*1.001;
   static constexpr int gridXsize = 44;
   static constexpr int gridYsize = 16;
   static constexpr int gridZsize = 21;

   //static constexpr CoordinatesType origin = {0, 0}; //.. I would like something like this
   static constexpr RealType gridXbegin = -0.1202f;
   static constexpr RealType gridYbegin = -0.1202f;
   static constexpr RealType gridZbegin = -0.1202f;

   using CoordinatesType = Containers::StaticVector< spaceDimension, int >;
   using CellIndexerType = SimpleCellIndex< spaceDimension, ParticleSystemConfig_boundary >;
   using NeighborListType = typename Algorithms::Segments::Ellpack< DeviceType, int >;
};

} //namespace ParticleSystem
} //namespace TNL

