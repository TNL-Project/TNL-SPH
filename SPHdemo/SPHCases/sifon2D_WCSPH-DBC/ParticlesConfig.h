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

   static constexpr int spaceDimension = 2;
   static constexpr int numberOfParticles = 0;
   static constexpr int numberOfAllocatedParticles = 80000;
   static constexpr float reallocationCoef = 1.5f;
   static constexpr int maxOfNeigborsPerParticle = 70;

   static constexpr RealType searchRadius = 0.0141422*1.001;
   static constexpr int gridXsize = 309;
   static constexpr int gridYsize = 49;

   //static constexpr CoordinatesType origin = {0, 0}; //.. I would like something like this
   static constexpr RealType gridXbegin = -0.0241922;
   static constexpr RealType gridYbegin = -0.0241922;

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

   static constexpr int spaceDimension = 2;
   static constexpr int numberOfParticles = 5645;
   static constexpr int numberOfAllocatedParticles = 5645;
   static constexpr int maxOfNeigborsPerParticle = 70;

   static constexpr RealType searchRadius = 0.0141422*1.001;
   static constexpr int gridXsize = 309;
   static constexpr int gridYsize = 49;

   //static constexpr CoordinatesType origin = {0, 0}; //.. I would like something like this
   static constexpr RealType gridXbegin = -0.0241922;
   static constexpr RealType gridYbegin = -0.0241922;

   using CoordinatesType = Containers::StaticVector< spaceDimension, int >;
   using CellIndexerType = SimpleCellIndex< spaceDimension, ParticleSystemConfig_boundary >; //?
   using NeighborListType = typename Algorithms::Segments::Ellpack< DeviceType, int >;
};

template< typename Device >
class ParticleSystemConfig_inletBuffer
{
   public:
   using DeviceType = Device;

   using GlobalIndexType = int;
   using LocalIndexType = int;
   using CellIndexType = int;
   using RealType = float;

   static constexpr int spaceDimension = 2;
   static constexpr int numberOfParticles = 72;
   static constexpr int numberOfAllocatedParticles = 72;
   static constexpr int maxOfNeigborsPerParticle = 70;

   static constexpr RealType searchRadius = 0.0141422*1.001;
   static constexpr int gridXsize = 309;
   static constexpr int gridYsize = 49;

   //static constexpr CoordinatesType origin = {0, 0}; //.. I would like something like this
   static constexpr RealType gridXbegin = -0.0241922;
   static constexpr RealType gridYbegin = -0.0241922;

   using CoordinatesType = Containers::StaticVector< spaceDimension, int >;
   using CellIndexerType = SimpleCellIndex< spaceDimension, ParticleSystemConfig_inletBuffer >;
   using NeighborListType = typename Algorithms::Segments::Ellpack< DeviceType, int >;
};

template< typename Device >
class ParticleSystemConfig_inlet2Buffer
{
   public:
   using DeviceType = Device;

   using GlobalIndexType = int;
   using LocalIndexType = int;
   using CellIndexType = int;
   using RealType = float;

   static constexpr int spaceDimension = 2;
   static constexpr int numberOfParticles = 120;
   static constexpr int numberOfAllocatedParticles = 120;
   static constexpr int maxOfNeigborsPerParticle = 70;

   static constexpr RealType searchRadius = 0.0141422*1.001;
   static constexpr int gridXsize = 309;
   static constexpr int gridYsize = 49;

   //static constexpr CoordinatesType origin = {0, 0}; //.. I would like something like this
   static constexpr RealType gridXbegin = -0.0241922;
   static constexpr RealType gridYbegin = -0.0241922;

   using CoordinatesType = Containers::StaticVector< spaceDimension, int >;
   using CellIndexerType = SimpleCellIndex< spaceDimension, ParticleSystemConfig_inlet2Buffer >;
   using NeighborListType = typename Algorithms::Segments::Ellpack< DeviceType, int >;
};

} //namespace ParticleSystem
} //namespace TNL

