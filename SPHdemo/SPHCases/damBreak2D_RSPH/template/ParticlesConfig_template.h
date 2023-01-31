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

   static constexpr int spaceDimension = placeholderDimension;
   static constexpr int numberOfParticles = placeholderFluidParticles;
   static constexpr int maxOfNeigborsPerParticle = 70;

   static constexpr RealType searchRadius = placeholderSearchRadius*1.001;
   static constexpr int gridXsize = placeholderGridXSize;
   static constexpr int gridYsize = placeholderGridYSize;

   //static constexpr CoordinatesType origin = {0, 0}; //.. I would like something like this
   static constexpr RealType gridXbegin = placeholderGridXBegin;
   static constexpr RealType gridYbegin = placeholderGridYBegin;

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

   static constexpr int spaceDimension = placeholderDimension;
   static constexpr int numberOfParticles = placeholderBoundaryParticles;
   static constexpr int maxOfNeigborsPerParticle = 70;

   static constexpr RealType searchRadius = placeholderSearchRadius*1.001;
   static constexpr int gridXsize = placeholderGridXSize;
   static constexpr int gridYsize = placeholderGridYSize;

   //static constexpr CoordinatesType origin = {0, 0}; //.. I would like something like this
   static constexpr RealType gridXbegin = placeholderGridXBegin;
   static constexpr RealType gridYbegin = placeholderGridYBegin;

   using CoordinatesType = Containers::StaticVector< spaceDimension, int >;
   using CellIndexerType = SimpleCellIndex< spaceDimension, ParticleSystemConfig_boundary >; //?
   using NeighborListType = typename Algorithms::Segments::Ellpack< DeviceType, int >;
};

} //namespace ParticleSystem
} //namespace TNL

