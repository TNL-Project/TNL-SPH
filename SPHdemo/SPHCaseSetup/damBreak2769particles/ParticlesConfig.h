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
   static constexpr int numberOfParticles = 2769;
   static constexpr int maxOfNeigborsPerParticle = 70;

   static constexpr RealType searchRadius = 0.02828426*1.005;
   static constexpr int gridXsize = 65 + 2;
   static constexpr int gridYsize = 35 + 2;

   //static constexpr CoordinatesType origin = {0, 0}; //.. I would like something like this
   static constexpr RealType gridXbegin = -0.04 - searchRadius * 1;
   static constexpr RealType gridYbegin = -0.04 - searchRadius * 1;

   using CoordinatesType = Containers::StaticVector< spaceDimension, int >;
   //using CellIndexerType = SimpleCellIndex<ParticleSystemConfig, DeviceType>; //?
   using CellIndexerType = ZOrderCurve<ParticleSystemConfig, DeviceType>; //?
   using NeighborListType = typename Algorithms::Segments::Ellpack< DeviceType, int >;
};

} //namespace ParticleSystem
} //namespace TNL

