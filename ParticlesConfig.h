#pragma once

#include <TNL/Containers/StaticVector.h>
#include <TNL/Algorithms/Segments/CSR.h>
#include <TNL/Algorithms/Segments/Ellpack.h>

#include "Particles/GenerateCellIndex.h"

namespace TNL {
namespace ParticleSystem {

class ParticleSystemConfig
{
  public:

  static constexpr int spaceDimension = 2;
  static constexpr int maxOfNeigborsPerParticle = 35;

  using CoordinatesType = Containers::StaticVector< 2, int >;

  static constexpr int searchRadius = 1;
  static constexpr int gridXsize = 8;
  static constexpr int gridYsize = 8;

  //static constexpr CoordinatesType origin = {0, 0}; //.. I would like something like this
  static constexpr int gridXbegin = 0;
  static constexpr int gridYbegin = 0;

  // ... set particle system ...

  using DeviceType = Devices::Cuda;
  using CellIndexerType = SimpleCellIndex<ParticleSystemConfig, DeviceType>; //?

  using GlobalIndexType = int;
  using LocalIndexType = int;
  using CellIndexType = int;
  using RealType = float;

  using NeighborListType = typename Algorithms::Segments::Ellpack< DeviceType, int >;
};

} //namespace ParticleSystem
} //namespace TNL

