#pragma once

#include <TNL/Containers/StaticVector.h>
#include <TNL/Algorithms/Segments/CSR.h>
#include <TNL/Algorithms/Segments/Ellpack.h>

#include "../../Particles/GenerateCellIndex.h"

namespace TNL {
namespace ParticleSystem {

class ParticleSystemConfig
{
  public:

  using GlobalIndexType = int;
  using LocalIndexType = int;
  using CellIndexType = int;
  using RealType = float;

  static constexpr int spaceDimension = 2;
  //static constexpr int numberOfParticles = 21001;
  static constexpr int numberOfParticles = 1961;
  static constexpr int maxOfNeigborsPerParticle = 150;

  using CoordinatesType = Containers::StaticVector< 2, int >;

  static constexpr RealType searchRadius = 0.07071068;
  static constexpr int gridXsize = 60; /* 20 */
  static constexpr int gridYsize = 60; /* 20 */

  //static constexpr CoordinatesType origin = {0, 0}; //.. I would like something like this
  /*
  static constexpr int gridXbegin = -1;
  static constexpr int gridYbegin = -1;
  */

  static constexpr RealType gridXbegin = -0.1;
  static constexpr RealType gridYbegin = -0.1;

  // ... set particle system ...

  using DeviceType = Devices::Host;
  using CellIndexerType = SimpleCellIndex<ParticleSystemConfig, DeviceType>; //?


  using NeighborListType = typename Algorithms::Segments::Ellpack< TNL::Devices::Host, int >;
};

} //namespace ParticleSystem
} //namespace TNL
