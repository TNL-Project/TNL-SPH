#pragma once

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Device>
class SPHFluidConfig
{
  public:

  static constexpr int spaceDimension = 2;

  using DeviceType = Device;

  using GlobalIndexType = int;
  using LocalIndexType = int;
  using CellIndexType = int;
  using RealType = float;

};

} // SPH
} // ParticleSystem
} // TNL
