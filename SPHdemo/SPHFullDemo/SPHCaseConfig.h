#pragma once

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Device>
class SPHCaseConfig
{
  public:
  using DeviceType = Device;

  static constexpr int spaceDimension = 2;

  using GlobalIndexType = int;
  using LocalIndexType = int;
  using CellIndexType = int;
  using RealType = float;

  static constexpr float mass = 0.625;
  static constexpr float speedOfSound = 88.033;
  static constexpr float coefB = 1107129;
  static constexpr float rho0 = 1000.;
  static constexpr float h = 0.035355;

};

} //SPH
} //namespace ParticleSystem
} //namespace TNL
