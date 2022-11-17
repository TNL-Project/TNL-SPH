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
  using RealType = double; //float

  static constexpr float mass = 0.1;
  static constexpr float speedOfSound = 34.3;
  static constexpr float coefB = 168070;
  static constexpr float rho0 = 1000.;
  static constexpr float h = 0.01414213;
  static constexpr float delta = 0.1;
  static constexpr float alpha = 0.02;
  static constexpr float eps = 0.001;

};

} //SPH
} //namespace ParticleSystem
} //namespace TNL

