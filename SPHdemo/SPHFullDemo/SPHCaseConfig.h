#pragma once

namespace TNL {
namespace ParticleSystem {
namespace SPH {

class SPHCaseConfig
{
  public:

  static constexpr float mass = 0.625;
  static constexpr float speedOfSound = 88.033;
  static constexpr float coefB = 1107129;
  static constexpr float rho0 = 1000.;
  static constexpr float h = 0.035355;

};

} //SPH
} //namespace ParticleSystem
} //namespace TNL
