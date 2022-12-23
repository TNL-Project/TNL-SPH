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

   static constexpr float mass = 0.001f;
   static constexpr float speedOfSound = 34.3f;
   static constexpr float coefB = 167610.859375f;
   static constexpr float rho0 = 1000.f;
   static constexpr float h = 0.001414f;
   static constexpr float delta = 0.1f;
   static constexpr float alpha = 0.02f;
   static constexpr float eps = 0.001f;
   static constexpr float dtInit = 0.00001f;
};

} //SPH
} //namespace ParticleSystem
} //namespace TNL

