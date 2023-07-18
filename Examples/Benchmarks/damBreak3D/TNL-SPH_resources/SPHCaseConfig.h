#pragma once

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Device>
class SPHCaseConfig
{
   public:
   using DeviceType = Device;

   using GlobalIndexType = int;
   using LocalIndexType = int;
   using CellIndexType = int;
   using RealType = float;

   static constexpr int spaceDimension = 3;

   static constexpr float mass = 0.008f;
   static constexpr float speedOfSound = 45.17f;
   static constexpr float coefB = 291475.6f;
   static constexpr float rho0 = 1000.0f;
   static constexpr float h = 0.04f;
   static constexpr float delta = 0.1f;
   static constexpr float alpha = 0.02f;
   static constexpr float eps = 0.001f;
   static constexpr float dtInit = 0.00013283f;
};

} //SPH
} //namespace ParticleSystem
} //namespace TNL

