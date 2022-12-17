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

   static constexpr float mass = 0.000250;
   static constexpr float speedOfSound = 34.3;
   static constexpr float coefB = 167891.140625;
   static constexpr float rho0 = 1000.;
   static constexpr float h = 0.000707;
   static constexpr float delta = 0.1;
   static constexpr float alpha = 0.02;
   static constexpr float eps = 0.001;
   static constexpr float dtInit = 0.000005;
};

} //SPH
} //namespace ParticleSystem
} //namespace TNL

