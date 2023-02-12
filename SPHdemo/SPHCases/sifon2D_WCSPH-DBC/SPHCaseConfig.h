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

   static constexpr int spaceDimension = 2;

   static constexpr float mass = 0.025f;
   static constexpr float speedOfSound = 34.3f;
   static constexpr float coefB = 168070.0f;
   static constexpr float rho0 = 1000.0f;
   static constexpr float h = 0.0070711f;
   static constexpr float delta = 0.1f;
   static constexpr float alpha = 0.02f;
   static constexpr float eps = 0.001f;
   static constexpr float dtInit = 2e-05f;

   struct INLET
   {
      static constexpr float orientation_x = 1.0f;
      static constexpr float orientation_y = 0.0f;
      static constexpr float velocity_x = 0.5f;
      static constexpr float velocity_y = 0.0f;
      static constexpr float position_x = 0.185f;
      static constexpr float position_y = 0.1f;
      static constexpr float inlet_density = 1000.0f;
      static constexpr float bufferWidth_x = 0.0175f; //ie 4 layers
      static constexpr float bufferWidth_y = 0.0f; //ie 4 layers
      static constexpr float bufferEdge = 0.2225f; //ie 4 layers
   };

};

} //SPH
} //namespace ParticleSystem
} //namespace TNL

