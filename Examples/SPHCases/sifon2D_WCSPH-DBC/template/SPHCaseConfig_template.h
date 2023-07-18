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

   static constexpr int spaceDimension = placeholderDimension;

   static constexpr float mass = placeholderMassf;
   static constexpr float speedOfSound = placeholderSpeedOfSoundf;
   static constexpr float coefB = placeholderCoefBf;
   static constexpr float rho0 = placeholderDensityf;
   static constexpr float h = placeholderSmoothingLengthf;
   static constexpr float delta = 0.1f;
   static constexpr float alpha = 0.02f;
   static constexpr float eps = 0.001f;
   static constexpr float dtInit = placeholderTimeStepf;

   struct INLET
   {
      static constexpr float orientation_x = placeholderOBP1Orientation_xf;
      static constexpr float orientation_y = placeholderOBP1Orientation_yf;
      static constexpr float velocity_x = placeholderOBP1Velocity_xf;
      static constexpr float velocity_y = placeholderOBP1Velocity_yf;
      static constexpr float position_x = placeholderOBP1Position_xf;
      static constexpr float position_y = placeholderOBP1Position_yf;
      static constexpr float inlet_density = placeholderOBP1Densityf;
      static constexpr float bufferWidth_x = placeholderOBP1Width_xf; //ie 4 layers
      static constexpr float bufferWidth_y = placeholderOBP1Width_yf; //ie 4 layers
      static constexpr float bufferEdge = placeholderOBP1BufferEdgef; //ie 4 layers
   };

};

} //SPH
} //namespace ParticleSystem
} //namespace TNL

