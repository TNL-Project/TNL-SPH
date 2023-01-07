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
   static constexpr float etaLimiter = 3.0f;
   static constexpr float dtInit = placeholderTimeStepf;
};

} //SPH
} //namespace ParticleSystem
} //namespace TNL

