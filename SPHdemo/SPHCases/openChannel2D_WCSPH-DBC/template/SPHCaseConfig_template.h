#pragma once

namespace TNL {
namespace ParticleSystem {
namespace SPH {

/**
 * PARAMETERS OF SPH SIMULATION AND SCHEMES (necessary)
 *
 * This class is used to store core parameters for simulation control and
 * simulation initialization. This includes variables such as paths for loading
 * and saving files or the length of the simulation and the frequency of saving outputs.
 *
 */
template< typename Device>
class SPHCaseConfig
{
   public:
   using DeviceType = Device;

   /**
    * Definition of basics data types for variables and indices.
    */
   using GlobalIndexType = int;
   using LocalIndexType = int;
   using CellIndexType = int;
   using RealType = float;

   /**
    * Define the space dimension of the problem.
    */
   static constexpr int spaceDimension = placeholderDimension;

   /**
    * Define SPH parameters connected to the resolution.
    * - h - smoothing length [m]
    * - dp - initial particle distance [m]
    */
   static constexpr float dp = placeholderInitParticleDistancef;
   static constexpr float h = placeholderSmoothingLengthf;

   /**
    * Define Basics SPH constants.
    * - mass - particle mass [kg]
    * - speedOfSound - numerical speed of sound used in the SPH calculations [m/s]
    * - coefB - coefficient of the Tait equation of state coefB = c^2 * rho0 / gamma
    * - rho0 - referential density of the fluid [kg/m^3]
    */
   static constexpr float mass = placeholderMassf;
   static constexpr float speedOfSound = placeholderSpeedOfSoundf;
   static constexpr float coefB = placeholderCoefBf;
   static constexpr float rho0 = placeholderDensityf;

   /**
    * Define coefficient of diffusive term (DT), [-].
    */
   static constexpr float delta = 0.1f;

   /**
    * Define coefficient of artificial viscosity.
    */
   static constexpr float alpha = 0.02f;

   /**
    * Define initial timestep [s].
    */
   static constexpr float dtInit = placeholderTimeStepf;

   /**
    * Define external forces [m^2 / s].
    */
   static constexpr float gravity[ 2 ] { 0.f, -9.81f };

   /**
    * Define constant to prevent zero in denominator [-].
    */
   static constexpr float eps = 0.001f;
};

//template< typename Device>
//class SPHCaseConfig
//{
//   public:
//   using DeviceType = Device;
//
//   using GlobalIndexType = int;
//   using LocalIndexType = int;
//   using CellIndexType = int;
//   using RealType = float;
//
//   static constexpr int spaceDimension = placeholderDimension;
//
//   static constexpr float mass = placeholderMassf;
//   static constexpr float speedOfSound = placeholderSpeedOfSoundf;
//   static constexpr float coefB = placeholderCoefBf;
//   static constexpr float rho0 = placeholderDensityf;
//   static constexpr float h = placeholderSmoothingLengthf;
//   static constexpr float delta = 0.1f;
//   static constexpr float alpha = 0.02f;
//   static constexpr float eps = 0.001f;
//   static constexpr float dtInit = placeholderTimeStepf;
//
//   struct INLET
//   {
//      static constexpr float orientation_x = placeholderOBP1Orientation_xf;
//      static constexpr float orientation_y = placeholderOBP1Orientation_yf;
//      static constexpr float velocity_x = placeholderOBP1Velocity_xf;
//      static constexpr float velocity_y = placeholderOBP1Velocity_yf;
//      static constexpr float position_x = placeholderOBP1Position_xf;
//      static constexpr float position_y = placeholderOBP1Position_yf;
//      static constexpr float inlet_density = placeholderOBP1Densityf;
//      static constexpr float bufferWidth_x = placeholderOBP1Width_xf; //ie 4 layers
//      static constexpr float bufferWidth_y = placeholderOBP1Width_yf; //ie 4 layers
//      static constexpr float bufferEdge = placeholderOBP1BufferEdgef; //ie 4 layers
//   };
//
//   struct INLET2
//   {
//      static constexpr float orientation_x = placeholderOBP2Orientation_xf;
//      static constexpr float orientation_y = placeholderOBP2Orientation_yf;
//      static constexpr float velocity_x = placeholderOBP2Velocity_xf;
//      static constexpr float velocity_y = placeholderOBP2Velocity_yf;
//      static constexpr float position_x = placeholderOBP2Position_xf;
//      static constexpr float position_y = placeholderOBP2Position_yf;
//      static constexpr float inlet_density = placeholderOBP2Densityf;
//      static constexpr float bufferWidth_x = placeholderOBP2Width_xf; //ie 4 layers
//      static constexpr float bufferWidth_y = placeholderOBP2Width_yf; //ie 4 layers
//      static constexpr float bufferEdge = placeholderOBP2BufferEdgef; //ie 4 layers
//   };
//};

} //SPH
} //namespace ParticleSystem
} //namespace TNL

