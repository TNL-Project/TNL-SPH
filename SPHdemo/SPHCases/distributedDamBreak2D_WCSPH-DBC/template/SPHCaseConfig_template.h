#pragma once

#include "../../../SPH/SPHTraits.h"
#include <limits>

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


} //SPH
} //namespace ParticleSystem
} //namespace TNL

