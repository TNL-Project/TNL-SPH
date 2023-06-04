#pragma once

#include "../../../SPH/Models/EquationOfState.h"
#include "../../../SPH/Models/DiffusiveTerms.h"
#include "../../../SPH/Models/VisousTerms.h"

#include "../../../SPH/SPHTraits.h"
#include <TNL/Devices/Cuda.h>
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
class SPHParamsConfig
{
   public:

   class SPHConfig
   {
      public:
      /**
       * Define the device on which the code should run.
       */
      using DeviceType = TNL::Devices::Cuda;

      /**
       * Definition of basics data types for fluid variables and indices.
       */
      using GlobalIndexType = int;
      using LocalIndexType = int;
      using CellIndexType = int;
      using RealType = float;

      /**
       * Define the space dimension of the problem.
       */
      static constexpr int spaceDimension = placeholderDimension;
   };

   /**
    * Define SPH parameters connected to the resolution.
    * - h - smoothing length [m]
    * - dp - initial particle distance [m]
    */
   float dp = placeholderInitParticleDistancef;
   float h = placeholderSmoothingLengthf;

   /**
    * Define Basics SPH constants.
    * - mass - particle mass [kg]
    */
   float mass = placeholderMassf;

   /**
    * Define coefficient of diffusive term (DT), [-].
    */
   using DiffusiveTerm = TNL::ParticleSystem::SPH::MolteniDiffusiveTerm< SPHConfig >;
   float delta = 0.1f;

   /**
    * Define coefficient of artificial viscosity.
    */
   using ViscousTerm = TNL::ParticleSystem::SPH::ArtificialViscosity< SPHConfig >;
   float alpha = 0.02f;

   /**
    * Define equation of state and its constants.
    * - speedOfSound - numerical speed of sound used in the SPH calculations [m/s]
    * - coefB - coefficient of the Tait equation of state coefB = c^2 * rho0 / gamma
    * - rho0 - referential density of the fluid [kg/m^3]
    */
   using EOS = TNL::ParticleSystem::SPH::TaitWeaklyCompressibleEOS< SPHConfig >;
   float speedOfSound = placeholderSpeedOfSoundf;
   float coefB = placeholderCoefBf;
   float rho0 = placeholderDensityf;

   /**
    * Define initial timestep [s].
    */
   float dtInit = placeholderTimeStepf;

   /**
    * Define external forces [m^2 / s].
    */
   float gravity[ 2 ] { 0.f, -9.81f };

   /**
    * Define constant to prevent zero in denominator [-].
    */
   float eps = 0.001f;
};


} //SPH
} //namespace ParticleSystem
} //namespace TNL

