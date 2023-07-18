#pragma once

#include "../../../../SPH/Models/EquationOfState.h"
#include "../../../../SPH/Models/DiffusiveTerms.h"
#include "../../../../SPH/Models/VisousTerms.h"

#include "../../../../SPH/SPHTraits.h"
#include <limits>

namespace TNL {
namespace ParticleSystem {
namespace SPH {
namespace SPHConfig {

/**
 * TYPES OF SPH SMULATION SYSTEM AND SCHEMES (necessary)
 *
 * This class is used to store parameters necessary for sph system,
 * i.e. data types for quantities and indices. It also defines dimension
 * and attributes of simulated system.
 *
 * It is necessary to enter TYPES for:
 * - GlobalIndexType
 * - LocalIndexType
 * - CellIndexType
 * - RealType
 */
template< typename Device >
class SPHConfig
{
   public:
   using DeviceType = Device;

   using GlobalIndexType = int;
   using LocalIndexType = int;
   using CellIndexType = int;
   using RealType = float;

   static constexpr int spaceDimension = placeholderDimension;
};


/**
 * PARAMETERS OF SPH SIMULATION AND SCHEMES (necessary)
 *
 * This class is used to store core parameters for simulation control and
 * simulation initialization. This includes variables such as paths for loading
 * and saving files or the length of the simulation and the frequency of saving outputs.
 *
 */
template< typename SPHConfig >
class SPHParamsConfig
{
   public:
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


} //namespace SPHConfig
} //SPH
} //namespace ParticleSystem
} //namespace TNL

