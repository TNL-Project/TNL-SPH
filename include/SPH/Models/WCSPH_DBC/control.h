#pragma once

#include <SPH/Models/EquationOfState.h>
#include <SPH/Models/DiffusiveTerms.h>
#include <SPH/Models/VisousTerms.h>
#include <SPH/Kernels.h>

#include <SPH/Models/WCSPH_DBC/BoundaryConditionsTypes.h>

#include <SPH/SPHTraits.h>
#include <SPH/TimeStep.h>
#include <limits>

namespace TNL {
namespace ParticleSystem {
namespace SPH {



template< typename SPHConfig >
bool SPHModelInit( TNL::Config::ConfigDescription config )
{
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   config.addEntry< RealType >( "dp", "Initial particle distance.", 0 );
   config.addEntry< RealType >( "h", "SPH method smoothing lentgh.", 0 );
   config.addEntry< RealType >( "mass", "Mass of particle, constant for all particles.", 0 );
   config.addEntry< RealType >( "delta", "Coefficient of artificial delta-WCSPH diffusive term.", 0 );
   config.addEntry< RealType >( "alpha", "Coefficient of artificial viscous term.", 0 );
   config.addEntry< RealType >( "dynamicViscosity", "Dynamic viscosity coefficient.", 0 );
   config.addEntry< RealType >( "speedOfSound", "Numerical speed of sound.", 0 );
   config.addEntry< RealType >( "rho0", "Referential density of the medium.", 0 );
   //float coefB = 168070.0f;
   config.addEntry< RealType >( "dtInit", "Initial time step.", 0 );
   config.addEntry< RealType >( "CFL", "CFL number.", 0 );
   config.addEntry< RealType >( "dtMin", "Minimal allowed time step.", 0 );
   config.addEntry< VectorType >( "gravity", "External bulk forces.", 0 );
   config.addEntry< RealType >( "eps", "Coefficient to prevent denominator from zero.", 0 );
}

/**
 * PARAMETERS OF SPH SIMULATION AND SCHEMES (necessary)
 *
 * This class is used to store core parameters for simulation control and
 * simulation initialization. This includes variables such as paths for loading
 * and saving files or the length of the simulation and the frequency of saving outputs.
 *
 */
template< typename Devices >
class WCSPH_DBCConfig
{
   public:
   using SPHConfig = SPHConfig< Device >;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;

   /**
    * Define SPH parameters connected to the resolution.
    * - h - smoothing length [m]
    * - dp - initial particle distance [m]
    */
   float dp = 0.002f;
   float h = 0.0028284f;

   /**
    * Define SPH weight function (kernel).
    * - Use "WendlandKernel" for 4th order Wendland kernel.
    */
   using KernelFunction = TNL::ParticleSystem::SPH::WendlandKernel< SPHConfig >;

   /**
    * Define Basics SPH constants.
    * - mass - particle mass [kg]
    */
   float mass = 0.004f;

   /**
    * Define coefficient of diffusive term (DT), [-].
    */
   using DiffusiveTerm = TNL::ParticleSystem::SPH::MolteniDiffusiveTerm< SPHConfig >;
   float delta = 0.1f;

   /**
    * Define viscous term and its coefficients.
    * - Use "ArtificialViscosity" and its parameter alpha.
    * - Use "PhysicalViscosity" defined by dynamic viscosity coeffition.
    */
   using ViscousTerm = TNL::ParticleSystem::SPH::ArtificialViscosity< SPHConfig >;
   float alpha = 0.02f;

   //using ViscousTerm = TNL::ParticleSystem::SPH::PhysicalViscosity< SPHConfig >;
   //float dynamicViscosity = 1e-3f;

   /**
    * Define equation of state and its constants.
    * - speedOfSound - numerical speed of sound used in the SPH calculations [m/s]
    * - coefB - coefficient of the Tait equation of state coefB = c^2 * rho0 / gamma
    * - rho0 - referential density of the fluid [kg/m^3]
    */
   using EOS = TNL::ParticleSystem::SPH::TaitWeaklyCompressibleEOS< SPHConfig >;
   float speedOfSound = 34.3f;
   float coefB = 168070.0f;
   float rho0 = 1000.0f;

   /**
    * Define type of boundary conditions.
    * - DBC - dynamic boundary conditions
    * - MDBC - modified dynamic boundary conditions {requires ghost nodes for boundary particles}
    */
   using BCType = TNL::ParticleSystem::SPH::WCSPH_BCTypes::DBC;
   //using BCType = TNL::ParticleSystem::SPH::WCSPH_BCTypes::MDBC;

   /**
    * Define initial timestep [s].
    * - Use "ConstantTimeStep" with dtInit representing the step [s].
    * - Use "VariableTimeStep" with CFL number [-], initial tiem steop [s] and minimum timestep [s].
    * - Use "VariableTimeStepWithReduction" with CFL number [-], initial tiem steop [s] and minimum timestep [s].
    */
   using TimeStepping = TNL::ParticleSystem::SPH::ConstantTimeStep< SPHConfig >;
   float dtInit = 1e-05f;

   //using TimeStepping = TNL::ParticleSystem::SPH::VariableTimeStep< SPHConfig >;
   //float CFL = 0.3f;
   //float dtInit = 0.25f * h / speedOfSound;
   //float dtMin = 0.05f * h / speedOfSound;

   /**
    * Define external forces [m^2 / s].
    */
   float gravity[ 2 ] { 0.f, 0.f };

   /**
    * Define constant to prevent zero in denominator [-].
    */
   float eps = 0.001f;
};


} //namespace SPH
} //namespace ParticleSystem
} //namespace TNL

