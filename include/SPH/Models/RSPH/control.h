#pragma once

#include <SPH/Models/EquationOfState.h>
#include <SPH/Models/RiemannSolvers.h>
#include <SPH/Kernels.h>
#include <SPH/Models/RSPH/IntegrationSchemes/VerletScheme.h>


#include <SPH/SPHTraits.h>
#include <SPH/TimeStep.h>
#include <complex>
#include <limits>

namespace TNL {
namespace SPH {

template< typename SPHConfig >
void
configSetupModel( TNL::Config::ConfigDescription& config )
{
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   config.addDelimiter( "WCSPH-DBC model parameters" );
   config.addEntry< float >( "dp", "Initial particle distance.", 0 );
   config.addEntry< float >( "h", "SPH method smoothing lentgh.", 0 );
   config.addEntry< float >( "mass", "Mass of particle, constant for all particles.", 0 );
   config.addEntry< float >( "limiterEta", "Coefficient of modified linearized Roe scheme.", 0 );
   config.addEntry< float >( "dynamicViscosity", "Dynamic viscosity coefficient.", 0 );
   config.addEntry< float >( "speedOfSound", "Numerical speed of sound.", 0 );
   config.addEntry< float >( "rho0", "Referential density of the medium.", 0 );
   config.addEntry< RealType >( "initial-time-step", "Initial time step.", 0 );
   config.addEntry< RealType >( "CFL", "CFL number.", 0 );
   config.addEntry< RealType >( "minimal-time-step", "Minimal allowed time step.", 0 );
   config.addEntry< RealType >( "external-force-x", "External bulk forces.", 0 );
   config.addEntry< RealType >( "external-force-y", "External bulk forces.", 0 );
   config.addEntry< RealType >( "external-force-z", "External bulk forces.", 0 );
   config.addEntry< RealType >( "eps", "Coefficient to prevent denominator from zero.", 0 );
}

/**
 * \brief Class used to store core parameters of SPH scheme.
 */
template< typename SPHDefs >
class RSPHConfig
{
public:

   using SPHConfig = typename SPHDefs::SPHConfig;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   void
   init( TNL::Config::ParameterContainer& parameters )
   {
      h = parameters.getParameter< RealType >( "h" );
      dp = parameters.getParameter< RealType >( "dp" );
      mass = parameters.getParameter< RealType >( "mass" );
      limiterEta = parameters.getParameter< RealType >( "limiterEta" );
      dynamicViscosity = parameters.getParameter< RealType >( "dynamicViscosity" );
      speedOfSound = parameters.getParameter< RealType >( "speedOfSound" );
      rho0 = parameters.getParameter< RealType >( "rho0" );
      dtInit = parameters.getParameter< RealType >( "initial-time-step" );
      cfl = parameters.getParameter< RealType >( "CFL" );
      dtMin = parameters.getParameter< RealType >( "minimal-time-step" );
      eps = parameters.getParameter< RealType >( "eps" );
      gravity = parameters.getXyz< VectorType >( "external-force" );

      coefB = speedOfSound * speedOfSound * rho0 / 7.f;
      dtMin = 0.05f * h / speedOfSound;
   }

   //dp - initial particle distance [m]
   RealType dp = 0.f;
   //h - smoothing length [m]
   RealType h = 0.f;
   //mass - particle mass [kg]
   RealType mass = 0.f;
   //SPH weight function (kernel).
   using KernelFunction = typename SPHDefs::KernelFunction;

   //Define Riemann solver and its coefficients.
   using RiemannSolver = typename SPHDefs::RiemannSolver;
   //limiterEta - constant of in-build limiter of modified Roe scheme [-]
   RealType limiterEta = 3.0f;

   //TODO: Add viscous term to RSPH scheme
   //Define viscous term and its coefficients.
   //using ViscousTerm = typename SPHDefs::ViscousTerm;
   //dynamicViscosity - value of dynamic viscosity [Pa/s];
   RealType dynamicViscosity = 1e-3f;

   // Define equation of state and its constants.
   using EOS = typename SPHDefs::EOS;
   //speedOfSound - numerical speed of sound used in the SPH calculations [m/s]
   RealType speedOfSound = 0.f;
   //coefB - coefficient of the Tait equation of state coefB = c^2 * rho0 / gamma
   RealType coefB = 0.f;
   //rho0 - referential density of the fluid [kg/m^3]
   RealType rho0 = 0.f;

   //TODO: Add boundary conditions types to RSPH scheme
   //Define type of boundary conditions.
   //using BCType = typename SPHDefs::BCType;

   //Type of integration scheme
   using IntegrationScheme = typename SPHDefs::IntegrationScheme;
   //Type of time stepping scheme
   using TimeStepping = typename SPHDefs::TimeStepping;
   //dtInit - initial time step [ s ]
   RealType dtInit = 0.f;
   //cfl - CFL number [-];
   RealType cfl = 0.f;
   //dtMin - minimal allowed time step [s];
   RealType dtMin = 0.f;

   //gravity - external forces [m^2 / s].
   VectorType gravity = 0.f;

   //eps - constant to prevent zero in denominator [-].
   RealType eps = 0.001f;
};

template< typename ModelParams >
void writePrologModel( TNL::Logger& logger, ModelParams& modelParams )
{
   logger.writeHeader( "TNL::SPH::RSPH model parameters" );
   logger.writeParameter( "Resolution parameters", "" );
   logger.writeParameter( "Initial particle distance (dp):", modelParams.dp, 1 );
   logger.writeParameter( "Smoothing length (h):", modelParams.h, 1 );
   logger.writeParameter( "Spatial resolution (dp/h):", modelParams.dp / modelParams.h, 1 );
   logger.writeParameter( "Particle mass (mass):", modelParams.mass, 1 );
   logger.writeParameter( "Model parameters", "" );
   if constexpr ( std::is_same_v< typename ModelParams::RiemannSolver, RiemannSolvers::RoeLinearized< typename ModelParams::SPHConfig> > ){
      logger.writeParameter( "Riemann solver:", "TNL::SPH::RoeLinearized", 1 );
      logger.writeParameter( "Roe scheme limiter constant:", modelParams.limiterEta, 1 );
   }
   if constexpr ( std::is_same_v< typename ModelParams::EOS, EquationsOfState::TaitWeaklyCompressibleEOS< typename ModelParams::SPHConfig> > ){
      logger.writeParameter( "Equation of state:", "TNL::SPH::TaitWeaklyCompressibleEOS", 1 );
      logger.writeParameter( "Coefficient of EOS (coefB): ", modelParams.coefB );
   }
   if constexpr ( std::is_same_v< typename ModelParams::EOS, EquationsOfState::TaitLinearizedWeaklyCompressibleEOS< typename ModelParams::SPHConfig> > )
      logger.writeParameter( "Equation of state:", "TNL::SPH::LinearizedTaitWeaklyCompressibleEOS", 1 );
   logger.writeParameter( "Speed of sound (speedOfSound):", modelParams.speedOfSound, 1 );
   logger.writeParameter( "Referentail density (rho0):", modelParams.rho0, 1 );
   logger.writeParameter( "Time integration", "" );
   if constexpr ( std::is_same_v< typename ModelParams::IntegrationScheme, IntegrationSchemes::VerletScheme< typename ModelParams::SPHConfig> > )
      logger.writeParameter( "Integration scheme:", "TNL::SPH::RSPH::VerletScheme", 1 );
   if constexpr ( std::is_same_v< typename ModelParams::TimeStepping, ConstantTimeStep< typename ModelParams::SPHConfig> > ){
      logger.writeParameter( "Time stepping:", "TNL::SPH::ConstantTimeStep", 1 );
      logger.writeParameter( "Initial time step (dtInit):", modelParams.dtInit, 1 );
   }
   if constexpr ( std::is_same_v< typename ModelParams::TimeStepping, VariableTimeStep< typename ModelParams::SPHConfig> > ){
      logger.writeParameter( "Time stepping:", "TNL::SPH::VariableTimeStep", 1 );
      logger.writeParameter( "Initial time step (dtInit):", modelParams.dtInit, 1 );
      logger.writeParameter( "Minimal time step (dtMin):", modelParams.dtMin, 1 );
      logger.writeParameter( "CFL number (CFL):", modelParams.cfl, 1 );
   }
   logger.writeParameter( "External bulk force:", modelParams.gravity );
}

} //namespace SPH
} //namespace TNL

