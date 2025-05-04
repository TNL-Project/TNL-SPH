#pragma once

#include <SPH/Models/DiffusiveTerms.h>
#include <SPH/Kernels.h>
#include <SPH/Models/SHTC/stress.h>
#include <SPH/Models/SHTC/IntegrationSchemes/ExplicitEulerScheme.h>

#include <SPH/SPHTraits.h>
#include <SPH/TimeStep.h>

//// FIXME: Not necessary
//#include "../../OpenBoundaryConfig.h"

namespace TNL {
namespace SPH {

/**
 * \brief Class used to store core parameters of SPH scheme.
 */
template< typename SPHDefs >
class SHTCConfig
{
public:

   using SPHConfig = typename SPHDefs::SPHConfig;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   static void
   configSetupModel( TNL::Config::ConfigDescription& config )
   {
      config.addDelimiter( "WCSPH-SHTC model parameters" );
      config.addEntry< double >( "dp", "Initial particle distance.", 0 );
      config.addEntry< double >( "h", "SPH method smoothing lentgh.", 0 );
      config.addEntry< double >( "mass", "Mass of particle, constant for all particles.", 0 );
      config.addEntry< double >( "delta", "Coefficient of artificial delta-WCSPH diffusive term.", 0 );
      config.addEntry< double >( "dynamicViscosity", "Dynamic viscosity coefficient.", 0 );
      config.addEntry< double >( "speedOfSound_bulk", "Numerical speed of sound.", 0 );
      config.addEntry< double >( "speedOfSound_shear", "Numerical speed of sound.", 0 );
      config.addEntry< double >( "rho0", "Referential density of the medium.", 0 );
      config.addEntry< double >( "tau", "Relaxation time.", 0 );
      config.addEntry< double >( "antiClumpFactor", "Anti-clumping factor..", 1e-3 );
      config.addEntry< double >( "initial-time-step", "Initial time step.", 0 );
      config.addEntry< double >( "CFL", "CFL number.", 0 );
      config.addEntry< double >( "minimal-time-step", "Minimal allowed time step.", 0 );
      config.addEntry< double >( "external-force-x", "External bulk forces.", 0 );
      config.addEntry< double >( "external-force-y", "External bulk forces.", 0 );
      config.addEntry< double >( "external-force-z", "External bulk forces.", 0 );
      config.addEntry< double >( "eps", "Coefficient to prevent denominator from zero.", 0 );
   }

   void
   init( TNL::Config::ParameterContainer& parameters )
   {
      h = parameters.getParameter< double >( "h" );
      dp = parameters.getParameter< double >( "dp" );
      mass = parameters.getParameter< double >( "mass" );
      delta = parameters.getParameter< double >( "delta" );
      dynamicViscosity = parameters.getParameter< double >( "dynamicViscosity" );
      speedOfSound = parameters.getParameter< double >( "speedOfSound_bulk" ); //FIXME: keep default name due to interaces
      speedOfSound_bulk = parameters.getParameter< double >( "speedOfSound_bulk" );
      speedOfSound_shear = parameters.getParameter< double >( "speedOfSound_shear" );
      rho0 = parameters.getParameter< double >( "rho0" );
      tau = parameters.getParameter< double >( "tau" );
      antiClumpFactor = parameters.getParameter< double >( "antiClumpFactor" );
      dtInit = parameters.getParameter< double >( "initial-time-step" );
      cfl = parameters.getParameter< double >( "CFL" );
      dtMin = parameters.getParameter< double >( "minimal-time-step" );
      eps = parameters.getParameter< double >( "eps" );
      gravity = parameters.getXyz< VectorType >( "external-force" );

      dtMin = 0.05f * dtInit;
   }

   //dp - initial particle distance [m]
   RealType dp = 0.f;
   //h - smoothing length [m]
   RealType h = 0.f;
   //mass - particle mass [kg]
   RealType mass = 0.f;
   //SPH weight function (kernel).
   using KernelFunction = typename SPHDefs::KernelFunction;

   //Diffusive term type
   using DiffusiveTerm = typename SPHDefs::DiffusiveTerm;
   // Define coefficient of diffusive term (DT), [-].
   RealType delta = 0.1f;

   //TODO: Add viscous term to RSPH scheme
   //Define viscous term and its coefficients.
   //using ViscousTerm = typename SPHDefs::ViscousTerm;
   //dynamicViscosity - value of dynamic viscosity [Pa/s];
   RealType dynamicViscosity = 1e-3f;

   // Define equation of state and its constants.
   using Stress = typename SPHDefs::Stress;
   //speedOfSound - numerical speed of sound used in the SPH calculations [m/s]
   RealType speedOfSound = 0.f; //FIXME - keep default name due to interfaces
   RealType speedOfSound_bulk = 0.f;
   RealType speedOfSound_shear = 0.f;
   //rho0 - referential density of the fluid [kg/m^3]
   RealType rho0 = 0.f;
   //tauFact - relaxation factor tau [?]
   RealType tau = 0.f;

   //TODO: Add boundary conditions types to RSPH scheme
   //Define type of boundary conditions.
   //using BCType = typename SPHDefs::BCType;

   // antiClumpFactor - anti clumping factor [-]
   RealType antiClumpFactor = 1e-3f;

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
   logger.writeHeader( "TNL::SPH::SHTC model parameters" );
   logger.writeParameter( "Resolution parameters", "" );
   logger.writeParameter( "Initial particle distance (dp):", modelParams.dp, 1 );
   logger.writeParameter( "Smoothing length (h):", modelParams.h, 1 );
   logger.writeParameter( "Spatial resolution (h/dp):", modelParams.h / modelParams.dp, 1 );
   logger.writeParameter( "Particle mass (mass):", modelParams.mass, 1 );
   logger.writeParameter( "Model parameters", "" );
   if constexpr ( std::is_same_v< typename ModelParams::DiffusiveTerm, DiffusiveTerms::MolteniDiffusiveTerm< typename ModelParams::SPHConfig> > ){
      logger.writeParameter( "Diffusive term:", "TNL::SPH::MolteniDiffusiveTerm", 1 );
      logger.writeParameter( "Diffusive term coefficient (delta):", modelParams.delta, 1 );
   }
   if constexpr ( std::is_same_v< typename ModelParams::Stress, Stress::FluidStress< typename ModelParams::SPHConfig> > ){
      logger.writeParameter( "Stress model:", "TNL::SPH::FluidStress", 1 );
   }
   logger.writeParameter( "Build speed of sound (speedOfSound_bulk):", modelParams.speedOfSound_bulk, 1 );
   logger.writeParameter( "Shear speed of sound (speedOfSound_shear):", modelParams.speedOfSound_shear, 1 );
   logger.writeParameter( "Relaxation time (tau):", modelParams.tau, 1 );
   logger.writeParameter( "Referentail density (rho0):", modelParams.rho0, 1 );
   logger.writeParameter( "Time integration", "" );
   if constexpr ( std::is_same_v< typename ModelParams::IntegrationScheme, IntegrationSchemes::ExplicitEulerScheme< typename ModelParams::SPHConfig> > )
      logger.writeParameter( "Integration scheme:", "TNL::SPH::SHTC::ExplicitEulerScheme", 1 );
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

