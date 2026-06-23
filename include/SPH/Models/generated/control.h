// control.h  –  AUTO-GENERATED
#pragma once

#include <SPH/SPHTraits.h>
#include <SPH/TimeStep.h>
#include <SPH/Kernels.h>
#include <SPH/Models/EquationOfState.h>
#include <SPH/Models/DiffusiveTerms.h>
#include <SPH/Models/VisousTerms.h>
#include <SPH/Models/WCSPH_DBC/BoundaryConditionsTypes.h>
#include <SPH/Models/WCSPH_DBC/IntegrationSchemes/VerletScheme.h>
#include <SPH/Models/WCSPH_DBC/OpenBoundaryConfig.h>

namespace TNL {
namespace SPH {

// Parameters for WCSPH_DBC
template< typename SPHDefs >
class WCSPH_DBCConfig
{
public:
   using SPHConfig         = typename SPHDefs::SPHConfig;
   using SPHTraitsType     = SPHFluidTraits< SPHConfig >;
   using RealType          = typename SPHTraitsType::RealType;
   using VectorType        = typename SPHTraitsType::VectorType;
   using KernelFunction    = typename SPHDefs::KernelFunction;
   using DiffusiveTerm     = typename SPHDefs::DiffusiveTerm;
   using ViscousTerm       = typename SPHDefs::ViscousTerm;
   using EOS               = typename SPHDefs::EOS;
   using BCType            = typename SPHDefs::BCType;
   using IntegrationScheme = typename SPHDefs::IntegrationScheme;
   using TimeStepping      = typename SPHDefs::TimeStepping;

   static void configSetupModel( TNL::Config::ConfigDescription& config )
   {
      config.addDelimiter( "WCSPH_DBC model parameters" );
      config.addEntry< double >( "h", "SPH method smoothing length.", 0 );
      config.addEntry< double >( "dp", "Initial particle distance.", 0 );
      config.addEntry< double >( "mass", "Mass of particle, constant for all particles.", 0 );
      config.addEntry< double >( "delta", "Coefficient of artificial delta-WCSPH diffusive term.", 0 );
      config.addEntry< double >( "alpha", "Coefficient of artificial viscous term.", 0 );
      config.addEntry< double >( "dynamicViscosity", "Dynamic viscosity coefficient.", 0 );
      config.addEntry< double >( "speedOfSound", "Numerical speed of sound.", 0 );
      config.addEntry< double >( "rho0", "Referential density of the medium.", 0 );
      config.addEntry< double >( "dtInit", "Initial time step.", 0 );
      config.addEntry< double >( "cfl", "CFL number.", 0 );
      config.addEntry< double >( "dtMin", "Minimal allowed time step.", 0 );
      config.addEntry< double >( "eps", "Coefficient to prevent denominator from zero.", 0 );
      config.addEntry< double >( "gravity-x", "External bulk forces. x", 0 );
      config.addEntry< double >( "gravity-y", "External bulk forces. y", 0 );
      config.addEntry< double >( "gravity-z", "External bulk forces. z", 0 );
   }

   void init( TNL::Config::ParameterContainer& parameters )
   {
      h = parameters.getParameter< double >( "h" );
      dp = parameters.getParameter< double >( "dp" );
      mass = parameters.getParameter< double >( "mass" );
      delta = parameters.getParameter< double >( "delta" );
      alpha = parameters.getParameter< double >( "alpha" );
      dynamicViscosity = parameters.getParameter< double >( "dynamicViscosity" );
      speedOfSound = parameters.getParameter< double >( "speedOfSound" );
      rho0 = parameters.getParameter< double >( "rho0" );
      dtInit = parameters.getParameter< double >( "dtInit" );
      cfl = parameters.getParameter< double >( "cfl" );
      dtMin = parameters.getParameter< double >( "dtMin" );
      eps = parameters.getParameter< double >( "eps" );
      gravity = parameters.getXyz< VectorType >( "gravity" );
      coefB = speedOfSound * speedOfSound * rho0 / 7.f;
      dtMin = 0.05f * dtInit;
   }

   RealType   h = 0.f;
   RealType   dp = 0.f;
   RealType   mass = 0.f;
   RealType   delta = 0.f;
   RealType   alpha = 0.f;
   RealType   dynamicViscosity = 0.f;
   RealType   speedOfSound = 0.f;
   RealType   rho0 = 0.f;
   RealType   dtInit = 0.f;
   RealType   cfl = 0.f;
   RealType   dtMin = 0.f;
   RealType   eps = 0.f;
   RealType   coefB = 0.f;
   RealType   p0 = 0.f;
   VectorType gravity = 0.f;

};

template< typename ModelParams >
void
writePrologModel( TNL::Logger& logger, ModelParams& modelParams )
{
   logger.writeHeader( "TNL::SPH::WCSPH_DBC model parameters" );
   logger.writeParameter( "Resolution parameters", "" );
   logger.writeParameter( "Initial particle distance (dp):", modelParams.dp, 1 );
   logger.writeParameter( "Smoothing length (h):", modelParams.h, 1 );
   logger.writeParameter( "Particle mass (mass):", modelParams.mass, 1 );
   logger.writeParameter( "Spatial resolution (h/dp):", modelParams.h / modelParams.dp, 1 );
   logger.writeParameter( "Model parameters", "" );
   logger.writeParameter( "Diffusive term coefficient (delta):", modelParams.delta, 1 );
   logger.writeParameter( "Artificial viscosity coefficient (alpha):", modelParams.alpha, 1 );
   logger.writeParameter( "Dynamic viscosity (dynamicViscosity):", modelParams.dynamicViscosity, 1 );
   logger.writeParameter( "Speed of sound (speedOfSound):", modelParams.speedOfSound, 1 );
   logger.writeParameter( "Referential density (rho0):", modelParams.rho0, 1 );
   logger.writeParameter( "Initial time step (dtInit):", modelParams.dtInit, 1 );
   logger.writeParameter( "CFL number (CFL):", modelParams.cfl, 1 );
   logger.writeParameter( "Minimal time step (dtMin):", modelParams.dtMin, 1 );
   logger.writeParameter( "Epsilon (eps):", modelParams.eps, 1 );
   logger.writeParameter( "External bulk force:", modelParams.gravity );
   logger.writeParameter( "External bulk force:", modelParams.gravity );
}

} // SPH
} // TNL

// end control.h