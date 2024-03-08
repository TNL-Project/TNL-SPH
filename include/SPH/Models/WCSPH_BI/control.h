#pragma once

#include <SPH/Models/EquationOfState.h>
#include <SPH/Models/DiffusiveTerms.h>
#include <SPH/Models/VisousTerms.h>
#include <SPH/Kernels.h>
#include <SPH/Models/WCSPH_BI/IntegrationSchemes/VerletScheme.h>

#include <SPH/Models/WCSPH_BI/BoundaryConditionsTypes.h>

#include <SPH/SPHTraits.h>
#include <SPH/TimeStep.h>
#include <complex>
#include <limits>

#include "OpenBoundaryConfig.h"

namespace TNL {
namespace SPH {

template< typename SPHConfig >
void
configSetupModel( TNL::Config::ConfigDescription& config )
{
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   config.addDelimiter( "WCSPH-BI model parameters" );

   config.addEntry< float >( "dp", "Initial particle distance.", 0 );
   config.addEntry< float >( "h", "SPH method smoothing lentgh.", 0 );
   config.addEntry< float >( "boundaryElementSize", "Size of bounadry inegrals element.", 0 );
   config.addEntry< float >( "mass", "Mass of particle, constant for all particles.", 0 );
   config.addEntry< float >( "delta", "Coefficient of artificial delta-WCSPH diffusive term.", 0 );
   config.addEntry< float >( "alpha", "Coefficient of artificial viscous term.", 0 );
   config.addEntry< float >( "dynamicViscosity", "Dynamic viscosity coefficient.", 0 );
   config.addEntry< float >( "speedOfSound", "Numerical speed of sound.", 0 );
   config.addEntry< float >( "rho0", "Referential density of the medium.", 0 );
   config.addEntry< RealType >( "dtInit", "Initial time step.", 0 );
   config.addEntry< RealType >( "CFL", "CFL number.", 0 );
   config.addEntry< RealType >( "dtMin", "Minimal allowed time step.", 0 );
   config.addEntry< RealType >( "external-force-x", "External bulk forces.", 0 );
   config.addEntry< RealType >( "external-force-y", "External bulk forces.", 0 );
   config.addEntry< RealType >( "external-force-z", "External bulk forces.", 0 );
   config.addEntry< RealType >( "eps", "Coefficient to prevent denominator from zero.", 0 );

   for( int i = 0; i < SPHConfig::numberOfBoundaryBuffers; i++ )
   {
      std::string prefix = "buffer-" + std::to_string( i + 1 ) + "-";
      configSetupOpenBoundaryModelPatch< SPHConfig >( config, prefix );
   }
   for( int i = 0; i < SPHConfig::numberOfPeriodicBuffers; i++ )
   {
      std::string prefix = "buffer-" + std::to_string( i + 1 ) + "-";
      configSetupOpenBoundaryModelPatch< SPHConfig >( config, prefix );
   }
}

/**
 * \brief Class used to store core parameters of SPH scheme.
 */
template< typename SPHDefs >
class WCSPH_BIConfig
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
      searchRadius = parameters.getParameter< RealType >( "searchRadius" );
      mass = parameters.getParameter< RealType >( "mass" );
      boundaryElementSize = parameters.getParameter< RealType >( "boundaryElementSize" );
      delta = parameters.getParameter< RealType >( "delta" );
      alpha = parameters.getParameter< RealType >( "alpha" );
      dynamicViscosity = parameters.getParameter< RealType >( "dynamicViscosity" );
      speedOfSound = parameters.getParameter< RealType >( "speedOfSound" );
      rho0 = parameters.getParameter< RealType >( "rho0" );
      dtInit = parameters.getParameter< RealType >( "dtInit" );
      CFL = parameters.getParameter< RealType >( "CFL" );
      dtMin = parameters.getParameter< RealType >( "dtMin" );
      eps = parameters.getParameter< RealType >( "eps" );
      gravity = parameters.getXyz< VectorType >( "external-force" );

      coefB = speedOfSound * speedOfSound * rho0 / 7.f;
   }

   //dp - initial particle distance [m]
   RealType dp = 0.f;
   //h - smoothing length [m]
   RealType h = 0.f;
   //searchRadius - radius of kernel support [m]
   RealType searchRadius = 0.f;
   //mass - particle mass [kg]
   RealType mass = 0.f;
   //boundaryElementSize - size of boundary element [m]
   RealType boundaryElementSize = 0.f;
   //SPH weight function (kernel).
   using KernelFunction = typename SPHDefs::KernelFunction;

   //Diffusive term type
   using DiffusiveTerm = typename SPHDefs::DiffusiveTerm;
   // Define coefficient of diffusive term (DT), [-].
   RealType delta = 0.1f;

   //Define viscous term and its coefficients.
   using ViscousTerm = typename SPHDefs::ViscousTerm;
   //alpha - parameter of artificial viscosity [-]
   RealType alpha = 0.02f;
   //dynamicViscosity - value of dynamic viscosity [Pa/s];
   RealType dynamicViscosity = 0.f;

   //Viscosity model for boundary interaction
   using BoundaryViscousTerm = typename SPHDefs::BoundaryViscousTerm;

   // Define equation of state and its constants.
   using EOS = typename SPHDefs::EOS;
   //speedOfSound - numerical speed of sound used in the SPH calculations [m/s]
   RealType speedOfSound = 0.f;
   //coefB - coefficient of the Tait equation of state coefB = c^2 * rho0 / gamma
   RealType coefB = 0.f;
   //rho0 - referential density of the fluid [kg/m^3]
   RealType rho0 = 0.f;

   //Define type of boundary conditions.
   using BCType = typename SPHDefs::BCType;

   //Type of integration scheme
   using IntegrationScheme = typename SPHDefs::IntegrationScheme;
   //Type of time stepping scheme
   using TimeStepping = typename SPHDefs::TimeStepping;
   //dtInit - initial time step [ s ]
   RealType dtInit = 0.f;
   //CFL - cfl number [-];
   RealType CFL = 0.f;
   //dtMin - minimal allowed time step [s];
   RealType dtMin = 0.05f * h / speedOfSound;

   //gravity - external forces [m^2 / s].
   VectorType gravity = 0.f;

   //eps - constant to prevent zero in denominator [-].
   RealType eps = 0.001f;


};

template< typename ModelParams >
void writePrologModel( TNL::Logger& logger, ModelParams& modelParams )
{
   logger.writeHeader( "TNL::SPH::WCSPH_BI (delta-WCSPH) model parameters" );
   logger.writeParameter( "Resolution parameters", "" );
   logger.writeParameter( "Initial particle distance (dp):", modelParams.dp, 1 );
   logger.writeParameter( "Smoothing length (h):", modelParams.h, 1 );
   logger.writeParameter( "Spatial resolution (dp/h):", modelParams.dp / modelParams.h, 1 );
   logger.writeParameter( "Search radius (searchRadius):", modelParams.searchRadius, 1 );
   logger.writeParameter( "Particle mass (mass):", modelParams.mass, 1 );
   logger.writeParameter( "Size of boundary elements (boundaryElementSize):", modelParams.boundaryElementSize, 1 );
   logger.writeParameter( "Model parameters", "" );
   if constexpr ( std::is_same_v< typename ModelParams::DiffusiveTerm, DiffusiveTerms::MolteniDiffusiveTerm< typename ModelParams::SPHConfig> > ){
      logger.writeParameter( "Diffusive term:", "TNL::SPH::MolteniDiffusiveTerm", 1 );
      logger.writeParameter( "Diffusive term coefficient (delta):", modelParams.delta, 1 );
   }
   if constexpr ( std::is_same_v< typename ModelParams::ViscousTerm, ViscousTerms::ArtificialViscosity< typename ModelParams::SPHConfig> > ){
      logger.writeParameter( "Viscous term:", "TNL::SPH::ArtificialViscosity", 1 );
      logger.writeParameter( "Artificial vicosity coefficient (alpha):", modelParams.alpha, 1 );
   }
   if constexpr ( std::is_same_v< typename ModelParams::ViscousTerm, ViscousTerms::PhysicalViscosity< typename ModelParams::SPHConfig> > ){
      logger.writeParameter( "Viscous term:", "TNL::SPH::PhysicalViscosity", 1 );
      logger.writeParameter( "Dynamic viscosity (dynamicViscosity):", modelParams.dynamicViscosity, 1 );
   }
   if constexpr ( std::is_same_v< typename ModelParams::EOS, EquationsOfState::TaitWeaklyCompressibleEOS< typename ModelParams::SPHConfig> > ){
      logger.writeParameter( "Equation of state:", "TNL::SPH::TaitWeaklyCompressibleEOS", 1 );
      logger.writeParameter( "Coefficient of EOS (coefB): ", modelParams.coefB );
   }
   if constexpr ( std::is_same_v< typename ModelParams::EOS, EquationsOfState::TaitLinearizedWeaklyCompressibleEOS< typename ModelParams::SPHConfig> > )
      logger.writeParameter( "Equation of state:", "TNL::SPH::LinearizedTaitWeaklyCompressibleEOS", 1 );
   logger.writeParameter( "Speed of sound (speedOfSound):", modelParams.speedOfSound, 1 );
   logger.writeParameter( "Referentail density (rho0):", modelParams.rho0, 1 );
   std::string boundaryConditionsTypes;
   if constexpr ( std::is_same_v< typename ModelParams::BCType, WCSPH_BCTypes::BI_numeric > )
      boundaryConditionsTypes = "TNL::SPH::WCSPH_BI::BI_numeric";
   logger.writeParameter( "Boundary condition type", boundaryConditionsTypes );
   logger.writeParameter( "Time integration", "" );
   if constexpr ( std::is_same_v< typename ModelParams::IntegrationScheme, IntegrationSchemes::VerletScheme< typename ModelParams::SPHConfig> > )
      logger.writeParameter( "Integration scheme:", "TNL::SPH::WCSPH_BI::VerletScheme", 1 );
   if constexpr ( std::is_same_v< typename ModelParams::TimeStepping, ConstantTimeStep< typename ModelParams::SPHConfig> > ){
      logger.writeParameter( "Time stepping:", "TNL::SPH::ConstantTimeStep", 1 );
      logger.writeParameter( "Initial time step (dtInit):", modelParams.dtInit, 1 );
   }
   if constexpr ( std::is_same_v< typename ModelParams::TimeStepping, VariableTimeStep< typename ModelParams::SPHConfig> > ){
      logger.writeParameter( "Time stepping:", "TNL::SPH::VariableTimeStep", 1 );
      logger.writeParameter( "Initial time step (dtInit):", modelParams.dtInit, 1 );
      logger.writeParameter( "Minimal time step (dtMin):", modelParams.dtMin, 1 );
      logger.writeParameter( "CFL number (CFL):", modelParams.dtMin, 1 );
   }
   logger.writeParameter( "External bulk force:", modelParams.gravity );
}

} //namespace SPH
} //namespace TNL

