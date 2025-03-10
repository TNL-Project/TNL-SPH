#pragma once

#include "../../SPHTraits.h"
#include "Variables.h"

/**
 * Modules used as default.
 **/
#include "control.h"

namespace TNL {
namespace SPH {

template< typename Particles, typename ModelConfig >
class RSPH
{
public:

   using Model = RSPH< Particles, ModelConfig >;
   using ModelParams = RSPHConfig< ModelConfig >;
   using ParticlesType = Particles;
   using ModelConfigType = ModelConfig;
   using SPHConfig = typename ModelConfig::SPHConfig;
   using DeviceType = typename SPHConfig::DeviceType;

   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using LocalIndexType = typename SPHTraitsType::LocalIndexType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using IndexVectorType = typename SPHTraitsType::IndexVectorType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   using FluidVariables = FluidVariables< SPHConfig >;
   using BoundaryVariables = BoundaryVariables< ModelConfig >;
   using OpenBoundaryVariables = Empty< SPHConfig >;
   using IntegrationSchemeType = typename ModelConfig::IntegrationScheme;
   using IntegrationSchemeVariables = typename IntegrationSchemeType::IntegrationSchemeVariablesType;
   using KernelFunction = typename ModelConfig::KernelFunction;
   using RiemannSolver = typename ModelConfig::RiemannSolver;
   using EOS = typename ModelConfig::EOS;

   using OpenBoundaryConfig = Empty< SPHConfig >;
   using OpenBoundaryModel = Empty< SPHConfig >;

   /**
    * Constructor.
    */
   RSPH() = default;

   /**
    * Print model identifier.
    */
   static std::string
   writeModelType()
   {
      return "TNL::SPH::RSPH";
   }

   /**
    * Compute pressure from density.
    */
   template< typename EquationOfState = EquationsOfState::TaitWeaklyCompressibleEOS< SPHConfig >,
             typename PhysicalObjectPointer >
   void
   computePressureFromDensity( PhysicalObjectPointer& variables, ModelParams& modelParams );

   /**
    * Function to realize fluid-fluid and fluid-boundary interaction.
    */
   template< typename FluidPointer, typename BoudaryPointer >
   void
   interaction( FluidPointer& fluid, BoudaryPointer& boundary, ModelParams& modelParams );

   /**
    * Function to realize boundary conditions for solid wall.
    */
   template< typename FluidPointer, typename BoudaryPointer >
   void
   updateSolidBoundary( FluidPointer& fluid, BoudaryPointer& boundary, ModelParams& modelParams );

   /**
    * Perform preinteraction procedures.
    */
   template< typename FluidPointer, typename BoundaryPointer >
   void
   initializeInteraction( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams ) {}

   /**
    * Perform postinteraction procedures.
    */
   template< typename FluidPointer, typename BoundaryPointer >
   void
   finalizeInteraction( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams ) {}

   /**
    * Perform postinteraction procedures.
    */
   template< typename FluidPointer, typename BoundaryPointer >
   void
   finalizeBoundaryInteraction( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams ) {}

};

} // SPH
} // TNL

#include "Interactions.hpp"

