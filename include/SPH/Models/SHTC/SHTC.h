#pragma once

#include "../../SPHTraits.h"
#include "Variables.h"

/**
 * Modules used as default.
 **/
#include "control.h"

// placeholder for empty open boundary conditions
#include "../emptyOpenBCPlaceholder.h"
#include "../SPHMultisetSolverTemplate.h"


namespace TNL {
namespace SPH {

template< typename Particles, typename ModelConfig >
class SHTC : public SPHMultisetSolverTemplate< Particles, ModelConfig >
{
public:

   using Model = SHTC< Particles, ModelConfig >;
   using ModelParams = SHTCConfig< ModelConfig >;
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
   using MatrixType = typename SPHTraitsType::MatrixType;

   using FluidVariables = SHTCVariables< SPHConfig >;
   using BoundaryVariables = SHTCBoundaryVariables< ModelConfig >;
   //using OpenBoundaryVariables = Empty< SPHConfig >;
   using OpenBoundaryVariables = SHTCVariables< SPHConfig >; //FIXME: Empty doesn't work since open bc are no constexpr anymore
   using IntegrationSchemeType = typename ModelConfig::IntegrationScheme;
   using IntegrationSchemeVariables = typename IntegrationSchemeType::IntegrationSchemeVariablesType;
   using KernelFunction = typename ModelConfig::KernelFunction;
   using DiffusiveTerm = typename ModelConfig::DiffusiveTerm;
   using Stress = typename ModelConfig::Stress;

   using OpenBoundaryConfig = NoOpenBC< SPHConfig >;
   using OpenBoundaryModel = Empty< SPHConfig >;

   /**
    * Constructor.
    */
   SHTC() = default;

   /**
    * Print model identifier.
    */
   static std::string
   writeModelType()
   {
      return "TNL::SPH::SHTC";
   }

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

   template< typename FluidPointer, typename TimeStepping >
   void
   relaxDistortion( FluidPointer& fluid, TimeStepping& timeStepping, ModelParams& modelParams );

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

   // Empty functions to match the solver interface

   //template< typename OpenBoundaryPointer, typename BoudaryPointer >
   //void
   //updateSolidBoundaryOpenBoundary( BoudaryPointer& boundary,
   //                                 OpenBoundaryPointer& openBoundaryPointer,
   //                                 ModelParams& modelParams )
   //{};

   //template< typename FluidPointer, typename OpenBoudaryPointer >
   //void
   //interactionWithOpenBoundary( FluidPointer& fluid, OpenBoudaryPointer& openBoundary, ModelParams& modelParams )
   //{}


};

} // SPH
} // TNL

#include "SHTC.hpp"

