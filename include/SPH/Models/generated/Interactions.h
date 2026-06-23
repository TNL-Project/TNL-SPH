// Interactions.h  –  AUTO-GENERATED
#pragma once

#include "../../SPHTraits.h"
#include "../WCSPH_DBC/BoundaryConditionsTypes.h"
#include "../WCSPH_DBC/OpenBoundaryConfig.h"
#include "../WCSPH_DBC/OpenBoundaryConditions.h"
#include "Variables.h"
#include <TNL/Matrices/StaticMatrix.h>
#include "../EquationOfState.h"
#include "../DiffusiveTerms.h"
#include "../VisousTerms.h"
#include "control.h"

namespace TNL {
namespace SPH {

template< typename Particles, typename ModelConfig >
class WCSPH_DBC
{
public:
   using Model = WCSPH_DBC< Particles, ModelConfig >;
   using ModelParams = WCSPH_DBCConfig< ModelConfig >;
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
   using Matrix = Matrices::StaticMatrix< RealType, SPHConfig::spaceDimension + 1, SPHConfig::spaceDimension + 1 >;
   using VectorExtendedType = Containers::StaticVector< SPHConfig::spaceDimension + 1, RealType >;

   using FluidVariables = FluidVariables< SPHConfig >;
   using BoundaryVariables = BoundaryVariables< ModelConfig >;
   using OpenBoundaryVariables = OpenBoundaryVariables< SPHConfig >;
   using IntegrationSchemeType = typename ModelConfig::IntegrationScheme;
   using IntegrationSchemeVariables = typename IntegrationSchemeType::IntegrationSchemeVariablesType;
   using KernelFunction = typename ModelConfig::KernelFunction;
   using DiffusiveTerm = typename ModelConfig::DiffusiveTerm;
   using ViscousTerm = typename ModelConfig::ViscousTerm;
   using EOS = typename ModelConfig::EOS;

   using OpenBoundaryConfig = DBCOpenBoundaryConfig< SPHConfig >;
   using OpenBoundaryModel = OpenBoundaryConditionsBuffers< SPHConfig, ModelConfig >;

   static std::string writeModelType() { return "TNL::SPH::WCSPH_DBC"; }
   WCSPH_DBC() = default;

   template< typename EOS_ = EOS, typename PhysicalObjectPointer >
   void computePressureFromDensity( PhysicalObjectPointer& physicalObject, ModelParams& modelParams );

   template< typename FluidPointer, typename BoundaryPointer >
   void interaction( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams );

   template< typename FluidPointer, typename BoundaryPointer >
   void updateSolidBoundary( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams );

   template< typename FluidPointer, typename BoundaryPointer >
   void initializeInteraction( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams ) {}

   template< typename FluidPointer, typename BoundaryPointer >
   void finalizeInteraction( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams ) {}

   template< typename FluidPointer, typename BoundaryPointer >
   void initializeBoundaryInteraction( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams ) {}

   template< typename FluidPointer, typename BoundaryPointer >
   void finalizeBoundaryInteraction( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams ) {}

   template< typename OpenBoundaryPointer, typename BoundaryPointer >
   void updateSolidBoundaryOpenBoundary( BoundaryPointer& boundary, OpenBoundaryPointer& openBoundary, ModelParams& modelParams ) {}

   template< typename FluidPointer, typename OpenBoundaryPointer >
   void interactionWithOpenBoundary( FluidPointer& fluid, OpenBoundaryPointer& openBoundary, ModelParams& modelParams ) {}

   template< typename FluidPointer, typename OpenBoundaryPointer >
   void interactionWithBoundaryPatches( FluidPointer& fluid, OpenBoundaryPointer& boundaryPatch, ModelParams& modelParams ) {}

   template< typename FluidPointer, typename OpenBoundaryPointer >
   void extrapolateOpenBoundaryData( FluidPointer& fluid, OpenBoundaryPointer& openBoundary, ModelParams& modelParams, OpenBoundaryConfig& openBoundaryParams ) {}

};

} // SPH
} // TNL

#include "Interactions.hpp"

// end Interactions.h