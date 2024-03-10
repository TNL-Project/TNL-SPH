#pragma once

#include "../../SPHTraits.h"
#include "BoundaryConditionsTypes.h"
#include "OpenBoundaryConfig.h"
#include "SPH/Models/WCSPH_DBC/OpenBoundaryConditions.h"
#include "Variables.h"
#include "../../../Particles/neighborSearchLoop.h"
#include <TNL/Matrices/StaticMatrix.h>

/**
 * Modules used as default.
 **/
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

   /**
    * Constructor.
    */
   WCSPH_DBC( ) = default;

   /**
    * Print model identifier.
    */
   static std::string
   writeModelType()
   {
      return "TNL::SPH::WCSPH_DBC";
   }

   /**
    * Compute pressure from density.
    */
   template< typename EquationOfState = EquationsOfState::TaitWeaklyCompressibleEOS< SPHConfig >,
             typename PhysicalObjectPointer >
   void
   computePressureFromDensity( PhysicalObjectPointer& physicalObject, ModelParams& modelParams );

   /**
    * Function to realize fluid-fluid and fluid-boundary interaction.
    */
   template< typename FluidPointer, typename BoudaryPointer >
   void
   interaction( FluidPointer& fluid, BoudaryPointer& boundary, ModelParams& modelParams );

   /**
    * Function to realize boundary conditions for solid wall.
    * Realized by Dynamic Boundary Conditions (DBC) - Crespo et. 2007
    */
   template< typename FluidPointer,
             typename BoudaryPointer,
             typename BCType = typename ModelConfig::BCType,
             std::enable_if_t< std::is_same_v< BCType, WCSPH_BCTypes::DBC >, bool > Enabled = true >
   void
   updateSolidBoundary( FluidPointer& fluid, BoudaryPointer& boundary, ModelParams& modelParams );

   template< typename OpenBoundaryPointer,
             typename BoudaryPointer,
             typename BCType = typename ModelConfig::BCType,
             typename std::enable_if_t< std::is_same_v< BCType, WCSPH_BCTypes::DBC >, bool > Enabled = true >
   void
   updateSolidBoundaryOpenBoundary( BoudaryPointer& boundary, OpenBoundaryPointer& openBoundary, ModelParams& modelParams );

   /**
    * Function to realize boundary conditions for solid wall.
    * Realized by Modified Dynamic Boundary Conditions (MDBC) - English et. al. 2021
    */
   template< typename FluidPointer,
             typename BoudaryPointer,
             typename BCType = typename ModelConfig::BCType,
             typename std::enable_if_t< std::is_same_v< BCType, WCSPH_BCTypes::MDBC >, bool > Enabled = true >
   void
   updateSolidBoundary( FluidPointer& fluid, BoudaryPointer& boundary, ModelParams& modelParams );

   template< typename OpenBoundaryPointer,
             typename BoudaryPointer,
             typename BCType = typename ModelConfig::BCType,
             typename std::enable_if_t< std::is_same_v< BCType, WCSPH_BCTypes::MDBC >, bool > Enabled = true >
   void
   updateSolidBoundaryOpenBoundary( BoudaryPointer& boundary, OpenBoundaryPointer& openBoundary, ModelParams& modelParams );

   /**
    * Function to realize boundary conditions for solid wall.
    * Realized by Generalized Wall Boundary Conditions (GWBC) - Adami, Hu 2012
    *
    * TODO: Not implemented.
    */
   template< typename FluidPointer,
             typename BoudaryPointer,
             typename BCType = typename ModelConfig::BCType,
             typename std::enable_if_t< std::is_same_v< BCType, WCSPH_BCTypes::GWBC >, bool > Enabled = true >
   void
   updateSolidBoundary( FluidPointer& fluid, BoudaryPointer& boundary, ModelParams& modelParams );


   //TODO: Where should be this placed
   template< typename FluidPointer, typename OpenBoudaryPointer >
   void
   interactionWithOpenBoundary( FluidPointer& fluid, OpenBoudaryPointer& openBoundary, ModelParams& modelParams );

   //TODO: Experiment:
   template< typename FluidPointer, typename OpenBoudaryPointer >
   void
   interactionWithBoundaryPatches( FluidPointer& fluid, OpenBoudaryPointer& boundaryPatch, ModelParams& modelParams );

   /**
    * Functions to extrapolate data on open boundary buffers in 2D.
    * Extrapolation function is provided in 3 alternatives - extrapolation only density or velocity
    * or extrapolating both together.
    * By default, all variables are extrapolated using corrected interpolation - Liu et. al. 2006
    */
   template< typename FluidPointer, typename OpenBoudaryPointer >
   void
   extrapolateOpenBoundaryDensity2D( FluidPointer& fluid,
                                     OpenBoudaryPointer& openBoundary,
                                     ModelParams& modelParams,
                                     OpenBoundaryConfig& openBoundaryParams );

   template< typename FluidPointer, typename OpenBoudaryPointer >
   void
   extrapolateOpenBoundaryVelocity2D( FluidPointer& fluid,
                                      OpenBoudaryPointer& openBoundary,
                                      ModelParams& modelParams,
                                      OpenBoundaryConfig& openBoundaryParams );

   template< typename FluidPointer, typename OpenBoudaryPointer >
   void
   extrapolateOpenBoundaryData2D( FluidPointer& fluid,
                                  OpenBoudaryPointer& openBoundary,
                                  ModelParams& modelParams,
                                  OpenBoundaryConfig& openBoundaryParams );

   /**
    * Functions to extrapolate data on open boundary buffers in 3D.
    * Extrapolation function is provided in 3 alternatives - extrapolation only density or velocity
    * or extrapolating both together.
    * By default, all variables are extrapolated using corrected interpolation - Liu et. al. 2006
    */
   template< typename FluidPointer, typename OpenBoudaryPointer >
   void
   extrapolateOpenBoundaryDensity3D( FluidPointer& fluid,
                                     OpenBoudaryPointer& openBoundary,
                                     ModelParams& modelParams,
                                     OpenBoundaryConfig& openBoundaryParams );

   template< typename FluidPointer, typename OpenBoudaryPointer >
   void
   extrapolateOpenBoundaryVelocity3D( FluidPointer& fluid,
                                      OpenBoudaryPointer& openBoundary,
                                      ModelParams& modelParams,
                                      OpenBoundaryConfig& openBoundaryParams );

   template< typename FluidPointer, typename OpenBoudaryPointer >
   void
   extrapolateOpenBoundaryData3D( FluidPointer& fluid,
                                  OpenBoudaryPointer& openBoundary,
                                  ModelParams& modelParams,
                                  OpenBoundaryConfig& openBoundaryParams );

   /**
    * General function to perform extrapolation of open boundary conditions.
    */
   template< typename FluidPointer, typename OpenBoudaryPointer >
   void
   extrapolateOpenBoundaryData( FluidPointer& fluid,
                                OpenBoudaryPointer& openBoundary,
                                ModelParams& modelParams,
                                OpenBoundaryConfig& openBoundaryParams );

   template< typename FluidPointer, typename BoundaryPointer >
   void
   initializeInteraction( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams ) {}

   template< typename FluidPointer,
             typename BoundaryPointer,
             typename BCType = typename ModelConfig::BCType,
             typename std::enable_if_t< std::is_same_v< BCType, WCSPH_BCTypes::DBC >, bool > Enabled = true >
   void
   finalizeInteraction( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams ) {}

   template< typename FluidPointer,
             typename BoundaryPointer,
             typename BCType = typename ModelConfig::BCType,
             typename std::enable_if_t< std::is_same_v< BCType, WCSPH_BCTypes::MDBC >, bool > Enabled = true >
   void
   finalizeInteraction( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams );

   // experiment
   template< typename FluidPointer, typename BoudaryPointer, typename OpenBoundaryPointer >
   void
   interactWithPeriodicBoundary( FluidPointer& fluid,
                                 BoudaryPointer& boundary,
                                 OpenBoundaryPointer& openBoundary,
                                 ModelParams& modelParams,
                                 const VectorType shift );

   template< typename OpenBoundaryPointer,
             typename BoudaryPointer,
             typename FluidPointer,
             typename BCType = typename ModelConfig::BCType,
             typename std::enable_if_t< std::is_same_v< BCType, WCSPH_BCTypes::MDBC >, bool > Enabled = true >
   void
   updateSolidBoundaryPeriodicBoundary( FluidPointer& fluidPointer,
                                        BoudaryPointer& boundary,
                                        OpenBoundaryPointer& openBoundary,
                                        ModelParams& modelParams,
                                        const VectorType shift );

}; // SPH

} // SPH
} // TNL

#include "Interactions.hpp"
#include "BoundaryConditions/DBC.h"
#include "BoundaryConditions/MDBC.h"
//#include "OpenBoundaryConditionsInteractions.h"

#include "OpenBoundaryConditionsDataExtrapolation.h"

