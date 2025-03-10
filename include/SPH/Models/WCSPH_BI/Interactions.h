#pragma once

#include "../../SPHTraits.h"
#include "BoundaryConditionsTypes.h"
#include "OpenBoundaryConfig.h"
#include "SPH/Models/WCSPH_BI/OpenBoundaryConditions.h"
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
class WCSPH_BI
{
public:
   using Model = WCSPH_BI< Particles, ModelConfig >;
   using ModelParams = WCSPH_BIConfig< ModelConfig >;
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
   using BoundaryViscousTerm = typename ModelConfig::BoundaryViscousTerm;
   using EOS = typename ModelConfig::EOS;
   using DensityFilter = typename ModelConfig::DensityFilter;

   using OpenBoundaryConfig = BIOpenBoundaryConfig< SPHConfig >;
   using OpenBoundaryModel = OpenBoundaryConditionsBuffers< SPHConfig, ModelConfig >;

   /**
    * Constructor.
    */
   WCSPH_BI() = default;

   /**
    * Print model identifier.
    */
   static std::string
   writeModelType()
   {
      return "TNL::SPH::WCSPH_BI";
   }

   /**
    * Compute pressure from density.
    * TODO: Move out.
    */
   template< typename EquationOfState = EOS,
             typename PhysicalObjectPointer >
   void
   computePressureFromDensity( PhysicalObjectPointer& physicalObject, ModelParams& modelParams );

   template< typename FluidPointer >
   void
   filterDensity( FluidPointer& fluid, ModelParams& modelParams );

   /**
    * Function to realize fluid-fluid and fluid-boundary interaction.
    */
   template< typename FluidPointer, typename BoudaryPointer >
   void
   interaction( FluidPointer& fluid, BoudaryPointer& boundary, ModelParams& modelParams );

   template< typename FluidPointer, typename BoudaryPointer >
   void
   updateSolidBoundary( FluidPointer& fluid, BoudaryPointer& boundary, ModelParams& modelParams );

   template< typename OpenBoundaryPointer, typename BoudaryPointer >
   void
   updateSolidBoundaryOpenBoundary( BoudaryPointer& boundary,
                                    OpenBoundaryPointer& openBoundaryPointer,
                                    ModelParams& modelParams );

   template< typename FluidPointer, typename OpenBoudaryPointer >
   void
   interactionWithOpenBoundary( FluidPointer& fluid, OpenBoudaryPointer& openBoundary, ModelParams& modelParams );

   template< typename FluidPointer, typename OpenBoudaryPointer >
   void
   interactionWithBoundaryPatches( FluidPointer& fluid, OpenBoudaryPointer& openBoundary, ModelParams& modelParams );

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
   initializeInteraction( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams )
   {}

   template< typename FluidPointer, typename BoundaryPointer >
   void
   finalizeInteraction( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams );

   template< typename FluidPointer, typename BoundaryPointer >
   void
   finalizeBoundaryInteraction( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams );
};

}  //namespace SPH
}  //namespace TNL

#include "Interactions.hpp"

#include "OpenBoundaryConditionsDataExtrapolation.h"

