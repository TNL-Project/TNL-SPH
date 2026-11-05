// SPDX-FileComment: This file is part of TNL-SPH (https://gitlab.com/tnl-project/tnl-sph)
// SPDX-License-Identifier: MIT

#pragma once

namespace TNL::SPH {

template< typename Particles, typename ModelConfig >
class SPHMultiSetSolverTemplate
{
public:
   /**
    * Constructor.
    */
   SPHMultiSetSolverTemplate() = default;

   /**
    * Print model identifier.
    */
   static std::string
   writeModelType()
   {
      return "TNL::SPH::NOT_DEFINED";
   }

   /**
    * Function to realize fluid-fluid and fluid-boundary interaction.
    */
   template< typename FluidPointer, typename BoundaryPointer, typename ModelParams >
   void
   interaction( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams )
   {}

   template< typename FluidPointer, typename BoundaryPointer, typename ModelParams >
   void
   updateSolidBoundary( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams )
   {}

   template< typename OpenBoundaryPointer, typename BoundaryPointer, typename ModelParams >
   void
   updateSolidBoundaryOpenBoundary(
      BoundaryPointer& boundary,
      OpenBoundaryPointer& openBoundaryPointer,
      ModelParams& modelParams )
   {}

   template< typename FluidPointer, typename OpenBoundaryPointer, typename ModelParams >
   void
   interactionWithOpenBoundary( FluidPointer& fluid, OpenBoundaryPointer& openBoundary, ModelParams& modelParams )
   {}

   template< typename FluidPointer, typename OpenBoundaryPointer, typename ModelParams >
   void
   interactionWithBoundaryPatches( FluidPointer& fluid, OpenBoundaryPointer& openBoundary, ModelParams& modelParams )
   {}

   ///**
   // * General function to perform extrapolation of open boundary conditions.
   // */
   //template< typename FluidPointer, typename OpenBoundaryPointer, typename ModelParams, typename OpenBoundaryConfig  >
   //void
   //extrapolateOpenBoundaryData( FluidPointer& fluid, OpenBoundaryPointer& openBoundary, ModelParams& modelParams,
   //OpenBoundaryConfig& openBoundaryParams )
   //{}

   template< typename FluidPointer, typename BoundaryPointer, typename ModelParams >
   void
   initializeInteraction( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams )
   {}

   template< typename FluidPointer, typename BoundaryPointer, typename ModelParams >
   void
   finalizeInteraction( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams )
   {}

   template< typename FluidPointer, typename BoundaryPointer, typename ModelParams >
   void
   finalizeBoundaryInteraction( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams )
   {}

   // mm?
   template< typename PhysicalObjectPointer, typename ModelParams >
   void
   computePressureFromDensity( PhysicalObjectPointer& physicalObject, ModelParams& modelParams )
   {}
};

}  // namespace TNL::SPH
