#pragma once

namespace TNL {
namespace SPH {

template< typename Particles, typename ModelConfig >
class SPHMultisetSolverTemplate
{
public:

   /**
    * Constructor.
    */
   SPHMultisetSolverTemplate() = default;

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
   template< typename FluidPointer, typename BoudaryPointer, typename ModelParams >
   void
   interaction( FluidPointer& fluid, BoudaryPointer& boundary, ModelParams& modelParams )
   {}

   template< typename FluidPointer, typename BoudaryPointer, typename ModelParams >
   void
   updateSolidBoundary( FluidPointer& fluid, BoudaryPointer& boundary, ModelParams& modelParams )
   {}

   template< typename OpenBoundaryPointer, typename BoudaryPointer, typename ModelParams >
   void
   updateSolidBoundaryOpenBoundary( BoudaryPointer& boundary, OpenBoundaryPointer& openBoundaryPointer, ModelParams& modelParams )
   {}

   template< typename FluidPointer, typename OpenBoudaryPointer, typename ModelParams >
   void
   interactionWithOpenBoundary( FluidPointer& fluid, OpenBoudaryPointer& openBoundary, ModelParams& modelParams )
   {}

   template< typename FluidPointer, typename OpenBoudaryPointer, typename ModelParams >
   void
   interactionWithBoundaryPatches( FluidPointer& fluid, OpenBoudaryPointer& openBoundary, ModelParams& modelParams )
   {}

   ///**
   // * General function to perform extrapolation of open boundary conditions.
   // */
   //template< typename FluidPointer, typename OpenBoudaryPointer, typename ModelParams, typename OpenBoundaryConfig  >
   //void
   //extrapolateOpenBoundaryData( FluidPointer& fluid, OpenBoudaryPointer& openBoundary, ModelParams& modelParams, OpenBoundaryConfig& openBoundaryParams )
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

}  //namespace SPH
}  //namespace TNL

