#pragma once

#include "../../SPHTraits.h"
#include "Variables.h"
#include "../../../Particles/neighborSearchLoop.h"
#include "BoundaryConditionsTypes.h"
#include "OpenBoundaryConfig.h"

/**
 * Modules used as default.
 **/
#include "../EquationOfState.h"
#include "../DiffusiveTerms.h"
#include "../VisousTerms.h"
#include "Integrator.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Particles, typename SPHState >
class WCSPH_BI
{
public:

   using SPHConfig = typename SPHState::SPHConfig;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using DeviceType = typename SPHConfig::DeviceType; //TODO:Resolve

   using LocalIndexType = typename SPHTraitsType::LocalIndexType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using ScalarType = typename SPHTraitsType::ScalarType;
   using VectorType = typename SPHTraitsType::VectorType;
   using IndexVectorType = typename SPHTraitsType::IndexVectorType;

   /* VARIABLES FIELDS */
   using EOS = TaitWeaklyCompressibleEOS< SPHConfig >;

   /* Integrator */
   using Model = WCSPH_BI< Particles, SPHConfig >;
   using Integrator = VerletIntegrator< SPHConfig >;
   using IntegratorVariables = IntegratorVariables< SPHConfig >;

   /*Swap variables*/
   using FluidVariables = SPHFluidVariables< SPHConfig >;
   using BoundaryVariables = SPHBoundaryVariables< SPHConfig >;
   using OpenBoundaryVariables = SPHOpenBoundaryVariables< SPHConfig >;

   using ParticlesType = Particles;

   //Open boundary
   using OpenBoundaryConfig = BIOpenBoundaryConfig< SPHConfig >;

   /**
    * Constructor.
    */
   WCSPH_BI( ) = default;

   /**
    * Compute pressure from density.
    * TODO: Move out.
    */
   template< typename EquationOfState = TaitWeaklyCompressibleEOS< SPHConfig >,
             typename PhysicalObjectPointer >
   void
   computePressureFromDensity( PhysicalObjectPointer& physicalObject, SPHState& sphState );

   /**
    * Function to realize fluid-fluid and fluid-boundary interaction.
    */
   template< typename FluidPointer,
             typename BoudaryPointer,
             typename SPHKernelFunction,
             typename DiffusiveTerm,
             typename ViscousTerm,
             typename EOS >
   void
   interaction( FluidPointer& fluid, BoudaryPointer& boundary, SPHState& sphState );

   template< typename FluidPointer,
             typename BoudaryPointer,
             typename SPHKernelFunction,
             typename DiffusiveTerm,
             typename ViscousTerm,
             typename EOS >
   void
   updateSolidBoundary( FluidPointer& fluid, BoudaryPointer& boundary, SPHState& sphState );

   //TODO: Experiment:
   template< typename FluidPointer,
             typename BoundaryPointer,
             typename OpenBoudaryPointer,
             typename SPHKernelFunction,
             typename DiffusiveTerm,
             typename ViscousTerm,
             typename EOS >
   void
   interactionWithOpenBoundary( FluidPointer& fluid,
                                BoundaryPointer& boundary,
                                OpenBoudaryPointer& openBoundary,
                                SPHState& sphState );

   template< typename FluidPointer, typename BoundaryPointer >
   void
   finalizeInteraction( FluidPointer& fluid, BoundaryPointer& boundary, SPHState& sphState );

};

} // SPH
} // ParticleSystem
} // TNL

#include "Interactions.hpp"

