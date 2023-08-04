#pragma once

#include "../../SPHTraits.h"
#include "Variables.h"
#include "../../../Particles/neighborSearchLoop.h"
#include <TNL/Matrices/StaticMatrix.h>

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

template< typename Particles, typename SPHFluidConfig, typename Variables = SPHFluidVariables< SPHFluidConfig> >
class WCSPH_DBC
{
public:

   using SPHConfig = SPHFluidConfig;
   using SPHTraitsType = SPHFluidTraits< SPHFluidConfig >;
   using DeviceType = typename SPHConfig::DeviceType; //TODO:Resolve

   using LocalIndexType = typename SPHTraitsType::LocalIndexType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using ScalarType = typename SPHTraitsType::ScalarType;
   using VectorType = typename SPHTraitsType::VectorType;
   using IndexVectorType = typename SPHTraitsType::IndexVectorType;

   /* VARIABLES FIELDS */
   using EOS = TaitWeaklyCompressibleEOS< SPHFluidConfig >;

   /* Integrator */
   using Model = WCSPH_DBC< Particles, SPHFluidConfig >;
   using Integrator = VerletIntegrator< typename Pointers::SharedPointer< Model, DeviceType >, SPHFluidConfig >;
   using IntegratorVariables = IntegratorVariables< SPHFluidConfig >;

   /*Swap variables*/
	using FluidVariables = Variables;
	using BoundaryVariables = Variables;
   using OpenBoundaryVariables = SPHOpenBoundaryVariables< SPHFluidConfig >;

   using ParticlesType = Particles;

   //Open boundary
   using Matrix = Matrices::StaticMatrix< RealType, SPHConfig::spaceDimension + 1, SPHConfig::spaceDimension + 1 >;
   //using Matrix = Matrices::StaticMatrix< RealType, SPHConfig::spaceDimension, SPHConfig::spaceDimension >;
   using VectorExtendedType = Containers::StaticVector< SPHConfig::spaceDimension + 1, RealType >;

   /**
    * Constructor.
    */
   WCSPH_DBC( ) = default; //THIS WORKS

   /**
    * Compute pressure from density.
    * TODO: Move out.
    */
   template< typename EquationOfState = TaitWeaklyCompressibleEOS< SPHFluidConfig >,
             typename PhysicalObjectPointer,
             typename SPHState >
   void
   computePressureFromDensity( PhysicalObjectPointer& physicalObject, SPHState& sphState );

   template< typename FluidPointer,
             typename BoudaryPointer,
             typename SPHKernelFunction,
             typename DiffusiveTerm,
             typename ViscousTerm,
             typename EOS,
             typename SPHState  >
   void
   interaction( FluidPointer& fluid, BoudaryPointer& boundary, SPHState& sphState );

   //TODO: Remove, testing.
   template< typename FluidPointer,
             typename BoudaryPointer,
             typename SPHKernelFunction,
             typename DiffusiveTerm,
             typename ViscousTerm,
             typename EOS,
             typename SPHState  >
   RealType
   interactionWithReduction( FluidPointer& fluid, BoudaryPointer& boundary, SPHState& sphState );

   //TODO: Where should be this placed
   template< typename FluidPointer,
             typename OpenBoudaryPointer,
             typename SPHKernelFunction,
             typename DiffusiveTerm,
             typename ViscousTerm,
             typename EOS,
             typename SPHState  >
   void
   interactionWithOpenBoundary( FluidPointer& fluid, OpenBoudaryPointer& openBoundary, SPHState& sphState );

   //TODO: Experiment:
   template< typename FluidPointer,
             typename BoundaryPointer,
             typename OpenBoudaryPointer,
             typename SPHKernelFunction,
             typename DiffusiveTerm,
             typename ViscousTerm,
             typename EOS,
             typename SPHState  >
   void
   interactionWithOpenBoundary( FluidPointer& fluid,
                                BoundaryPointer& boundary,
                                OpenBoudaryPointer& openBoundary,
                                SPHState& sphState );

   //TODO: Figure this out:
   template< typename FluidPointer,
             typename OpenBoudaryPointer,
             typename SPHKernelFunction,
             typename EOS,
             typename SPHState  >
   void
   extrapolateOpenBoundaryData( FluidPointer& fluid,
                                OpenBoudaryPointer& openBoundary,
                                SPHState& sphState );

};

} // SPH
} // ParticleSystem
} // TNL

#include "Interactions.hpp"
#include "OpenBoundaryContirionsInteractions.h"

