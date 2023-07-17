#pragma once

#include "../../SPHTraits.h"
#include "Variables.h"
#include "../../../Particles/neighborSearchLoop.h"

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
   using VariablesPointer = typename Pointers::SharedPointer< Variables, DeviceType >;

   using ParticlesType = Particles;

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
   ComputePressureFromDensity( PhysicalObjectPointer& physicalObject, SPHState& sphState );

   template< typename FluidPointer,
             typename BoudaryPointer,
             typename SPHKernelFunction,
             typename DiffusiveTerm,
             typename ViscousTerm,
             typename EOS,
             typename SPHState  >
   void
   Interaction( FluidPointer& fluid, BoudaryPointer& boundary, SPHState& sphState );

};

} // SPH
} // ParticleSystem
} // TNL

#include "Interactions_impl.h"

