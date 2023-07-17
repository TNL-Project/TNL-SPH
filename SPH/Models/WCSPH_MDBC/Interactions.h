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

//#include "../../myStaticMatrix.h"
#include <TNL/Matrices/StaticMatrix.h>

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Particles, typename SPHFluidConfig, typename Variables = SPHFluidVariables< SPHFluidConfig> >
class WCSPH_MDBC
{
public:

   using SPHConfig = SPHFluidConfig;
   using SPHFluidTraitsType = SPHFluidTraits< SPHFluidConfig >;
   using DeviceType = typename SPHConfig::DeviceType; //TODO:Resolve

   using LocalIndexType = typename SPHFluidTraitsType::LocalIndexType;
   using GlobalIndexType = typename SPHFluidTraitsType::GlobalIndexType;
   using RealType = typename SPHFluidTraitsType::RealType;
   using ScalarType = typename SPHFluidTraitsType::ScalarType;
   using VectorType = typename SPHFluidTraitsType::VectorType;
   using IndexVectorType = typename SPHFluidTraitsType::IndexVectorType;

   //Model specific types
   using Matrix = Matrices::StaticMatrix< RealType, SPHConfig::spaceDimension + 1, SPHConfig::spaceDimension + 1 >;
   using VectorExtendedType = Containers::StaticVector< SPHConfig::spaceDimension + 1, RealType >;


   /* VARIABLES FIELDS */
   using EOS = TaitWeaklyCompressibleEOS< SPHFluidConfig >;

   /* Integrator */
   using Model = WCSPH_MDBC< Particles, SPHFluidConfig >;
   using Integrator = VerletIntegrator< typename Pointers::SharedPointer< Model, DeviceType >, SPHFluidConfig >;
   using IntegratorVariables = IntegratorVariables< SPHFluidConfig >;

   /*Swap variables*/
	using FluidVariables = Variables;
	using BoundaryVariables = SPHBoundaryVariables< SPHConfig >;
   using VariablesPointer = typename Pointers::SharedPointer< Variables, DeviceType >;

   using ParticlesType = Particles;

   /**
    * Constructor.
    */
   WCSPH_MDBC() = default; //THIS WORKS

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

};

} // SPH
} // ParticleSystem
} // TNL

#include "Interactions.hpp"

