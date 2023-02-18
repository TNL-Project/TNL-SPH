#pragma once

#include <TNL/Containers/Array.h>
#include <TNL/Containers/ArrayView.h>

/**
 * Use thrust for sorting.
 **/
#include <thrust/sort.h>
#include <thrust/gather.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/zip_iterator.h>

#include "../../SPHTraits.h"
#include "Variables.h"

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
class RSPHSimple
{
public:

   using SPHConfig = SPHFluidConfig;
   using SPHFluidTraitsType = SPHFluidTraits< SPHFluidConfig >;
   using DeviceType = typename Particles::Device;

   using LocalIndexType = typename SPHFluidTraitsType::LocalIndexType;
   using GlobalIndexType = typename SPHFluidTraitsType::GlobalIndexType;
   using RealType = typename SPHFluidTraitsType::RealType;
   using PointType = typename Particles::PointType;
   using ScalarType = typename SPHFluidTraitsType::ScalarType;
   using VectorType = typename SPHFluidTraitsType::VectorType;
   using IndexVectorType = typename SPHFluidTraitsType::IndexVectorType;

   using ParticlePointer = typename Pointers::SharedPointer< Particles, DeviceType >;

   /* VARIABLES FIELDS */
   using ScalarArrayType = typename SPHFluidTraitsType::ScalarArrayType;
   using VectorArrayType = typename SPHFluidTraitsType::VectorArrayType;
   using EOS = TaitWeaklyCompressibleEOS< SPHFluidConfig >;

   /* Thrust sort */
   using IndexArrayType = Containers::Array< GlobalIndexType, DeviceType >;
   using IndexArrayTypePointer = typename Pointers::SharedPointer< IndexArrayType, DeviceType >;

   /* Integrator */
   using Model = RSPHSimple< Particles, SPHFluidConfig >;
   using Integrator = VerletIntegrator< typename Pointers::SharedPointer< Model, DeviceType >, SPHFluidConfig >;
   using IntegratorVariables = IntegratorVariables< SPHFluidConfig >;

   /*Swap variables*/
   using ModelVariables = Variables;
   using VariablesPointer = typename Pointers::SharedPointer< ModelVariables, DeviceType >;

   /**
    * Constructor.
    */
   RSPHSimple() = default;

   /**
    * Compute pressure from density.
    * TODO: Move out.
    */
   template< typename EquationOfState = TaitWeaklyCompressibleEOS< SPHFluidConfig > >
   void
   ComputePressureFromDensity( VariablesPointer& variables, GlobalIndexType numberOfParticles );

   template< typename FluidPointer, typename BoudaryPointer, typename NeighborSearchPointer, typename RiemannSolver, typename EOS >
   void
   Interaction( FluidPointer& fluid, BoudaryPointer& boundary );

   /* Constants */ //Move to protected
   RealType h, m, speedOfSound, coefB, rho0, delta, alpha;
   VectorType g;

};

} // SPH
} // ParticleSystem
} // TNL

#include "Interactions_impl.h"
#include "Interactions_implLambda.h"

