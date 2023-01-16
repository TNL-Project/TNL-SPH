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
class WCSPH_DBC
{
public:

   using SPHConfig = SPHFluidConfig;
   using SPHFluidTraitsType = SPHFluidTraits< SPHFluidConfig >;
   using DeviceType = typename Particles::Device; //?

   using LocalIndexType = typename SPHFluidTraitsType::LocalIndexType; //Particles::?
   using GlobalIndexType = typename SPHFluidTraitsType::GlobalIndexType; //Particles::?
   using RealType = typename SPHFluidTraitsType::RealType; //Particles::?

   using PointType = typename Particles::PointType;
   using PointArrayType = typename Particles::PointArrayType;

   using ScalarType = typename SPHFluidTraitsType::ScalarType;
   using VectorType = typename SPHFluidTraitsType::VectorType;

   using DiffusiveTerm = MolteniDiffusiveTerm< SPHFluidConfig >; //-> template
   using ViscousTerm = ArtificialViscosity< SPHFluidConfig >; //-> template

   using ParticlePointer = typename Pointers::SharedPointer< Particles, DeviceType >;

   /* VARIABLES FIELDS */
   using ScalarArrayType = typename SPHFluidTraitsType::ScalarArrayType;
   using VectorArrayType = typename SPHFluidTraitsType::VectorArrayType;
   using EOS = TaitWeaklyCompressibleEOS< SPHFluidConfig >;

   /* Thrust sort */
   using IndexArrayType = Containers::Array< GlobalIndexType, DeviceType >;
   using IndexArrayTypePointer = typename Pointers::SharedPointer< IndexArrayType, DeviceType >;

   /* Integrator */
   using Model = WCSPH_DBC< Particles, SPHFluidConfig >;
   using Integrator = VerletIntegrator< typename Pointers::SharedPointer< Model, DeviceType >, SPHFluidConfig >;
   using IntegratorVariables = IntegratorVariables< SPHFluidConfig >;

   /*Swap variables*/
   using SwapVariables = SWAPFluidVariables< SPHFluidConfig, PointArrayType >; //REDO
   using ModelVariables = Variables;
   using VariablesPointer = typename Pointers::SharedPointer< ModelVariables, DeviceType >;

   /**
    * Constructor.
    */
   WCSPH_DBC( GlobalIndexType size,
              GlobalIndexType size_boundary,
              GlobalIndexType size_inlet ) //THIS WORKS
   : swapFluid( size ), swapBoundary( size_boundary ), swapInlet( size_inlet ) {} //THIS WORKS


   /**
    * Get fileds with variables.
    */
   const IndexArrayType&
   getIndicesForReoder() const;

   IndexArrayType&
   getIndicesForReoder();

   /**
    * Sort particles and all variables based on particle cell index.
    * TODO: Move this on the side of nbsearch/particles.
    */
   void
   sortParticlesAndVariablesThrust( ParticlePointer& particles, VariablesPointer& variables, SwapVariables& variables_swap );

   void
   sortBoundaryParticlesAndVariablesThrust( ParticlePointer& particles, VariablesPointer& variables, SwapVariables& variables_swap );

   void
   sortVariables( ParticlePointer& particles, VariablesPointer& variables, IndexArrayTypePointer& map );

   void
   sortVariablesInPlace( ParticlePointer& particles, VariablesPointer& variables );

   //void
   //sortInletParticlesAndVariablesThrust( ParticlePointer& particles, Variables& variables, SwapVariables& variables_swap );

   /**
    * Compute pressure from density.
    * TODO: Move out.
    */
   //USETHIS:template< typename EquationOfState = TaitWeaklyCompressibleEOS< SPHFluidConfig > >
   //USETHIS:void
   //USETHIS:ComputePressureFromDensity();

   template< typename FluidPointer, typename BoudaryPointer, typename OpenBoundaryPointer, typename NeighborSearchPointer, typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm, typename EOS  >
   void
   Interaction( FluidPointer& fluid, BoudaryPointer& boundary, OpenBoundaryPointer& openBoundary );

   /* Constants */ //Move to protected
   RealType h, m, speedOfSound, coefB, rho0, delta, alpha;

   SwapVariables swapFluid;
   SwapVariables swapBoundary;
   SwapVariables swapInlet;



#ifdef PREFER_SPEED_OVER_MEMORY
   /* TEMP - Indices for thrust sort. */

#endif

};

} // SPH
} // ParticleSystem
} // TNL

#include "Interactions_impl.h"
#include "Interactions_implLambda.h"

