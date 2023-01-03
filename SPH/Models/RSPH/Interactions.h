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
   using ParticleTypeArrayType = typename SPHFluidTraitsType::ParticleTypeArrayType;
   using EOS = TaitWeaklyCompressibleEOS< SPHFluidConfig >;

   /* Thrust sort */
   using IndexArrayType = Containers::Array< GlobalIndexType, DeviceType >;

   /* Integrator */
   using Model = RSPHSimple< Particles, SPHFluidConfig >;
   using Integrator = VerletIntegrator< typename Pointers::SharedPointer< Model, DeviceType >, SPHFluidConfig >;

   /*Swap variables*/
   using SwapVariables = SWAPFluidVariables< SPHFluidConfig, PointArrayType >; //REDO

   /**
    * Constructor.
    */
   RSPHSimple( GlobalIndexType size, ParticlePointer& particles, GlobalIndexType size_boundary, ParticlePointer& particles_boundary )
   : FluidVariables( size ), BoundaryVariables( size_boundary ),
#ifdef PREFER_SPEED_OVER_MEMORY
     swapFluid( size ), swapBoundary( size_boundary ),
#endif
   particles( particles ), boundaryParticles( particles_boundary ) {
   }

   /**
    * Get fileds with variables.
    */
   const Variables&
   getFluidVariables() const;

   Variables&
   getFluidVariables();

   const Variables&
   getBoundaryVariables() const;

   Variables&
   getBoundaryVariables();

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
   void sortParticlesAndVariablesThrust( ParticlePointer& particles, Variables& variables, SwapVariables& variables_swap );

   /**
    * Compute pressure from density.
    * TODO: Move out.
    */
   template< typename EquationOfState = TaitWeaklyCompressibleEOS< SPHFluidConfig > >
   void
   ComputePressureFromDensity();

   template< typename NeighborSearchPointer, typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm, typename EOS >
   void
   Interaction( NeighborSearchPointer& neighborSearch, NeighborSearchPointer& neighborSearch_bound );

   /* Constants */ //Move to protected
   RealType h, m, speedOfSound, coefB, rho0, delta, alpha;

//protected:

   /* Variables - Fields */
   Variables FluidVariables;
   Variables BoundaryVariables;

   ParticlePointer particles;
   ParticlePointer boundaryParticles;

   SwapVariables swapFluid;
   SwapVariables swapBoundary;

#ifdef PREFER_SPEED_OVER_MEMORY
   /* TEMP - Indices for thrust sort. */

#endif

};

} // SPH
} // ParticleSystem
} // TNL

#include "Interactions_impl.h"
#include "Interactions_implLambda.h"

