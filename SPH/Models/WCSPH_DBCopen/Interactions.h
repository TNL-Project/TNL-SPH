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

template< typename Particles, typename OpenBoundary, typename SPHFluidConfig, typename Variables = SPHFluidVariables< SPHFluidConfig> >
class WCSPH_DBC
{
public:

   //Test
   using ParticleConfig = typename Particles::Config;
   using NeighborSearch = typename TNL::ParticleSystem::NeighborSearch< ParticleConfig, Particles >;
   using SPHSim = SPHOpenSystem< WCSPH_DBC< Particles, OpenBoundary, SPHFluidConfig >, Particles, NeighborSearch >;
   using SPHSimPointer = typename Pointers::SharedPointer< SPHSim, typename Particles::Device >;
   using Ptr = std::shared_ptr< SPHSim >;
   //Test

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
   using OpenBoundaryPointer = typename Pointers::SharedPointer< OpenBoundary, DeviceType >;

   /* VARIABLES FIELDS */
   using ScalarArrayType = typename SPHFluidTraitsType::ScalarArrayType;
   using VectorArrayType = typename SPHFluidTraitsType::VectorArrayType;
   using ParticleTypeArrayType = typename SPHFluidTraitsType::ParticleTypeArrayType;
   using EOS = TaitWeaklyCompressibleEOS< SPHFluidConfig >;

   /* Thrust sort */
   using IndexArrayType = Containers::Array< GlobalIndexType, DeviceType >;

   /* Integrator */
   using Model = WCSPH_DBC< Particles, OpenBoundary, SPHFluidConfig >;
   using Integrator = VerletIntegrator< typename Pointers::SharedPointer< Model, DeviceType >, SPHFluidConfig >;

   /*Swap variables*/
   using SwapVariables = SWAPFluidVariables< SPHFluidConfig, PointArrayType >; //REDO

   /**
    * Constructor.
    */
   WCSPH_DBC( GlobalIndexType size,
              ParticlePointer& particles,
              GlobalIndexType size_boundary,
              ParticlePointer& particles_boundary,
              GlobalIndexType size_inlet,
              OpenBoundaryPointer& openBoundary ) //THIS WORKS
   : fluidVariables( size ), boundaryVariables( size_boundary ), inletVariables( size_inlet ),
#ifdef PREFER_SPEED_OVER_MEMORY
     swapFluid( size ), swapBoundary( size_boundary ), swapInlet( size_inlet ),
#endif
   particles( particles ), particles_bound( particles_boundary ), openBoundary( openBoundary ) {} //THIS WORKS

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

   const Variables&
   getInletVariables() const;

   Variables&
   getInletVariables();

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
   sortParticlesAndVariablesThrust( ParticlePointer& particles, Variables& variables, SwapVariables& variables_swap );

   void
   sortBoundaryParticlesAndVariablesThrust( ParticlePointer& particles, Variables& variables, SwapVariables& variables_swap );

   //void
   //sortInletParticlesAndVariablesThrust( ParticlePointer& particles, Variables& variables, SwapVariables& variables_swap );

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
   Variables fluidVariables;
   Variables boundaryVariables;
   Variables inletVariables;

   ParticlePointer particles;
   ParticlePointer particles_bound;

   SwapVariables swapFluid;
   SwapVariables swapBoundary;
   SwapVariables swapInlet;

   OpenBoundaryPointer openBoundary;
   //SPHSim simulation;
   //std::shared_ptr< SPHSim > simulation;
   //std::weak_ptr< SPHSim > simulation;


#ifdef PREFER_SPEED_OVER_MEMORY
   /* TEMP - Indices for thrust sort. */

#endif

};

} // SPH
} // ParticleSystem
} // TNL

#include "Interactions_impl.h"
#include "Interactions_implLambda.h"

