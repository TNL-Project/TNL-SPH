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
//#include "Integrator.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Particles, typename SPHFluidConfig, typename Variables = SPHFluidVariables< SPHFluidConfig> >
class WCSPH_DBC
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

   /**
    * Constructor.
    */
   WCSPH_DBC( GlobalIndexType size, ParticlePointer& particles )
   : FluidVariables( size ) , rhoO( size ), vO( size ),
#ifdef PREFER_SPEED_OVER_MEMORY
   indicesMap( size ), rho_swap( size ), v_swap( size ), rhoO_swap( size ), vO_swap( size ), points_swap( size ),
#endif
   particles( particles )
   {
      vO = 0.; //remove
      rhoO = 1000.; //remove
   }

   /**
    * Get fileds with variables.
    */
   const Variables&
   getFluidVariables() const;

   Variables&
   getFluidVariables();

   /**
    * Sort particles and all variables based on particle cell index.
    * TODO: Move this on the side of nbsearch/particles.
    */
   void sortParticlesAndVariables();
   void sortParticlesAndVariablesThrust();

   /**
    * Compute pressure from density.
    * TODO: Move out.
    */
   template< typename EquationOfState = TaitWeaklyCompressibleEOS< SPHFluidConfig > >
   void
   ComputePressureFromDensity();

   /* TEMP INTEGRATORS, MOVE OUT */
   void
   IntegrateVerlet( RealType dt );

   void
   IntegrateEuler( RealType dt );

   template< typename NeighborSearchPointer, typename ModelPointer, typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm, typename EOS >
   void
   Interaction( NeighborSearchPointer& neighborSearch, NeighborSearchPointer& neighborSearch_bound, ModelPointer& boundary );


   /* Constants */ //Move to protected
   RealType h, m, speedOfSound, coefB, rho0, delta, alpha;

protected:

   /* Variables - Fields */
   Variables FluidVariables;
   //Variables BoundaryVariables

   /* TEMP INTEGRATORS, MOVE OUT */
   ScalarArrayType rhoO;
   VectorArrayType vO;


   ParticlePointer particles;
   //ParticlePointer boundaryParticles;

   //Integrator integrator; //temp

#ifdef PREFER_SPEED_OVER_MEMORY
   /* TEMP - Indices for thrust sort. */
   IndexArrayType indicesMap;

   /* Copy of arrays needed for resort */
   ScalarArrayType rho_swap;
   VectorArrayType v_swap;
   PointArrayType points_swap;

   ScalarArrayType rhoO_swap;
   VectorArrayType vO_swap;
#endif

};

} // SPH
} // ParticleSystem
} // TNL

#include "Interactions_impl.h"
#include "Interactions_implLambda.h"

