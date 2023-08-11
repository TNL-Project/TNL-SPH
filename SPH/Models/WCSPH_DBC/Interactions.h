#pragma once

#include "../../SPHTraits.h"
#include "BoundaryConditionsTypes.h"
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

//#include "BoundaryConditionsTypes.h"
//#include "./BoundaryConditions/DBC.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Particles, typename SPHState >
class WCSPH_DBC
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
   using Model = WCSPH_DBC< Particles, SPHState >;
   using Integrator = VerletIntegrator< SPHConfig >;
   using IntegratorVariables = IntegratorVariables< SPHConfig >;

   /*Swap variables*/
	using FluidVariables = SPHFluidVariables< SPHConfig >;
	using BoundaryVariables = SPHBoundaryVariables< SPHState >;
   using OpenBoundaryVariables = SPHOpenBoundaryVariables< SPHConfig >;

   using ParticlesType = Particles;


   //Open boundary
   using Matrix = Matrices::StaticMatrix< RealType, SPHConfig::spaceDimension + 1, SPHConfig::spaceDimension + 1 >;
   //using Matrix = Matrices::StaticMatrix< RealType, SPHConfig::spaceDimension, SPHConfig::spaceDimension >;
   using VectorExtendedType = Containers::StaticVector< SPHConfig::spaceDimension + 1, RealType >;

   //using DefinedBCType = typename SPHState::BCType;


   /**
    * Constructor.
    */
   WCSPH_DBC( ) = default;

   /**
    * Compute pressure from density.
    * TODO: Move out.
    */
   template< typename EquationOfState = TaitWeaklyCompressibleEOS< SPHConfig >,
             typename PhysicalObjectPointer >
   void
   computePressureFromDensity( PhysicalObjectPointer& physicalObject, SPHState& sphState );

   template< typename FluidPointer,
             typename BoudaryPointer,
             typename SPHKernelFunction,
             typename DiffusiveTerm,
             typename ViscousTerm,
             typename EOS >
   void
   interaction( FluidPointer& fluid, BoudaryPointer& boundary, SPHState& sphState );

   /**
    * Function to realize boundary conditions for solid wall.
    * Realized by Dynamic Boundary Conditions (DBC) - Crespo et. 2007
    */
   template< typename FluidPointer,
             typename BoudaryPointer,
             typename SPHKernelFunction,
             typename DiffusiveTerm,
             typename ViscousTerm,
             typename EOS,
             typename BCType = typename SPHState::BCType,
             std::enable_if_t< std::is_same_v< BCType, WCSPH_BCTypes::DBC >, bool > Enabled = true >
   void
   updateSolidBoundary( FluidPointer& fluid, BoudaryPointer& boundary, SPHState& sphState );

   /**
    * Function to realize boundary conditions for solid wall.
    * Realized by Modified Dynamic Boundary Conditions (MDBC) - English et. al. 2021
    */
   template< typename FluidPointer,
             typename BoudaryPointer,
             typename SPHKernelFunction,
             typename DiffusiveTerm,
             typename ViscousTerm,
             typename EOS,
             typename BCType = typename SPHState::BCType,
             typename std::enable_if_t< std::is_same_v< BCType, WCSPH_BCTypes::MDBC >, bool > Enabled = true >
   void
   updateSolidBoundary( FluidPointer& fluid, BoudaryPointer& boundary, SPHState& sphState );

   /**
    * Function to realize boundary conditions for solid wall.
    * Realized by Generalized Wall Boundary Conditions (GWBC) - Adami, Hu 2012
    *
    * TODO: Not implemented.
    */
   template< typename FluidPointer,
             typename BoudaryPointer,
             typename SPHKernelFunction,
             typename DiffusiveTerm,
             typename ViscousTerm,
             typename EOS,
             typename BCType = typename SPHState::BCType,
             typename std::enable_if_t< std::is_same_v< BCType, WCSPH_BCTypes::GWBC >, bool > Enabled = true >
   void
   updateSolidBoundary( FluidPointer& fluid, BoudaryPointer& boundary, SPHState& sphState );


   //TODO: Where should be this placed
   template< typename FluidPointer,
             typename OpenBoudaryPointer,
             typename SPHKernelFunction,
             typename DiffusiveTerm,
             typename ViscousTerm,
             typename EOS >
   void
   interactionWithOpenBoundary( FluidPointer& fluid, OpenBoudaryPointer& openBoundary, SPHState& sphState );

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

   //TODO: Figure this out:
   template< typename FluidPointer,
             typename OpenBoudaryPointer,
             typename SPHKernelFunction,
             typename EOS >
   void
   extrapolateOpenBoundaryData( FluidPointer& fluid,
                                OpenBoudaryPointer& openBoundary,
                                SPHState& sphState );

};

} // SPH
} // ParticleSystem
} // TNL

#include "Interactions.hpp"

#include "BoundaryConditions/DBC.h"
//#include "BoundaryConditions/MDBC.h"

//#include "BoundaryConditions/DBC_with_enable_if_t.h"
//#include "BoundaryConditions/MDBC_with_enable_if_t.h"


#include "OpenBoundaryContirionsInteractions.h"

