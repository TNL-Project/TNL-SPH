// SPDX-FileComment: This file is part of TNL-SPH (https://gitlab.com/tnl-project/tnl-sph)
// SPDX-License-Identifier: MIT

#pragma once

#include <thrust/sort.h>
#include <thrust/gather.h>

#include <TNL/Containers/Array.h>
#include <TNL/Containers/ArrayView.h>
#include <TNL/Pointers/SharedPointer.h>

#include <SPH/shared/thrustExecPolicySelector.h>
#include <SPH/Models/WCSPH_DBC/BoundaryConditionsTypes.h>
#include <SPH/Models/WCSPH_DBC/OpenBoundaryConfig.h>
#include <SPH/SPHTraits.h>

namespace TNL::SPH {

template< typename SPHConfig, typename ModelConfig >
class OpenBoundaryConditionsBuffers
{
public:

   using DeviceType = typename SPHConfig::DeviceType;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using PairIndexType = Containers::StaticVector< 2, GlobalIndexType >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   using OpenBoundaryConfig = DBCOpenBoundaryConfig< SPHConfig >;

   /**
    * New set of functions to realize open boundary conditions.
    */
   template< typename FluidPointer,
             typename OpenBoundaryPointer >
   void
   applyOpenBoundary( RealType dt,
                      FluidPointer& fluid,
                      OpenBoundaryPointer& openBoundary,
                      OpenBoundaryConfig& openBoundaryParams );

   template< typename FluidPointer,
             typename OpenBoundaryPointer >
   void
   applyInletBoundaryCondition( RealType dt,
                                FluidPointer& fluid,
                                OpenBoundaryPointer& openBoundary,
                                OpenBoundaryConfig& openBoundaryParams );

   template< typename FluidPointer,
             typename OpenBoundaryPointer >
   void
   applyOuletBoundaryCondition( RealType dt,
                                FluidPointer& fluid,
                                OpenBoundaryPointer& openBoundary,
                                OpenBoundaryConfig& openBoundaryParams );

   template< typename OpenBoundaryPointer >
   GlobalIndexType
   moveInletBufferParticles( RealType dt, OpenBoundaryPointer& openBoundary );

   template< typename OpenBoundaryPointer >
   GlobalIndexType
   moveOutletBufferParticles( RealType dt, OpenBoundaryPointer& openBoundary );

   template< typename OpenBoundaryPointer >
   void
   sortBufferParticlesByMark( OpenBoundaryPointer& openBoundary );

   template< typename FluidPointer, typename OpenBoundaryPointer >
   void
   convertBufferToFluid( FluidPointer& fluid,
                         OpenBoundaryPointer& openBoundary,
                         OpenBoundaryConfig& openBoundaryParams,
                         const GlobalIndexType numberOfRetyped );

   template< typename FluidPointer, typename OpenBoundaryPointer >
   GlobalIndexType
   getFluidParticlesEnteringOutlet( FluidPointer& fluid, OpenBoundaryPointer& openBoundary );

   template< typename FluidPointer, typename OpenBoundaryPointer >
   void
   convertFluidToBuffer( FluidPointer& fluid,
                         OpenBoundaryPointer& openBoundary,
                         const GlobalIndexType fluidToBufferCount );

   /**
    * Functions to realize periodic boundary conditions.
    */

   template< typename FluidPointer, typename PeriodicBoundaryPatch >
   void
   periodicityParticleTransfer( FluidPointer& fluid,
                                PeriodicBoundaryPatch& periodicPatch );

};

}  // namespace TNL::SPH

#include <SPH/Models/WCSPH_DBC/OpenBoundaryConditions.hpp>

