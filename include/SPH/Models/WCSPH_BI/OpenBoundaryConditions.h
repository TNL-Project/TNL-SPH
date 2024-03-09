#pragma once

#include <TNL/Containers/Array.h>
#include <TNL/Containers/ArrayView.h>
#include <TNL/Pointers/SharedPointer.h>

#include "../../shared/thrustExecPolicySelector.h"
#include <thrust/sort.h>
#include <thrust/gather.h>
#include "BoundaryConditionsTypes.h"
#include "OpenBoundaryConfig.h"

#include "../../SPHTraits.h"

#ifdef HAVE_MPI
#include "../../shared/utils.h"
#endif

namespace TNL {
namespace SPH {

template< typename SPHConfig >
class OpenBoundaryConditionsBuffers
{
public:

   using DeviceType = typename SPHConfig::DeviceType;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using PairIndexType = Containers::StaticVector< 2, GlobalIndexType >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   using OpenBoundaryConfig = BIOpenBoundaryConfig< SPHConfig >;

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
   template< typename FluidPointer, typename OpenBoundaryPointer >
   void
   applyPeriodicBoundary( FluidPointer& fluid,
                          OpenBoundaryPointer& periodicBoundary1,
                          OpenBoundaryPointer& periodicBoundary2,
                          OpenBoundaryConfig& periodicBoundary1Params,
                          OpenBoundaryConfig& periodicBoundary2Params );

   template< typename BoundaryPointer, typename OpenBoundaryPointer >
   void
   applyPeriodicBoundaryOnBoundary( BoundaryPointer& boundary,
                                    OpenBoundaryPointer& periodicBoundary1,
                                    OpenBoundaryPointer& periodicBoundary2,
                                    OpenBoundaryConfig& periodicBoundary1Params,
                                    OpenBoundaryConfig& periodicBoundary2Params );

   template< typename FluidPointer, typename OpenBoundaryPointer >
   void
   copyGhostParticles( FluidPointer& fluid,
                       OpenBoundaryPointer& sendingBuffer,
                       OpenBoundaryPointer& receivingBuffer,
                       VectorType shift );

   template< typename BoundaryPointer, typename OpenBoundaryPointer >
   void
   copyBoundaryGhostParticles( BoundaryPointer& boundary,
                               OpenBoundaryPointer& sendingBuffer,
                               OpenBoundaryPointer& receivingBuffer,
                               VectorType shift );

   template< typename FluidPointer, typename OpenBoundaryPointer >
   void
   periodicityParticleTransfer( FluidPointer& fluid,
                                OpenBoundaryPointer& periodicBuffer,
                                OpenBoundaryConfig& periodicBoundaryParams );

};

} // SPH
} // TNL

#include "OpenBoundaryConditions.hpp"

