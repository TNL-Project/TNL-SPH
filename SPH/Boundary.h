#pragma once

#include "PhysicalObject.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename ParticleSystem,
          typename NeighborSearch,
          typename SPHCaseConfig,
          typename Variables,
          typename IntegratorVariables >
class Boundary : public PhysicalObject< ParticleSystem, NeighborSearch, SPHCaseConfig, Variables, IntegratorVariables >
{
   public:
   using BaseType = PhysicalObject< ParticleSystem, NeighborSearch, SPHCaseConfig, Variables, IntegratorVariables >;
   using GlobalIndexType = typename BaseType::GlobalIndexType;
   using RealType = typename BaseType::RealType;
   using VariablesPointerType = typename BaseType::VariablesPointerType;

   Boundary( GlobalIndexType size, GlobalIndexType sizeAllocated, RealType h, GlobalIndexType numberOfCells )
   : PhysicalObject< ParticleSystem, NeighborSearch, SPHCaseConfig, Variables, IntegratorVariables >( size, sizeAllocated, h, numberOfCells ) {};

   VariablesPointerType&
   getBoundaryVariables()
   {
      return this->variables;
   }

   const VariablesPointerType&
   getBoundaryVariables() const
   {
      return this->variables;
   }
};

}
}
}

