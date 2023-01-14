#pragma once

#include "SPHTraits.h"
namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename ParticleSystem, typename NeighborSearch, typename SPHCaseConfig, typename Variables >
class OpenBoundaryBuffer_orthogonal
{
   public:
   using DeviceType = typename ParticleSystem::Device;
   using ParticlePointer = typename Pointers::SharedPointer< ParticleSystem, DeviceType >;
   using NeighborSearchPointer = typename Pointers::SharedPointer< NeighborSearch, DeviceType >;
   using VariablesPointer = typename Pointers::SharedPointer< Variables, DeviceType >;

   using SPHTraitsType = SPHFluidTraits< SPHCaseConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   OpenBoundaryBuffer_orthogonal( GlobalIndexType size, GlobalIndexType sizeAllocated, RealType h, GlobalIndexType numberOfCells )
   : particles( size, sizeAllocated, h ), neighborSearch( particles, numberOfCells ), variables( sizeAllocated ) {};

   VariablesPointer&
   getOpenBoundaryVariables()
   {
      return this->variables;
   }

   const VariablesPointer&
   getOpenBoundaryVariables() const
   {
      return this->variables;
   }

   struct Parameters
   {
      VectorType orientation;
      VectorType velocity;
      RealType density;

      RealType bufferEdge;
      VectorType bufferWidth;
   };
   Parameters parameters;

   ParticlePointer particles;
   NeighborSearchPointer neighborSearch;
   VariablesPointer variables;
};

}
}
}

