#pragma once

#include "SPHTraits.h"
namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename ParticleSystem, typename NeighborSearch, typename SPHCaseConfig, typename Variables >
class Fluid
{
   public:
   using DeviceType = typename ParticleSystem::Device;
   using ParticlePointer = typename Pointers::SharedPointer< ParticleSystem, DeviceType >;
   using NeighborSearchPointer = typename Pointers::SharedPointer< NeighborSearch, DeviceType >;
   using VariablesPointer = typename Pointers::SharedPointer< Variables, DeviceType >;

   using SPHTraitsType = SPHFluidTraits< SPHCaseConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;

   Fluid( GlobalIndexType size, GlobalIndexType sizeAllocated, RealType h, GlobalIndexType numberOfCells )
   : particles( size, sizeAllocated, h ), neighborSearch( particles, numberOfCells ), variables( sizeAllocated ) {};

   VariablesPointer&
   getFluidVariables()
   {
      return this->variables;
   }

   const VariablesPointer&
   getFluidVariables() const
   {
      return this->variables;
   }

   ParticlePointer particles;
   NeighborSearchPointer neighborSearch;
   VariablesPointer variables;
   //Variables variables;
};

}
}
}

