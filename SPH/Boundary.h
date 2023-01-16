#pragma once

#include "SPHTraits.h"
#include <thrust/sort.h>
#include <thrust/execution_policy.h>

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename ParticleSystem, typename NeighborSearch, typename SPHCaseConfig, typename Variables, typename IntegratorVariables >
class Boundary
{
   public:
   using DeviceType = typename ParticleSystem::Device;
   using ParticlePointer = typename Pointers::SharedPointer< ParticleSystem, DeviceType >;
   using NeighborSearchPointer = typename Pointers::SharedPointer< NeighborSearch, DeviceType >;
   using VariablesPointer = typename Pointers::SharedPointer< Variables, DeviceType >;
   using IntegratorVariablesPointer = typename Pointers::SharedPointer< IntegratorVariables, DeviceType >;

   using SPHTraitsType = SPHFluidTraits< SPHCaseConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;

   using IndexArrayType = typename SPHTraitsType::IndexArrayType;

   Boundary( GlobalIndexType size, GlobalIndexType sizeAllocated, RealType h, GlobalIndexType numberOfCells )
   : particles( size, sizeAllocated, h ), neighborSearch( particles, numberOfCells ), variables( sizeAllocated ),
     integratorVariables( size ), sortPermutations( size ), points_swap( size ) {};

   void sortParticles()
   {
      GlobalIndexType numberOfParticle = particles->getNumberOfParticles();
      auto view_particleCellIndices = particles->getParticleCellIndices().getView();
      auto view_points = particles->getPoints().getView();
      auto view_sortPermutations = sortPermutations.getView();

      sortPermutations.forAllElements( [] __cuda_callable__ ( int i, int& value ) { value = i; } );
      //TODO: tempalte thrust device with TNL::Devices
      thrust::sort_by_key( thrust::device, view_particleCellIndices.getArrayData(), view_particleCellIndices.getArrayData() + numberOfParticle, view_sortPermutations.getArrayData() );
   }

   VariablesPointer&
   getBoundaryVariables()
   {
      return this->variables;
   }

   const VariablesPointer&
   getBoundaryVariables() const
   {
      return this->variables;
   }

   //void
   //sortParticles
   //{

   //}

   ParticlePointer particles;
   NeighborSearchPointer neighborSearch;
   VariablesPointer variables;
   IntegratorVariablesPointer integratorVariables;

   IndexArrayType sortPermutations;

#ifdef PREFER_SPEED_OVER_MEMORY
   using PointArrayType = typename ParticleSystem::PointArrayType;
   PointArrayType points_swap;
#endif
};

}
}
}

