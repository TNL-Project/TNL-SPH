#pragma once

#include "SPHTraits.h"
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#include <thrust/gather.h>

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

   using IndexArrayType = typename SPHTraitsType::IndexArrayType;
   using IndexArrayTypePointer = typename Pointers::SharedPointer< IndexArrayType, DeviceType >;

   OpenBoundaryBuffer_orthogonal( GlobalIndexType size, GlobalIndexType sizeAllocated, RealType h, GlobalIndexType numberOfCells )
   : particles( size, sizeAllocated, h ), neighborSearch( particles, numberOfCells ), variables( sizeAllocated ),
     particleMark( sizeAllocated ), sortPermutations( sizeAllocated ), points_swap( sizeAllocated ) {};

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

   void sortParticles()
   {
      GlobalIndexType numberOfParticle = particles->getNumberOfParticles();
      auto view_particleCellIndices = particles->getParticleCellIndices().getView();
      auto view_map = sortPermutations->getView();

      sortPermutations->forAllElements( [] __cuda_callable__ ( int i, int& value ) { value = i; } );
      //TODO: tempalte thrust device with TNL::Devices
      thrust::sort_by_key( thrust::device, view_particleCellIndices.getArrayData(), view_particleCellIndices.getArrayData() + numberOfParticle, view_map.getArrayData() );

      auto view_points = particles->getPoints().getView();
#ifdef PREFER_SPEED_OVER_MEMORY
      auto view_points_swap = points_swap.getView();
      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticle, view_points.getArrayData(), view_points_swap.getArrayData() );
      particles->getPoints().swap( points_swap );
#else
      //TODO: Error or implement.
#endif
      variables->sortVariables( sortPermutations, particles->getNumberOfParticles() );
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

   IndexArrayTypePointer sortPermutations;
   IndexArrayTypePointer particleMark;

#ifdef PREFER_SPEED_OVER_MEMORY
   using PointArrayType = typename ParticleSystem::PointArrayType;
   PointArrayType points_swap;
#endif
};

}
}
}

