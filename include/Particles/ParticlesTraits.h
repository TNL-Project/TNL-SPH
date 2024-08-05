#pragma once

#include <TNL/Containers/StaticVector.h>
#include <TNL/Containers/Array.h>
#include <TNL/Meshes/Grid.h>

//NeighborListType
#include <TNL/Algorithms/Segments/CSR.h>
#include <TNL/Algorithms/Segments/Ellpack.h>

namespace TNL {
namespace ParticleSystem {

template< typename ParticleConfig, typename Device >
class ParticlesTraits
{
   public:
   static constexpr int spaceDimension = ParticleConfig::spaceDimension;

   /* Particles.h related */
   using DeviceType = Device;
   using GlobalIndexType = typename ParticleConfig::GlobalIndexType;
   using LocalIndexType = typename ParticleConfig::LocalIndexType;
   using RealType = typename ParticleConfig::RealType;

   using PointType = Containers::StaticVector< spaceDimension, RealType >;
   using PointArrayType = Containers::Array< PointType, Device, GlobalIndexType >;

   /* ParticlesLinkedList.h related */
   using IndexVectorType = Containers::StaticVector< spaceDimension, GlobalIndexType >;
   using CellIndexType = typename ParticleConfig::CellIndexType;
   using CellIndexArrayType = Containers::Array< CellIndexType, Device, GlobalIndexType >;
   using PairIndexType = Containers::StaticVector< 2, GlobalIndexType >;
   using PairIndexArrayType = Containers::Array< PairIndexType, DeviceType, GlobalIndexType >;

   ///* ParticlesLinkedListWithList.h related */
   using IndexArrayType = typename Containers::Array< GlobalIndexType, DeviceType, GlobalIndexType >;
   //using NeighborListType = typename ParticleConfig::NeighborListType;
   using NeighborListType = typename Algorithms::Segments::Ellpack< DeviceType, int >;
   using NeighborsCountArrayType = Containers::Array< LocalIndexType, Device, GlobalIndexType >;
   using NeighborsArrayType = Containers::Array< GlobalIndexType, Device, GlobalIndexType >;
};

} // namepsace ParticleSystem
} // namespace TNL

