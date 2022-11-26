#pragma once

#include <TNL/Containers/StaticVector.h>
#include <TNL/Containers/Array.h>
#include <TNL/Meshes/Grid.h>

namespace TNL {
namespace ParticleSystem {

template< typename ParticleConfig, typename Device = Devices::Cuda >
class ParticlesTraits
{
  public:
  static constexpr int spaceDimension = ParticleConfig::spaceDimension;

  using DeviceType = Device;
  using GlobalIndexType = typename ParticleConfig::GlobalIndexType;
  using LocalIndexType = typename ParticleConfig::LocalIndexType;
  using CellIndexType = typename ParticleConfig::CellIndexType;
  using RealType = typename ParticleConfig::RealType;

  /* particle related */
  using PointType = Containers::StaticVector< spaceDimension, RealType >;
  using PointArrayType = Containers::Array< PointType, Device, GlobalIndexType >;

  using NeighborListType = typename ParticleConfig::NeighborListType;
  using NeighborsCountArrayType = Containers::Array< LocalIndexType, Device, GlobalIndexType >;
  using NeighborsArrayType = Containers::Array< GlobalIndexType, Device, GlobalIndexType >;

  /* grid related */
  using GridType = Meshes::Grid< 2, RealType, DeviceType, GlobalIndexType >;
  using GridPointer = Pointers::SharedPointer< GridType >;
  using CellIndexArrayType = Containers::Array< CellIndexType, Device, GlobalIndexType >;

};

} // namepsace ParticleSystem
} // namespace TNL

